# モジュールのインポート
import extract
import main
import openai
import pubchempy as pcp
import pandas as pd
import streamlit as st

from streamlit_ketcher import st_ketcher
from rdkit import rdBase, Chem, DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, Draw, rdMHFPFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols

# 化合物A,Bを入力する関数の定義
def enter_reactants():

    DEFAULT_A = r"CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
    DEFAULT_B = r"CONC(=O)C1=CC=CC=C1N"

    col_A, col_B = st.columns(2, gap="large")

    with col_A:
        entered_A = st.text_input("Enter reactant A 'SMILES'", DEFAULT_A)
        reactant_A_smiles = st_ketcher(entered_A)
        reactant_A_mol = Chem.MolFromSmiles(reactant_A_smiles)

    with col_B:
        entered_B = st.text_input("Enter reactant B 'SMILES'", DEFAULT_B)
        reactant_B_smiles = st_ketcher(entered_B)
        reactant_B_mol = Chem.MolFromSmiles(reactant_B_smiles)

    return reactant_A_smiles, reactant_B_smiles, reactant_A_mol, reactant_B_mol

# 反応式を表示させる関数を定義
def show_rxn_formula(A, B, Y):

    rxn = f"{A}.{B}>>{Y}"
    st_ketcher(rxn)

    return show_rxn_formula

# アプリケーションタイトル
st.set_page_config(layout="wide")
st.title("React: A + B → Y")

# データセットのインポート
ord_dataset_id = \
    "ord_dataset-00005539a1e04c809a9a78647bea649c"

df_smiles_ABY,\
df_mol_ABY,\
df_smiles_mol\
= extract.import_dataset(ord_dataset_id)

# df_molの各要素に対して maccs_fps を生成
df_maccs_fps_ABY,\
df_smiles_maccs_fps\
= extract.generate_maccs_fps_df(df_mol_ABY, df_smiles_ABY)

# api_keyの入力
api_key = st.text_input("API keyを入力してください", \
                        type='password')
openai.api_key = api_key

# テスト化合物ABを変数に格納(関数の呼び出し)
reactant_A_smiles,\
reactant_B_smiles,\
reactant_A_mol, \
reactant_B_mol \
= enter_reactants()

# テスト化合物ABのmaccs fpsを生成
reactant_A_maccs_fps,\
reactant_A_count_one \
= extract.generate_maccs_fps(reactant_A_mol)

reactant_B_maccs_fps,\
reactant_B_count_one \
= extract.generate_maccs_fps(reactant_B_mol)

# トレーニングデータの取得
training_dataset\
= extract.extract_training_data(df_smiles_maccs_fps,\
                                df_maccs_fps_ABY,\
                                reactant_A_maccs_fps,\
                                reactant_B_maccs_fps,\
                                35)

# 化合物ABを反応させる (Yを予測させる)
if api_key:

    col_A, col_B = st.columns(2, gap="large")
    df_Y =\
    main.get_prodY_SMILES(reactant_A_smiles, \
                            reactant_B_smiles, \
                            training_dataset)
    
    best_Y = df_Y.iloc[0, 1]

    with col_A:
        st.dataframe(df_Y)

    with col_B:
        show_rxn_formula(reactant_A_smiles, \
                        reactant_B_smiles, \
                        best_Y) 

else:
    pass

