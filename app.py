# モジュールのインポート
import main
import numpy as np
import openai
import pubchempy as pcp
import pandas as pd
import pickle
import streamlit as st
import warnings

from streamlit_ketcher import st_ketcher
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# 化合物A,Bを入力する関数を定義する
def enter_reactants():

    # デフォルトで表示させておく化合物定数(SMILES)を定義する
    DEFAULT_A = r"CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
    DEFAULT_B = r"CONC(=O)C1=CC=CC=C1N"

    # pcpで取得する情報リストを定義する
    properties = ['IUPACName',  
                  'MolecularFormula',
                  'MolecularWeight',
                  'XLogP',
                  'TPSA',
                  'CanonicalSMILES']

    # streamlit appの表示を２分割するためのカラムを定義する
    col_A, col_B = st.columns(2, gap="medium")

    # 各カラムに入力画面を表示させる
    with col_A:
        # 入力部分を作成する
        entered_A = st.text_input("Enter reactant A 'SMILES'",
                                   DEFAULT_A)
        reactant_A_smiles = st_ketcher(entered_A,
                                       height = 400)

        # 入力された化合物情報をpubchemから入手する
        # データベースに情報がない化合物の場合の処理を後で考えておく
        df_reactant_A_pcp = pcp.get_properties(properties,
                                               reactant_A_smiles,
                                               'smiles',
                                               as_dataframe=True)
        df_reactant_A_pcp_transposed =\
              df_reactant_A_pcp.transpose() #.reset_index()
        # 情報データフレームを表示させる
        st.write("'Reactant A' Info from PubChem:")
        st.table(df_reactant_A_pcp_transposed) #dfを表示させる方法はいくつかあるみたい

    with col_B:
        # 入力部分を作成する
        entered_B = st.text_input("Enter reactant B 'SMILES'",
                                   DEFAULT_B)
        reactant_B_smiles = st_ketcher(entered_B,
                                       height = 400)

        # 入力された化合物情報をpubchemから入手する
        # データベースに情報がない化合物の場合の処理を後で考えておく
        df_reactant_B_pcp = pcp.get_properties(properties,
                                               reactant_B_smiles,
                                               'smiles',
                                               as_dataframe=True)
        df_reactant_B_pcp_transposed =\
              df_reactant_B_pcp.transpose() #.reset_index()
        # 情報データフレームを表示させる
        st.write("'Reactant B' Info from PubChem:")
        st.table(df_reactant_B_pcp_transposed)


    return reactant_A_smiles,\
            reactant_B_smiles,\

# 反応式を表示させる関数を定義する
def show_rxn_formula(smiles_A,
                     smiles_B,
                     smiles_Y):

    rxn_smiles = f"{smiles_A}.{smiles_B}>>{smiles_Y}"
    st_ketcher(rxn_smiles)

    return show_rxn_formula

# データセットを読み込む関数を定義する
@st.cache_data    #(hash_funcs={pandas.core.frame.DataFrame: my_hash_func}) #　２回目以降キャッシュから取り出す
def load_data(path):
    with open(path,'rb') as file:
        df_smiles_mol_maccsfps = pickle.load(file)
    nd_Amaccs = df_smiles_mol_maccsfps.iloc[:, 4].values
    nd_Bmaccs = df_smiles_mol_maccsfps.iloc[:, 5].values

    return df_smiles_mol_maccsfps,\
            nd_Amaccs,\
            nd_Bmaccs

def app_info():
    st.markdown(f"""
	# React ABY
	version {__version__}

	Prediction chemical product 'Y' from reactant 'A' and 'B' by GPT-3.5.
	""")

    st.write("Made by [Katsumi Yamashita](https://katsumiyamashita.github.io/).", unsafe_allow_html=True)

__version__ = "0.0.0"
app_name = "React_ABY"

# タブに表示させるタイトルを作成する
st.set_page_config(page_title=f'{app_name} {__version__}',
                   page_icon="⚗️",
                   initial_sidebar_state = "collapsed",
                   layout="wide")

# アプリケーションのセッション状態を保持
# 途中
ss = st.session_state
if 'debug' not in ss: ss['debug'] = {}

# アプリケーションタイトル
st.markdown("# React: A + B → Y")
st.sidebar.header("Report")
st.write(
    """report"""
)

# データセットパスを定義する
path = './ord_datasets/ord-data/df_SmilesMACCSFps.pickle'
# データセットをロードする
df_smiles_mol_maccsfps,\
nd_Amaccs,\
nd_Bmaccs\
= load_data(path)

# api_keyの入力
# 無くす予定
api_key = st.text_input("API keyを入力してください", \
                        type='password')
openai.api_key = api_key

# テスト化合物ABを変数に格納(関数の呼び出し)
reactant_A_smiles,\
reactant_B_smiles,\
= enter_reactants()

# テスト化合物ABのmaccs fpsを生成
reactant_A_mol = Chem.MolFromSmiles(reactant_A_smiles)
reactant_A_maccsfps = AllChem.GetMACCSKeysFingerprint(reactant_A_mol)

reactant_B_mol = Chem.MolFromSmiles(reactant_B_smiles)
reactant_B_maccsfps = AllChem.GetMACCSKeysFingerprint(reactant_B_mol)

# テスト分子とのTANIMOTO係数を計算する関数を定義する
def tnmt_similarity(nd_Amaccs,
                    nd_Bmaccs):

    nd_TNMT_A = DataStructs.TanimotoSimilarity(nd_Amaccs, 
                                               reactant_A_maccsfps)
    nd_TNMT_B = DataStructs.TanimotoSimilarity(nd_Bmaccs, 
                                               reactant_B_maccsfps)

    return nd_TNMT_A, nd_TNMT_B

# テスト分子とのTANIMOTO係数を計算する
uf_TNMTSimilarity = np.frompyfunc(tnmt_similarity, 2, 2)
nd_tnmt_A, nd_tnmt_B = uf_TNMTSimilarity(nd_Amaccs, nd_Bmaccs)

# トレーニングデータの取得
# トレーニングデータ数は選択できるようにしたい
str_training_dataset,\
df_training_dataset,\
= main.extract_training_data(nd_tnmt_A,
                             nd_tnmt_B,
                             df_smiles_mol_maccsfps,
                             10,
                             )

# 化合物ABを反応させる (Yを予測させる)
if api_key:

    df_Y =\
    main.get_prodY_SMILES(reactant_A_smiles,
                          reactant_B_smiles,
                          str_training_dataset)
    
    best_Y = df_Y.iloc[0,0]
  
    show_rxn_formula(reactant_A_smiles,
                     reactant_B_smiles,
                     best_Y) 
    
    st.dataframe(df_Y)

else:
    pass




