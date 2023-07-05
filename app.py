# モジュールのインポート
import main
import numpy as np
import openai
import pubchempy as pcp
import pandas as pd
import pickle
import streamlit as st
import warnings

from config import OPENAI_API_KEY
from streamlit_ketcher import st_ketcher
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

__version__ = "0.0.0"
app_name = "React_ABY"

# CONSTANTS
DEFAULT_A = r"CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
DEFAULT_B = r"CONC(=O)C1=CC=CC=C1N"
PATH = './ord_datasets/ord_datasets_csv/df_SmilesMACCSFps.pickle'

# タブに表示させるタイトルを作成する
st.set_page_config(page_title=f'{app_name} {__version__}',
                   page_icon="⚗️",
                   initial_sidebar_state = "collapsed",
                   layout="wide")

# session state
ss = st.session_state

ss.openai_api_key = OPENAI_API_KEY
ss.path = PATH

if "entered_A" not in ss:
    ss.entered_A = DEFAULT_A

if "entered_B" not in ss:
    ss.entered_B = DEFAULT_B

# Application title and description
st.markdown("# React: A + B → Y")
st.sidebar.header("Report")
st.write(
    """report"""
)

# データセットを読み込む関数を定義する
@st.cache_data   
def load_data(path):
    with open(path,'rb') as file:
        df_smiles_maccsfps = pickle.load(file)
    nd_Amaccs = df_smiles_maccsfps.iloc[:, 4].values
    nd_Bmaccs = df_smiles_maccsfps.iloc[:, 5].values

    return df_smiles_maccsfps,\
            nd_Amaccs,\
            nd_Bmaccs

# 化合物A,Bを入力する関数を定義する
def enter_reactants():

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
                                   ss.entered_A)
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
                                   ss.entered_B)
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

# 結果を表示させる関数を定義する
def show_report(smiles_A,
                smiles_B,
                smiles_Y,
                df_training_dataset):
    
    t1,t2 = st.tabs(['Prediction by GPT-3.5','Training Data from The ORD'])

    with t1:
        rxn_smiles = f"{smiles_A}.{smiles_B}>>{smiles_Y}"
        st_ketcher(rxn_smiles)

    with t2:
        st.write("## Reference reaction from The Open Reaction Databese")

        df = df_training_dataset.copy()

        li_rxn_smiles = [f"{df.loc[i, 'A']}.{df.loc[i, 'B']}>>{df.loc[i, 'Y']}" for i in range(len(df))]

        df["rxn_smiles"] = li_rxn_smiles

        def draw_rxn(rxn_smiles):
            drawer = rdMolDraw2D.MolDraw2DSVG(660,200)
            rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            svg_rxn = drawer.GetDrawingText()  

            return svg_rxn
        
        li_svg_rxn = [draw_rxn(rxn_smiles) for rxn_smiles in li_rxn_smiles]

        for i in range(len(df)):
            st.image(li_svg_rxn[i], use_column_width=True)

    return show_report


def app_info():
    st.markdown(f"""
	# React ABY
	version {__version__}

	Prediction chemical product 'Y' from reactant 'A' and 'B' by GPT-3.5.
	""")

    st.write("Made by [Katsumi Yamashita](https://katsumiyamashita.github.io/).", unsafe_allow_html=True)

# データセットをロードする
df_smiles_mol_maccsfps,\
nd_Amaccs,\
nd_Bmaccs\
= load_data(ss.path)

openai.api_key = ss.openai_api_key

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
# 関数を使わないで内包表記でいいんじゃないか
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
                             6,
                             )

predict_button = st.button("Pretict !", key=1)

# 化合物ABを反応させる (Yを予測させる)
if predict_button:
    
    st.write("predict_button pushed")
    #try:
    y =\
    main.get_prodY_SMILES(reactant_A_smiles,
                         reactant_B_smiles,
                         str_training_dataset)
    best_Y = y 
            
    show_report(reactant_A_smiles,
                reactant_B_smiles,
                best_Y,
                df_training_dataset)

    #except:

        #st.markdown("### Sorry...")
        #st.markdown("### Possible 'Y' could not be predicted.")  
        #st.markdown("### 👈 Adjust the parameters and try again!!")
    
else:
    pass






