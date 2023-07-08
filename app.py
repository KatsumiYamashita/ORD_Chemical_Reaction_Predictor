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

PATH = './ord_datasets/ord_datasets_csv/df_SmilesMACCSFpsID.pickle'

PROPERTIES = ['IUPACName',
              'MolecularFormula',
              'MolecularWeight',
              'XLogP',
              'TPSA',
              'CanonicalSMILES']

# Application tub
st.set_page_config(page_title=f'{app_name} {__version__}',
                   page_icon="⚗️",
                   initial_sidebar_state = "collapsed",
                   layout="wide")

# session state
ss = st.session_state

ss.openai_api_key = OPENAI_API_KEY

ss.path = PATH

ss.properties = PROPERTIES

if "reactant_A" not in ss:
    ss.reactant_A = DEFAULT_A

if "reactant_B" not in ss:
    ss.reactant_B = DEFAULT_B


# データセットを読み込む関数を定義する
@st.cache_data   
def load_data(path):
    with open(path,'rb') as file:
        df_smiles_maccsfps_id = pickle.load(file)

    return df_smiles_maccsfps_id

# 化合物A,Bを入力する関数を定義する
def enter_reactant(DEFAULTCOMPOUND):

    entered_compound = st.text_input("Enter reactant A 'SMILES'",
                                    DEFAULTCOMPOUND)
    
    reactant = st_ketcher(entered_compound,
                            height = 400)
    
    reactant_mol = Chem.MolFromSmiles(reactant)
    reactant_smiles = Chem.MolToSmiles(reactant_mol,
                                       isomericSmiles=False,
                                       kekuleSmiles=False,
                                       allBondsExplicit=False,
                                       allHsExplicit=False,
                                       canonical=True)
    
    return reactant_smiles

# pcpで取得する情報リストを検索する関数を定義する
def get_info_reactants(reactant_smiles):

    df_reactant_pcp = pcp.get_properties(ss.properties,
                                         reactant_smiles,
                                         'smiles',
                                         as_dataframe=True)
    
    df_reactant_pcp_tr = df_reactant_pcp.transpose() 

    return df_reactant_pcp_tr

# トレーニングデータを描画する関数を定義する
def draw_rxn(rxn_smiles):

    drawer = rdMolDraw2D.MolDraw2DSVG(660,200)

    rxn = AllChem.ReactionFromSmarts(rxn_smiles, useSmiles=True)

    drawer.DrawReaction(rxn)

    drawer.FinishDrawing()

    svg_rxn = drawer.GetDrawingText() 

    return svg_rxn

# 結果を表示させる関数を定義する
def show_report(smiles_A,
                smiles_B,
                df_Y,
                df_training_dataset):
    
    t1,t2 = st.tabs(['Prediction by GPT-3.5','Training Data from the ORD'])

    with t1:
        st.dataframe(df_Y)
        
        y_candidates = list(df_Y["Y_candidates"])
        
        li_rxn_smiles = [f"{smiles_A}.{smiles_B}>>{i}" for i in y_candidates]
        
        st_ketcher(li_rxn_smiles[0])

    with t2:
    
        li_rxn_smiles = [f"{df_training_dataset.loc[i, 'A']}.\
                           {df_training_dataset.loc[i, 'B']}>>\
                           {df_training_dataset.loc[i, 'Y']}"\
                            for i in range(len(df_training_dataset))]

        li_svg_rxn = [draw_rxn(rxn_smiles) for rxn_smiles in li_rxn_smiles]

        #li_id = list(df_training_dataset["ID"])

        st.write("## Training Data Reaction from the Open Reaction Databese")

        for i in range(len(li_svg_rxn)):
    
            li_id = df_training_dataset.loc[i, "ID"].split(', ')
            
            st.write("------------------------------------------------------")
                         
            for j in li_id:
                
                st.write(f"https://open-reaction-database.org/client/id/{j}")

            st.image(li_svg_rxn[i], use_column_width=False)
            
            st.write("------------------------------------------------------")

    return show_report


def app_info():
    st.markdown(f"""
	# React ABY
	version {__version__}

	Prediction chemical product 'Y' from reactant 'A' and 'B' by GPT-3.5.
	""")

    st.write("Made by [Katsumi Yamashita](https://katsumiyamashita.github.io/).", unsafe_allow_html=True)

# Application title and description
st.markdown("# React: A + B → Y")

st.sidebar.header("Report")

st.write(
    """report"""
)

# データセットをロードする
df_smiles_maccsfps_id = load_data(ss.path)
nd_Amaccs = df_smiles_maccsfps_id.loc[:, "maccs_A"].values
nd_Bmaccs = df_smiles_maccsfps_id.loc[:, "maccs_B"].values

# OpenAI API Keyの認証を行なう
openai.api_key = ss.openai_api_key

# streamlit appの表示を２分割するためのカラムを定義する
col_A, col_B = st.columns(2, gap="medium")

with col_A:

    reactant_A_smiles = enter_reactant(ss.reactant_A)
    ss.reactant_A = reactant_A_smiles
        
    df_reactant_A_pcp_tr = get_info_reactants(reactant_A_smiles)
    ss.df_A_pcp = df_reactant_A_pcp_tr

    st.write("'Reactant A' Info from PubChem:")
    st.table(ss.df_A_pcp)

with col_B:
        
    reactant_B_smiles = enter_reactant(ss.reactant_B)
    ss.reactant_B = reactant_B_smiles
        
    df_reactant_B_pcp_tr = get_info_reactants(reactant_B_smiles)
    ss.df_B_pcp = df_reactant_B_pcp_tr

    st.write("'Reactant B' Info from PubChem:")
    st.table(ss.df_B_pcp)

# テスト化合物ABのmaccs fpsを生成
reactant_A_mol = Chem.MolFromSmiles(ss.reactant_A)
reactant_A_maccsfps = AllChem.GetMACCSKeysFingerprint(reactant_A_mol)
ss.reactant_A_maccsfps = reactant_A_maccsfps

reactant_B_mol = Chem.MolFromSmiles(ss.reactant_B)
reactant_B_maccsfps = AllChem.GetMACCSKeysFingerprint(reactant_B_mol)
ss.reactant_B_maccsfps = reactant_B_maccsfps

# テスト分子とのTANIMOTO係数を計算する関数を定義する
# 関数を使わないで内包表記でいいんじゃないか
def tnmt_similarity(nd_Amaccs,
                    nd_Bmaccs):

    nd_TNMT_A = DataStructs.TanimotoSimilarity(nd_Amaccs, 
                                               ss.reactant_A_maccsfps)
    nd_TNMT_B = DataStructs.TanimotoSimilarity(nd_Bmaccs, 
                                               ss.reactant_B_maccsfps)

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
                             df_smiles_maccsfps_id,
                             20,
                             )

predict_button = st.button("Pretict !", key=1)

# 化合物ABを反応させる (Yを予測させる)
if predict_button:
    
    response, df_Y =\
    main.get_prodY_SMILES(ss.reactant_A,
                          ss.reactant_B,
                          str_training_dataset)
    
    #st.dataframe(df_Y)
    #st.dataframe(response)
    #ss.best_Y = df_Y 
            
    show_report(ss.reactant_A,
                ss.reactant_B,
                df_Y,
                df_training_dataset)

    #except:

        #st.markdown("### Sorry...")
        #st.markdown("### Possible 'Y' could not be predicted.")  
        #st.markdown("### 👈 Adjust the parameters and try again!!")
    
else:
    pass






