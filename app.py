import main
import numpy as np
import openai
import os
import pubchempy as pcp
import pandas as pd
import pickle
import streamlit as st
import warnings

from config import OPENAI_API_KEY
from streamlit_ketcher import st_ketcher
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

__version__ = "0.0.0"
app_name = "ORD Chemical Reaction Predictor"

# Application tub
st.set_page_config(
    page_title=f'{app_name} {__version__}',
    page_icon="⚗️",
    initial_sidebar_state = "collapsed",
    layout="wide"
)

# CONSTANTS
DEFAULT_A = r"Cc1nn(C)cc1Nc1cc(I)c(C(F)(F)F)cn1"
DEFAULT_B = r"CONC(=O)c1ccccc1N"
PATH = './ord_datasets/df_SmilesMACCSFpsID.pickle'
PROPERTIES = ['IUPACName', 'MolecularFormula', 'MolecularWeight', 'XLogP', 'TPSA', 'CanonicalSMILES']

# OPENAI_API_KEY
# os.environ['OPENAI_API_KEY'] = st.secrets['OPENAI_API_KEY']

# session state
ss = st.session_state

if "reactant_A" not in ss:
    ss.reactant_A = DEFAULT_A

if "reactant_B" not in ss:
    ss.reactant_B = DEFAULT_B

ss.path = PATH
ss.openai_api_key = OPENAI_API_KEY #os.environ['OPENAI_API_KEY']
ss.properties = PROPERTIES

def spacer(n=2, line=False, next_n=0):
	for _ in range(n):
		st.write('')
	if line:
		st.tabs([' '])
	for _ in range(next_n):
		st.write('')

def app_info():
    st.markdown(f"""
                
	# {app_name}
	version {__version__}
 
	""")
    spacer(1)
    st.write("Made by [Katsumi Yamashita](https://katsumiyamashita.github.io/)", unsafe_allow_html=True)
    spacer(1)
    st.write("Source code in [Github Repo](https://github.com/KatsumiYamashita/React_ABY)", unsafe_allow_html=True)
    spacer(1)
    st.write("-------------------------------")
    spacer(1)
    st.markdown("# References")
    spacer(1)
    st.write("・[The Open Reaction Database](https://docs.open-reaction-database.org/en/latest/)", unsafe_allow_html=True)
    st.write("・[streamlit-ketcher](https://github-com.translate.goog/mik-laj/streamlit-ketcher?ref=blog.streamlit.io&_x_tr_sl=en&_x_tr_tl=ja&_x_tr_hl=ja&_x_tr_pto=sc)", unsafe_allow_html=True)
    st.write("・[PubChem](https://pubchem.ncbi.nlm.nih.gov/)", unsafe_allow_html=True)
    st.write("・[OpenAI](https://platform.openai.com/overview)", unsafe_allow_html=True)
    spacer(1)
    st.write("-------------------------------")
    spacer(1)
    st.markdown("# Acknowledgment")
    spacer(1)
    st.write("[Suguru Tanaka](https://suguru-tanaka.com/)", unsafe_allow_html=True)
    st.write("I created my first St. App as part of his coaching program. I am truly grateful for his passionate guidance.")
    
    return app_info

def initialize_session_state():
    ss = st.session_state

    if "reactant_A" not in ss:
        ss.reactant_A = DEFAULT_A

    if "reactant_B" not in ss:
        ss.reactant_B = DEFAULT_B

    if "df_A_pcp" not in ss:
        ss.df_A_pcp = None

    if "df_B_pcp" not in ss:
        ss.df_B_pcp = None

# データセットを読み込む関数を定義する
@st.cache_data   
def load_data(path):
    with open(path,'rb') as file:
        df_smiles_maccsfps_id = pickle.load(file)
    return df_smiles_maccsfps_id

# 化合物A,Bを入力する関数を定義する
def enter_reactant(X, DEFAULTCOMPOUND):

    entered_compound = st.text_input(f"Reactant {X}: ", DEFAULTCOMPOUND)
    
    reactant = st_ketcher(entered_compound, height = 400)
    
    reactant_mol = Chem.MolFromSmiles(reactant)
    reactant_smiles = Chem.MolToSmiles(
        reactant_mol,
        isomericSmiles=False,
        kekuleSmiles=False,
        allBondsExplicit=False,
        allHsExplicit=False,
        canonical=True
    )
    return reactant_smiles

# pcpで取得する情報リストを検索する関数を定義する
def get_info_reactants(reactant_smiles):

    df_reactant_pcp = pcp.get_properties(
        ss.properties,
        reactant_smiles,
        'smiles',
        as_dataframe=True
    )
    
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
def show_report(smiles_A, smiles_B, df_Y, df_training_dataset):
    
    t1,t2 = st.tabs(['Prediction by GPT-3.5','Training Data from the ORD'])

    with t1:
        
        st.write("## Prediction product 'Y' from reactant 'A' and 'B' by GPT-3.5")
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

        st.write("## Training Data Reaction from 'the Open Reaction Databese'")
        st.write("Github pages → https://docs.open-reaction-database.org/en/latest/")
        st.write("ORD → https://open-reaction-database.org/client/browse")

        for i in range(len(li_svg_rxn)):
            li_id = df_training_dataset.loc[i, "ID"].split(', ')
            st.write("------------------------------------------------------")
            for j in li_id:
                st.write(f"https://open-reaction-database.org/client/id/{j}")
            st.image(li_svg_rxn[i], use_column_width=False)   
            st.write("------------------------------------------------------")

    return show_report


with st.sidebar:
    app_info()

# データセットをロードする
ss.df_smiles_maccsfps_id = load_data(ss.path)
ss.nd_Amaccs = ss.df_smiles_maccsfps_id.loc[:, "maccs_A"].values
ss.nd_Bmaccs = ss.df_smiles_maccsfps_id.loc[:, "maccs_B"].values

# OpenAI API Keyの認証を行なう
# openai.api_key = ss.openai_api_key

# Application title and description
st.title(f"{app_name}")

st.write("This App searches organic reaction data containing your input compounds from the dataset of over 300,000 entries in [the Open Reaction Database](https://open-reaction-database.org/client/browse).",
         "Even when the compound is not present in the dataset, app attempts products prediction using GPT-3.5 trained on reaction SMILES.")

spacer(2)

st.markdown("### 1.  Draw your compounds ! ")

# streamlit appの表示を２分割するためのカラムを定義する
col_A, col_B = st.columns(2, gap="medium")

with col_A:

    ss.reactant_A = enter_reactant("A", ss.reactant_A)
    ss.df_A_pcp = get_info_reactants(ss.reactant_A)
    cid_A = ss.df_A_pcp.columns
    ss.url_A =f"https://pubchem.ncbi.nlm.nih.gov/compound/{str(cid_A[0])}"
    st.write(f"'Reactant A' Info from [PubChem]({ss.url_A}):")
    st.table(ss.df_A_pcp)

with col_B:
        
    ss.reactant_B = enter_reactant("B", ss.reactant_B)
    ss.df_B_pcp = get_info_reactants(ss.reactant_B)
    cid_B = ss.df_B_pcp.columns
    ss.url_B =f"https://pubchem.ncbi.nlm.nih.gov/compound/{str(cid_B[0])}"
    st.write(f"'Reactant A' Info from [PubChem]({ss.url_B}):")
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
ss.nd_tnmt_A, ss.nd_tnmt_B = uf_TNMTSimilarity(ss.nd_Amaccs, ss.nd_Bmaccs)

# トレーニングデータの取得
# トレーニングデータ数は選択できるようにしたい
str_training_dataset,\
df_training_dataset,\
= main.extract_training_data(ss.nd_tnmt_A,
                             ss.nd_tnmt_B,
                             ss.df_smiles_maccsfps_id,
                             20,
                             )

spacer(2)

st.markdown("### 2.  Press the 'Prediction' button !! ")

predict_button = st.button("Prediction !", key=1)

# 化合物ABを反応させる (Yを予測させる)
if predict_button:
    try:
        response, df_Y =\
        main.get_prodY_SMILES(ss.reactant_A,
                                ss.reactant_B,
                                str_training_dataset) 
                
        show_report(ss.reactant_A,
                        ss.reactant_B,
                        df_Y,
                        df_training_dataset)

    except:

        st.markdown("### Sorry...")
        st.markdown("### Possible 'Y' could not be predicted.")  
        #st.markdown("### 👈 Adjust the parameters and try again!!")
    
else:
    pass






