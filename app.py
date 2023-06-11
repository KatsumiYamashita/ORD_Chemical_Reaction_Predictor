# モジュールのインポート
import main
import openai
import pubchempy as pcp
import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from context import training_dataset

# アプリケーションタイトル
st.title("React: A + B → Y")

# 化合物Aを入力する関数の定義
def enter_reactant_A():

    type_A = st.selectbox("type reactant A", ['Name', 'SMILES'])
    reactant_A = None
    mol_A = None

    if type_A == "SMILES":
        entered_A = st.text_input("Enter reactant A", key="reactant_A") # AのSMILESを入力させる
        mol_A = Chem.MolFromSmiles(entered_A)
        if mol_A == None:
            st.write('<span style="color: red;">Error</span>', \
                     ": Please check entered 'reactant A' type.", \
                    unsafe_allow_html=True)
        else:
            reactant_A = entered_A
        
    elif type_A == "Name":
        entered_A = st.text_input("Enter reactant A", key="reactant_A") # AのIUPAC名を入力させる
        if entered_A:    
            entered_A_pcp_li = pcp.get_compounds(entered_A, 'name') # AのIUPAC名からPCPで化合物情報を取得
            if len(entered_A_pcp_li) == 1:
                entered_A_pcp = entered_A_pcp_li[0]
                reactant_A = entered_A_pcp.canonical_smiles   # AのSMILESを取得
                mol_A = Chem.MolFromSmiles(reactant_A) # molオブジェクトの生成
            else: 
                st.write('<span style="color: red;">Error</span>', \
                     ": Please check entered 'reactant A' type.", \
                    unsafe_allow_html=True)

    if mol_A:
        # コンテナの準備
        drawer = rdMolDraw2D.MolDraw2DSVG(100,100)
        tmol_A = rdMolDraw2D.PrepareMolForDrawing(mol_A)

        drawer.DrawMolecule(tmol_A)
        drawer.FinishDrawing()
        # 描画
        svg_A = drawer.GetDrawingText()
        st.image(svg_A, use_column_width=True)

    return reactant_A

# 化合物Bを入力する関数の定義
def enter_reactant_B():

    type_B = st.selectbox("type reactant B", ['Name', 'SMILES'])
    reactant_B = None
    mol_B = None

    if type_B == "SMILES":
        entered_B = st.text_input("Enter reactant B", key="reactant_B") # BのSMILESを入力させる
        mol_B = Chem.MolFromSmiles(entered_B)
        if mol_B == None:
            st.write('<span style="color: red;">Error</span>', \
                     ": Please check entered 'reactant B' type.", \
                    unsafe_allow_html=True)
        else:
            reactant_B = entered_B
        
    elif type_B == "Name":
        entered_B = st.text_input("Enter reactant B", key="reactant_B") # AのIUPAC名を入力させる
        if entered_B:    
            entered_B_pcp_li = pcp.get_compounds(entered_B, 'name') # AのIUPAC名からPCPで化合物情報を取得
            if len(entered_B_pcp_li) == 1:
                entered_B_pcp = entered_B_pcp_li[0]
                reactant_B = entered_B_pcp.canonical_smiles   # AのSMILESを取得
                mol_A = Chem.MolFromSmiles(reactant_B) # molオブジェクトの生成
            else: 
                st.write('<span style="color: red;">Error</span>', \
                     ": Please check entered 'reactant B' type.", \
                    unsafe_allow_html=True)

    if mol_B:
        # コンテナの準備
        drawer = rdMolDraw2D.MolDraw2DSVG(100,100)
        tmol_B = rdMolDraw2D.PrepareMolForDrawing(mol_B)

        drawer.DrawMolecule(tmol_B)
        drawer.FinishDrawing()
        # 描画
        svg_B = drawer.GetDrawingText()
        st.image(svg_B, use_column_width=True)

    return reactant_B

# 反応物Yを探索する関数を定義
def react_AB():

    if entered_reactant_A and entered_reactant_B:
        generated_product_Y = main.get_prodY_SMILES(entered_reactant_A,\
                                                    entered_reactant_B,\
                                                    training_dataset)
        st.write(generated_product_Y)
        #product_Y = gpt のモデルによって回答の形式が異なる. SMILESだけ取得されるわけではないから難しい
    else:
        st.write("Waiting...")
    return generated_product_Y

# molオブジェクトからimageを生成
def show_rxn_formula(A, B, Y):

    drawer = rdMolDraw2D.MolDraw2DSVG(660,300)

    try:
        rxn = Reactions.ReactionFromSmarts(f'{A}.{B}>>{Y}',\
                                        useSmiles=True)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        svg_rxn = drawer.GetDrawingText()
        st.write("A + B → Y")
        st.image(svg_rxn, use_column_width=True)

    except ValueError as e:
        st.write("Impossible product Y was generated")
        st.write(str(e)) 

    return show_rxn_formula

#api_keyの入力
api_key = st.text_input("API keyを入力してください", type='password')
openai.api_key = api_key

# 化合物ABを変数に格納(関数の呼び出し)
entered_reactant_A = enter_reactant_A()
entered_reactant_B = enter_reactant_B()

# 化合物ABを反応させる (Yを予測させる)
if entered_reactant_A and entered_reactant_B:
    generated_product_Y = react_AB()

# A, B, およびYの反応式を表示させる
if entered_reactant_A and entered_reactant_B and generated_product_Y:
    show_rxn_formula(entered_reactant_A, entered_reactant_B, generated_product_Y)