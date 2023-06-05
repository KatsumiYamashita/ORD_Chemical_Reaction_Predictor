# モジュールのインポート
import streamlit as st
import main

from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

# アプリケーションタイトル
st.title("React: A + B → Y")

# 化合物Aを入力する関数の定義
def reactant_A():

    reactant_A = st.text_input("Enter reactant A")

    if reactant_A:
        pass #st.write(f"reactant A is :{reactant_A}")
    else:
        st.write("Please enter reactant A.")
    return reactant_A

# 化合物Bを入力する関数の定義
def reactant_B():

    reactant_B = st.text_input("Enter reactant B")

    if reactant_B:
        pass #st.write(f"reactant B is :{reactant_B}")
    else:
        st.write("Please enter reactant B.")
    return reactant_B

# 化合物ABを変数に格納
reactant_A = reactant_A()
reactant_B = reactant_B()

# 化合物ABをSMILESからMolオブジェクトへ変換
mol_A = Chem.MolFromSmiles(reactant_A)
mol_B = Chem.MolFromSmiles(reactant_B)

# molオブジェクトからimageを生成
img_A = Draw.MolToImage(mol_A)
img_B = Draw.MolToImage(mol_B)
st.image(img_A, use_column_width=True)
st.image(img_B, use_column_width=True)

def react_ab():

    if reactant_A and reactant_B:
        product_Y = main.react_ai(reactant_A, reactant_B)
        st.write(product_Y)
    else:
        st.write("Waiting...")

react_ab()


