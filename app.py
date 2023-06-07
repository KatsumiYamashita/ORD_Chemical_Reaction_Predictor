# モジュールのインポート
import streamlit as st
import main

from rdkit import Chem
from rdkit.Chem import Draw, rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
from context import training_dataset

# アプリケーションタイトル
st.title("React: A + B → Y")

# 化合物Aを入力する関数の定義
def enter_reactant_A():

    reactant_A = st.text_input("Enter reactant A", key="reactant_A")

    if reactant_A:
        st.write(f"reactant A is :{reactant_A}")
    else:
        st.write("Please enter reactant A.")
    return reactant_A

# 化合物Bを入力する関数の定義
def enter_reactant_B():

    reactant_B = st.text_input("Enter reactant B", key="reactant_B")

    if reactant_B:
        st.write(f"reactant B is :{reactant_B}")
    else:
        st.write("Please enter reactant B.")
    return reactant_B

# 化合物ABを変数に格納
entered_reactant_A = enter_reactant_A()
entered_reactant_B = enter_reactant_B()

# 反応物Yを探索する関数を定義
def react_AB():

    if entered_reactant_A and entered_reactant_B:
        generated_product_Y = main.get_prodY_SMILES(entered_reactant_A, \
                                        entered_reactant_B, \
                                        training_dataset)
        st.write(generated_product_Y)
        #product_Y = gpt のモデルによって回答の形式が異なる. SMILESだけ取得されるわけではないから難しい
    else:
        st.write("Waiting...")
    return generated_product_Y

product_Y = react_AB()

# 化合物ABYをSMILESからMolオブジェクトへ変換
#mol_A = Chem.MolFromSmiles(reactant_A2)
#mol_B = Chem.MolFromSmiles(reactant_B2)
#mol_Y = Chem.MolFromSmiles(product_Y2)

# molオブジェクトからimageを生成
def show_rxn_formula(A, B, Y):

    drawer = rdMolDraw2D.MolDraw2DSVG(660,350)
    rxn = Reactions.ReactionFromSmarts(f'{A}.{B}>>{Y}',\
                                        useSmiles=True)
    drawer.DrawReaction(rxn)
    drawer.FinishDrawing()
    svg_rxn = drawer.GetDrawingText()
    st.write("A + B → Y")
    st.image(svg_rxn, use_column_width=True)

    return show_rxn_formula

if entered_reactant_A and entered_reactant_B and product_Y:
    show_rxn_formula(entered_reactant_A, entered_reactant_B, product_Y)
   






