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
reactant_A2 = reactant_A()
reactant_B2 = reactant_B()

# 反応物Yを探索する関数を定義
def react_ab():

    if reactant_A2 and reactant_B2:
        answer = main.react_ai(reactant_A2, reactant_B2)
        st.write(answer)
        #product_Y = 
    else:
        st.write("Waiting...")

react_ab()
product_Y2 = react_ab()

# 化合物ABYをSMILESからMolオブジェクトへ変換
#mol_A = Chem.MolFromSmiles(reactant_A2)
#mol_B = Chem.MolFromSmiles(reactant_B2)
#mol_Y = Chem.MolFromSmiles(product_Y2)

# molオブジェクトからimageを生成
def A_B_Y():

    if reactant_A2 and reactant_B2 and product_Y2:
        drawer = rdMolDraw2D.MolDraw2DSVG(660,350)
        rxn = Reactions.ReactionFromSmarts(f'{reactant_A2}.{reactant_B2}>>{product_Y2}', useSmiles=True)
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        svg_rxn = drawer.GetDrawingText()
        st.write("A + B → Y")
        st.image(svg_rxn, use_column_width=True)
    else:
        pass
   
A_B_Y()





