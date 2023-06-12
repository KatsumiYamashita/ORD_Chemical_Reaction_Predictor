# モジュールのインポート
import main
import openai
import pubchempy as pcp
import streamlit as st

from context import training_dataset
from streamlit_ketcher import st_ketcher

# 化合物A,Bを入力する関数の定義
def enter_reactants():

    DEFAULT_A = r"CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
    entered_A = st.text_input("Enter reactant A 'SMILES'", DEFAULT_A)
    reactant_A = st_ketcher(entered_A)
    st.markdown(f"reactant A is ``{reactant_A}``")

    DEFAULT_B = r"CONC(=O)C1=CC=CC=C1N"
    entered_B = st.text_input("Enter reactant B 'SMILES'", DEFAULT_B)
    reactant_B = st_ketcher(entered_B)
    st.markdown(f"reactant B is: ``{reactant_B}``")

    return reactant_A, reactant_B

# 反応物Yを探索する関数を定義
def react_AB():

    if entered_reactant_A and entered_reactant_B:
        generated_product_Y = main.get_prodY_SMILES(entered_reactant_A,\
                                                    entered_reactant_B,\
                                                    training_dataset)
        
        #product_Y = gpt のモデルによって回答の形式が異なる. SMILESだけ取得されるわけではないから難しい
    else:
        st.write("Waiting...")
    return generated_product_Y

def show_rxn_formula(A, B, Y):

    rxn = f"{A}.{B}>>{Y}"
    st_ketcher(rxn)

    return show_rxn_formula

# アプリケーションタイトル
st.title("React: A + B → Y")

#api_keyの入力
api_key = st.text_input("API keyを入力してください", type='password')
openai.api_key = api_key

# 化合物ABを変数に格納(関数の呼び出し)
entered_reactant_A, entered_reactant_B = enter_reactants()

# 化合物ABを反応させる (Yを予測させる)
if api_key:
    generated_product_Y = react_AB()
else:
    pass

# A, B, およびYの反応式を表示させる
if entered_reactant_A and entered_reactant_B and generated_product_Y:
    show_rxn_formula(entered_reactant_A, entered_reactant_B, generated_product_Y)    
else:
    pass