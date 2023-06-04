import streamlit as st
import main

st.title("React: A + B → Y")

def reactant_A():

    reactant_A: str = st.text_input("Enter reactant A")

    if reactant_A:
        st.write(f"reactant A is :{reactant_A}")
    else:
        st.write("Please enter reactant A.")
    return reactant_A

def reactant_B():

    reactant_B: str = st.text_input("Enter reactant B")

    if reactant_B:
        st.write(f"reactant B is :{reactant_B}")
    else:
        st.write("Please enter reactant B.")
    return reactant_B

reactant_A2 = reactant_A()
reactant_B2 = reactant_B()

def react_ab():

    if reactant_A2 and reactant_B2:
        product_Y = main.react_ai(reactant_A2, reactant_B2)
        st.write(product_Y)
    else:
        st.write("Waiting...")
    
    return react_ab()

react_ab()
