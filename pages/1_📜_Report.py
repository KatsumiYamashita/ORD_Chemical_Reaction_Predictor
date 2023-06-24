import extract
import main
import openai
import pubchempy as pcp
import pandas as pd
import streamlit as st

from app import reactant_A_smiles, reactant_B_smiles,df_training_dataset
from streamlit_ketcher import st_ketcher
from st_aggrid import AgGrid
from rdkit import rdBase, Chem, DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, Draw, rdMHFPFingerprint, PandasTools
from rdkit.Chem.Fingerprints import FingerprintMols

st.set_page_config(page_title="Report",
                   page_icon="📜",
                   layout="wide")

st.markdown("# Reaport")
st.sidebar.header("Report")
st.write(
    """report"""
)

df_training_dataset_structure = df_training_dataset.copy()
df_training_dataset_structure["A"] = \
PandasTools.AddMoleculeColumnToFrame(df_training_dataset, "A")
df_training_dataset_structure["B"] = \
PandasTools.AddMoleculeColumnToFrame(df_training_dataset, "B")
df_training_dataset_structure["Y"] = \
PandasTools.AddMoleculeColumnToFrame(df_training_dataset, "Y")

AgGrid(df_training_dataset_structure,
       fit_columns_on_grid_load=True,
       height=400)

properties = ['IUPACName',
              'MolecularFormula',
              'MolecularWeight',
              'XLogP',
              'TPSA',
              'CanonicalSMILES']

col_A, col_B = st.columns(2, gap="large")

with col_A:

    df_reactant_A_pcp = pcp.get_properties(properties,
                                           reactant_A_smiles,
                                           'smiles',
                                           as_dataframe=True)
    
    df_reactant_A_pcp = df_reactant_A_pcp.transpose()

    AgGrid(df_reactant_A_pcp,
           fit_columns_on_grid_load=True,
           height=400)

with col_B:

    df_reactant_B_pcp = pcp.get_properties(properties,
                                           reactant_B_smiles,
                                           'smiles',
                                           as_dataframe=True)

    df_reactant_B_pcp = df_reactant_B_pcp.transpose()

    AgGrid(df_reactant_B_pcp,
           fit_columns_on_grid_load=True,
           height=400)
    


