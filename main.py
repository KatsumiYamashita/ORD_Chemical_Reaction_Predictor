import numpy as np
import openai
import pandas as pd
import re
import warnings


from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

# トレーニングデータを抽出する関数を定義する
def extract_training_data(nd_tnmt_A,
                          nd_tnmt_B,
                          df_smiles_mol_maccsfps,
                          td_number
                         ):

    # データセットからsmiles部分だけ取り出す
    df_smiles = df_smiles_mol_maccsfps.iloc[:, 0:3]
    # 計算したタニモト係数をdf_smilesに合わせる
    sr_tnmtA = pd.Series(nd_tnmt_A, name="tnmt_A")
    sr_tnmtB = pd.Series(nd_tnmt_B, name="tnmt_B")
    df_smiles_tnmt = pd.concat([df_smiles, sr_tnmtA, sr_tnmtB], axis=1)

    str_training_dataset= ""
    half_number = int(td_number/2)

    # データフレームを化合物Aのタニモト係数の降順で並び変える
    df_smiles_tnmt_A = df_smiles_tnmt.sort_values("tnmt_A", ascending=False)
    # ソートしたdfから化合物Aに対するタニモト係数上位トレーニングデータ数の半数を抜き取る
    df_training_data_A = df_smiles_tnmt_A.iloc[:half_number, 0:3]
    # 抜き取ったdfからstr型のトレーニングデータを作る
    for _, row in df_training_data_A.iterrows():
        template = "A: " + row['A'] + "\\" + "\n" + "B: " + row['B'] + "\\" + "\n" + "Y: " + row['Y'] + "\\" + "\n" + "\\" + "\n"
        str_training_dataset += template
    
    # データフレームを化合物Bのタニモト係数の降順で並び変える
    df_smiles_tnmt_B = df_smiles_tnmt.sort_values("tnmt_B", ascending=False)
    # ソートしたdfから化合物Aに対するタニモト係数上位トレーニングデータ数の半数を抜き取る
    df_training_data_B = df_smiles_tnmt_B.iloc[:half_number, 0:3]
    # 抜き取ったdfからstr型のトレーニングデータを作る
    for _, row in df_training_data_B.iterrows():
        template = "A: " + row['A'] + "\\" + "\n" + "B: " + row['B'] + "\\" + "\n" + "Y: " + row['Y'] + "\\" + "\n" + "\\" + "\n"
        str_training_dataset += template
    
    df_training_dataset = pd.concat([df_training_data_A, 
                                     df_training_data_B], 
                                     ignore_index=True
                                    )
    
    return str_training_dataset, df_training_dataset

def get_prodY_SMILES(test_A_smiles,\
                     test_B_smiles,\
                     training_dataset):

    question =\
    f"Synthesize compound 'Y' generated from the corresponding compound 'A' and 'B' below.\
    \
    {training_dataset}\
    \
    A: {test_A_smiles}\
    B: {test_B_smiles}\
    Y:\
    "

    response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-16k",
    messages=[
        {"role": "system", "content": "You are an chemist."},
        {"role": "user", "content": f"{question}"}
        ],
        max_tokens=300,
        temperature=0.5,
        )
    
    # gptのresponseからYの候補が含まれている部分を抜き出す
    y = response["choices"][0]["message"]["content"]
    
    return y