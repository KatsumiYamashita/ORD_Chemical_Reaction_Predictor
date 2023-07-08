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
                          df_smiles_maccsfps_id,
                          td_number
                          ):

    # データセットからsmiles部分だけ取り出す
    df_smiles_id = df_smiles_maccsfps_id.loc[:, ["A", "B", "Y", "ID"]]
    # 計算したタニモト係数をdf_smilesに合わせる
    df_smiles_id["tnmt_A"] = pd.DataFrame(nd_tnmt_A)
    df_smiles_id["tnmt_B"] = pd.DataFrame(nd_tnmt_B)
    df_smiles_tnmt_id = df_smiles_id

    str_training_dataset= ""
    half_number = int(td_number/2)

    # データフレームを化合物Aのタニモト係数の降順で並び変える
    df_smiles_tnmt_id_A = df_smiles_tnmt_id.sort_values("tnmt_A", ascending=False)
    # ソートしたdfから化合物Aに対するタニモト係数上位トレーニングデータ数の半数を抜き取る
    df_training_data_A = df_smiles_tnmt_id_A.iloc[:half_number, 0:4]
    # 抜き取ったdfからstr型のトレーニングデータを作る
    for _, row in df_training_data_A.iterrows():
        template = "A: " + row['A'] + "\\" + "\n" + "B: " + row['B'] + "\\" + "\n" + "Y: " + row['Y'] + "\\" + "\n" + "\\" + "\n"
        str_training_dataset += template
    
    # データフレームを化合物Bのタニモト係数の降順で並び変える
    df_smiles_tnmt_id_B = df_smiles_tnmt_id.sort_values("tnmt_B", ascending=False)
    # ソートしたdfから化合物Aに対するタニモト係数上位トレーニングデータ数の半数を抜き取る
    df_training_data_B = df_smiles_tnmt_id_B.iloc[:half_number, 0:4]
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
    Answer 5 different each other candidates of 'y1 to y5' in this format, 'Y1:y1, Y2:y2, Y3:y3, Y4:y4, Y5:y5'.\
    If {test_A_smiles} and {test_B_smiles} are in {training_dataset}, the corresponding 'Y' should be included in y1-y5.\
    Compound Y must have a chemically flawless structure.\
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
        {"role": "system", "content": "Synthesize compound Y from test data compounds A and B."},
        {"role": "user", "content": f"{question}"},
        ],
        max_tokens=1000,
        temperature=0.25,
        )

    # gptのresponseからYの候補が含まれている部分を抜き出す
    response_text = response["choices"][0]["message"]["content"]
    #テンプレートの作成
    pattern = r":\s(.+)" 
    #テンプレートをもとにSMILESを抜き出す→リスト型
    product_Y_candidates = re.findall(pattern, response_text)
    #リスト型→データフレーム
    df_product_Y_candidates = pd.DataFrame({"Y_candidates":product_Y_candidates})
    #Molオブジェクトの生成
    df_product_Y_candidates["Y_candidates_mol"] =\
        df_product_Y_candidates["Y_candidates"].\
        apply(lambda smiles: Chem.MolFromSmiles(smiles))
    #maccs_fpsの生成
    df_product_Y_candidates["Y_candidates_maccs_fps"] =\
        df_product_Y_candidates["Y_candidates_mol"].\
        apply(lambda mol: AllChem.GetMACCSKeysFingerprint(mol)\
        if mol is not None else None
        )
    y = response["choices"][0]["message"]["content"]

    test_A_mol = Chem.MolFromSmiles(test_A_smiles)
    test_A_maccs_fps = AllChem.GetMACCSKeysFingerprint(test_A_mol)

    #テスト化合物Aを基準にタニモト係数を計算
    df_product_Y_candidates["Y_candidates_tnmt"] = \
        df_product_Y_candidates["Y_candidates_maccs_fps"].\
        apply(lambda maccs_fps: DataStructs.TanimotoSimilarity(test_A_maccs_fps, maccs_fps)\
        if maccs_fps is not None else None
        )

    #計算したタニモト係数を降順に並べる
    df_Y1 = df_product_Y_candidates.sort_values("Y_candidates_tnmt", ascending=False)
    df_Y2 = df_Y1.loc[:, ["Y_candidates", "Y_candidates_tnmt"]]

    #df_Y = df_product_Y_candidates.iloc[:, 0].reset_index()
    
    return response, df_Y2