import numpy as np
import openai
import pandas as pd
import re
import warnings

from rdkit import rdBase, Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

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

    str_training_dataset= \
    "This is the reaction (A + B → Y) training dataset :\n\
    "
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

    test = f"A: {test_A_smiles}\
             B: {test_B_smiles}"
    answer_format = "{Y1:y1, Y2:y2, Y3:y3, Y4:y4, Y5:y5}"
    question = "Answer 5 candidates of 'y1 to y5' in \
                {answer_format}"
    condition_1 = "Don't use '\n' in your ansewer" #これは微妙
    condition_2 = "Exclude unclosed ring structures from y1-y5 \
                    #(example:'C1=CC=C(CNC2=NC=C(C=C2)C(=O)N)C(F)(F)F')." #これはなし
    condition_3 = f"If {test_A_smiles} and {test_B_smiles} are in {training_dataset}, \
                    the corresponding 'Y' should be included in y1-y5." #これは効果あった?
    condition_4 = f"Count the number of atoms of the compounds A, B, and Y \
                    in the {training_dataset} set and learn the changes before and after the reaction \
                    example: \
                    A: c1ccc(N)cc1 (number of atoms is 14) \
                    B: c1ccc(Br)cc1 (number of atoms is 12) \
                    Y: c1ccc(NC2C=CC=CC=2)cc1 (number of atoms is 24)."
    condition_5 = f"Atoms and functional groups of Y that not involved \
                    in the reaction remain unchanged from {test_A_smiles} and {test_B_smiles}."

    response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-16k",
    messages=[
        {"role": "system", "content": "Synthesize compound Y from test data compounds A and B."},
        {"role": "assistant", "content": f"{training_dataset}+{test}"},
        {"role": "user", "content": f"{question}"},
        {"role": "assistant", "content": f"{condition_1}"},
        {"role": "assistant", "content": f"{condition_2}"},
        {"role": "assistant", "content": f"{condition_3}"},
        #{"role": "assistant", "content": f"{condition_4}"},
        #{"role": "assistant", "content": f"{condition_5}"}
        ],
        max_tokens=300,
        temperature=0.5,
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
    
    test_A_mol = Chem.MolFromSmiles(test_A_smiles)
    test_A_maccs_fps = AllChem.GetMACCSKeysFingerprint(test_A_mol)

    #テスト化合物Aを基準にタニモト係数を計算
    df_product_Y_candidates["Y_candidates_tnmt"] = \
        df_product_Y_candidates["Y_candidates_maccs_fps"].\
        apply(lambda maccs_fps: DataStructs.TanimotoSimilarity(test_A_maccs_fps, maccs_fps)\
        if maccs_fps is not None else None
        )
    
    #計算したタニモト係数を降順に並べる
    df_Y = df_product_Y_candidates.sort_values("Y_candidates_tnmt", ascending=False)

    #df_product_Y_candidates.iloc[:, 0].reset_index()

    return df_Y