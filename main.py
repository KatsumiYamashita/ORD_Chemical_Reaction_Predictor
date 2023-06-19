import openai
import pandas as pd
import re

from rdkit import rdBase, Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

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
                    (example:'C1=CC=C(CNC2=NC=C(C=C2)C(=O)N)C(F)(F)F')." #これはなし
    condition_3 = f"If {test_A_smiles} and {test_B_smiles} are in {training_dataset}, \
                    the corresponding 'Y' should be included in y1-y5." #これは効果あった

    # エラーになった例を貯めていって学習させるっていうのはあり?

    response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-16k",
    messages=[
        {"role": "system", "content": "Synthesize compound Y from test data compounds A and B."},
        {"role": "assistant", "content": f"{training_dataset}+{test}"},
        {"role": "user", "content": f"{question}"},
        {"role": "assistant", "content": f"{condition_1}"},
        {"role": "assistant", "content": f"{condition_2}"},
        {"role": "assistant", "content": f"{condition_3}"}
        ],
        max_tokens=300,
        temperature=0,
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
    df_product_Y_candidates.sort_values("Y_candidates_tnmt", ascending=False)


    df_Y = df_product_Y_candidates.iloc[:, 0].reset_index()

    return df_Y