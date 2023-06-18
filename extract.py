"""Extraction of Training Data from Open Reaction Datasets"""

import pandas as pd
import numpy as np
 
from ord_schema import message_helpers, validations
from ord_schema.proto import dataset_pb2
from rdkit import rdBase, Chem, DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem, Draw, rdMHFPFingerprint
from rdkit.Chem.Fingerprints import FingerprintMols

def import_dataset(ord_dataset_id):
    
    # データセットのパス
    dir_dataset = f"./ord_datasets/{ord_dataset_id}.pb.gz" 
    # データセットの読み込み
    dataset_raw = message_helpers.load_message\
                                    (dir_dataset, dataset_pb2.Dataset)
    # pandasデータフレームへの変換 
    df_dataset_raw = message_helpers.messages_to_dataframe\
                    (dataset_raw.reactions,drop_constant_columns=False)

    # dfからA+B→Yとなる部分だけ抜き出してDataFrameにする
    # 最大の問題点:データセットによって以下の取り出しキーワードが異なる
    # 最初にすべてのデータセットを読み込んでおいてごちゃ混ぜしてしおけばok?
    # とりあえずBuchwald-Hartwig Reaction dataset
    df_raw_ABY = df_dataset_raw[[ \
        'inputs["aryl halide"].components[0].identifiers[0].value', \
        'inputs["amine"].components[0].identifiers[0].value', \
        'outcomes[0].products[0].identifiers[0].value' \
        ]]
    
    df_raw_ABY.columns = list('ABY') #列ラベルをわかりやすく A,B,Yに変換
    
    # 重複データの消去
    # インデックスのリセット
    df_smiles_ABY = df_raw_ABY.drop_duplicates()\
                                .reset_index()\
                                .drop(columns=["index"])
    
    # rdkit Molオブジェクトを生成する
    # applymapからnumpy計算へ変更予定
    df_mol_ABY = df_smiles_ABY.applymap\
                (lambda smiles: Chem.MolFromSmiles(smiles))\
                .rename(columns={"A": "mol_A", "B": "mol_B", "Y": "mol_Y"})
    
    df_smiles_mol = pd.concat([df_smiles_ABY, df_mol_ABY], axis=1)

    return df_smiles_ABY, df_mol_ABY, df_smiles_mol

def generate_maccs_fps_df(df_mol_ABY, df_smiles_ABY):

    df_maccs_fps_ABY = \
    df_mol_ABY.applymap(lambda mol: AllChem.GetMACCSKeysFingerprint(mol))\
                .rename(columns={"mol_A": "maccs_fps_A",\
                                "mol_B": "maccs_fps_B",\
                                "mol_Y": "maccs_fps_Y"}
                        )

    df_smiles_maccs_fps = pd.concat([df_smiles_ABY, df_maccs_fps_ABY], axis=1)

    return df_maccs_fps_ABY, df_smiles_maccs_fps

def extract_training_data(df_smiles_maccs_fps,\
                            df_maccs_fps_ABY,\
                            reactant_A_maccs_fps,\
                            reactant_B_maccs_fps,\
                            number):
    
    # df_maccs_fps の各要素に対してタニモト係数を生成
    df_smiles_maccs_fps["tnmt_A"] = \
    df_maccs_fps_ABY["maccs_fps_A"].apply(lambda maccs_fps: \
    DataStructs.TanimotoSimilarity(reactant_A_maccs_fps, maccs_fps))

    df_smiles_maccs_fps["tnmt_B"] = \
    df_maccs_fps_ABY["maccs_fps_B"].apply(lambda maccs_fps: \
    DataStructs.TanimotoSimilarity(reactant_B_maccs_fps, maccs_fps))

    # 化合物Aのタニモト係数でソートする
    df_smiles_maccs_fps_tnmt = df_smiles_maccs_fps.sort_values("tnmt_A", ascending=False)
    # ソートしたdfから化合物Aに対するタニモト係数上位35個を抜き取る
    df_rank_tnmt_A_35 = df_smiles_maccs_fps_tnmt.iloc[:number, 0:3]
    # 抜き取ったdfからstr型のトレーニングデータを作る
    training_dataset_rank_tnmt_A = \
    "This is training dataset (A + B → Y):\n\
    "
    for _, row in df_rank_tnmt_A_35.iterrows():
        template = "A: " + row['A'] + "\\" + "\n" + "B: " + row['B'] + "\\" + "\n" + "Y: " + row['Y'] + "\\" + "\n" + "\\" + "\n"
        training_dataset_rank_tnmt_A += template
    
    return training_dataset_rank_tnmt_A



