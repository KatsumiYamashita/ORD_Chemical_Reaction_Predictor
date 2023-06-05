# モジュールのインポート
import streamlit as st
import main

from rdkit import rdBase, Chem
from rdkit.Chem import Draw, rdChemReactions as Reactions
from rdkit.Chem.Draw import rdMolDraw2D, IPythonConsole
from IPython.display import SVG
from matplotlib.colors import ColorConverter

# アプリケーションタイトル
st.title("React: A + B → Y")

# とりあえずAとBが入力されたという設定で
reactant_A = "CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
reactant_B = "CONC(=O)C1=CC=CC=C1N"
product_Y = "CC1=NN(C=C1NC2=NC=C(C(=C2)NC3=CC=CC=C3C(=O)NC)C(F)(F)F)C"
mol_A = Chem.MolFromSmiles(reactant_A)
mol_B = Chem.MolFromSmiles(reactant_B)
mol_Y = Chem.MolFromSmiles(product_Y)

# 複数構造の描画はMolsToGridImage
#A_and_B2 = Draw.MolsToGridImage([mol_A, mol_B], molsPerRow=2, subImgSize=(300,300))
#st.write("複数構造の描画はMolsToGridImage")
#st.image(A_and_B2, use_column_width=True)

# より細かい設定をするならrdMolDraw2D
# コンテナと分子の準備
view = rdMolDraw2D.MolDraw2DSVG(440,350,220,350) #コンテナの作成 (440,350)の大きな箱の中に(220,350)の小さな箱
t_mol_A = rdMolDraw2D.PrepareMolForDrawing(mol_A)
t_mol_B = rdMolDraw2D.PrepareMolForDrawing(mol_B)
# オプションの設定
white = ColorConverter().to_rgb('white')
black = ColorConverter().to_rgb('black')

option = view.drawOptions()
option.padding=0.13
option.legendFontSize=18
#option.useBWAtomPalette() # 白黒構造で描画
option.setBackgroundColour(white) # 背景の色を変更

# コンテナに分子を3つ登録
view.DrawMolecules([t_mol_A,t_mol_B])   # コンテナに分子を格納 
view.FinishDrawing()                    # コンテナのファイナライズ
svg = view.GetDrawingText()             # 書き込んだデータをstr形式で取り出し

st.write("より細かい設定をするならrdMolDraw2D")
st.image(svg, use_column_width=True)

# rdMolDraw2D.MolDraw2DCairo
# MolDraw2DCairoオブジェクトを作成
#drawer = rdMolDraw2D.MolDraw2DCairo(220,350)  # 描画サイズを設定

# 分子を描画
#drawer.DrawMolecule(mol_A)

# 描画を保存
#cairo = drawer.GetDrawingText()

#st.write("rdMolDraw2D.MolDraw2DCairo")
#st.image(cairo, use_column_width=True)

# Reaction
rxn = Reactions.ReactionFromSmarts(f'{reactant_A}.{reactant_B}>>{product_Y}', useSmiles=True)
rxn_image = Draw.ReactionToImage(rxn)

st.image(rxn_image, use_column_width=True)

# Draw.Reaction
drawer = rdMolDraw2D.MolDraw2DSVG(660,350)
rxn = Reactions.ReactionFromSmarts(f'{reactant_A}.{reactant_B}>>{product_Y}', useSmiles=True)
drawer.DrawReaction(rxn)
drawer.FinishDrawing()
svg_rxn = drawer.GetDrawingText()   

st.write("rdMolDraw2D + svg + DrawReaction これが最強かも")
st.image(svg_rxn, use_column_width=True)