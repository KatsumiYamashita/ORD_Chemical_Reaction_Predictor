# ORD Chemical Reaction Finder: 

アプリは [hugging face](https://huggingface.co/spaces/kumasan681104/React_St) からご覧ください!

You can check out the app on [hugging face](https://huggingface.co/spaces/kumasan681104/React_St) !

![image](https://github.com/KatsumiYamashita/React_ABY/blob/main/img/ord_finder_title_white.jpg?raw=true)

# このアプリについて (About)

このアプリでは, [the Open Reaction Database](https://open-reaction-database.org/client/browse)の300,000を超える有機合成反応データから,
あなたが入力した化合物を含む反応を検索できます.
化合物がデータセットに含まれなくても, 類似化合物の反応SMILESをトレーニングしたGPT-3.5を使って生成物の予測を試みました.

This app can find organic reaction data containing your input compounds from the dataset of over 300,000 entries in [the Open Reaction Database](https://open-reaction-database.org/client/browse). 
Even when the compound is not present in the dataset, app attempts products prediction using GPT-3.5 trained on reaction SMILES.

# 使用方法 (Usage)

アプリの簡単な使用方法は以下の通りです;
1. まずは, [アプリ](https://huggingface.co/spaces/kumasan681104/React_St) にアクセスしましょう!
2. 反応を検索したい二つの化合物の構造を描きましょう! 化合物の情報はPubChemのデータベースから取得できます.
3. "Find"ボタンを押してthe Open Reaction Databaseから反応データを探索しましょう!

A simple usage of the app is as follows;
1. Let's access [app](https://huggingface.co/spaces/kumasan681104/React_St)!
2. Draw the structures of the two compounds whose reactions you want to search for! 
   Your compound information can be retrieved from PubChem database!
3. Press the "Find" button to explore reaction data from the Open Reaction Database!

# 将来的機能　(Future features)

現在は2つの反応物間の反応データのみを表示しています.
さらに機能を追加していきます。

- [ ] より高速な応答を実現します。
- [ ] ORD データセットからすべての反応データを検索します。
- [ ] 反応による検索が可能になります。
- [ ] 製品予測パフォーマンスを向上させます。

For now, only two reactant reactions are shown.  
I will be adding more function.

- [ ] Achieve faster response.
- [ ] Find all reaction data from the ORD dataset.
- [ ] Enables searching by reaction.
- [ ] Improve product prediction performance.

# Packages

This app uses the following packages:

・[The Open Reaction Database](https://docs.open-reaction-database.org/en/latest/)

・[streamlit-ketcher](https://github-com.translate.goog/mik-laj/streamlit-ketcher?ref=blog.streamlit.io&_x_tr_sl=en&_x_tr_tl=ja&_x_tr_hl=ja&_x_tr_pto=sc)

・[PubChem](https://pubchem.ncbi.nlm.nih.gov/)

・[OpenAI](https://platform.openai.com/overview)

# Acknowledgment

[田中　統](https://suguru-tanaka.com/)

このアプリケーションは田中氏のコーチングプログラムに参加して制作しました. 
素人の私に対しても

[Suguru Tanaka](https://suguru-tanaka.com/)


# License
