# ORD Chemical Reaction Predictor: 

\[[English ver.](https://github.com/KatsumiYamashita/React_ABY/blob/main/README.md)\

アプリは [hugging face](https://huggingface.co/spaces/kumasan681104/React_St) からご覧ください!

![image](https://github.com/KatsumiYamashita/React_ABY/blob/main/img/ord_finder_title_white.jpg?raw=true)

### このアプリについて
　
このアプリでは, [the Open Reaction Database](https://open-reaction-database.org/client/browse)の300,000を超える有機合成反応データの中　から,あなたが入力した化合物を含む反応を検索できます.
化合物がデータセットに含まれなくても, 類似化合物の反応SMILESをトレーニングしたGPT-3.5を使って生成物の予測を試みました.

### 使用方法

アプリの簡単な使用方法は以下の通りです;

1. まずは, [アプリ](https://huggingface.co/spaces/kumasan681104/React_St) にアクセスしましょう!

2. 反応を検索したい二つの化合物の構造を描きましょう! 化合物の情報はPubChemのデータベースから取得できます.

3. "Find"ボタンを押してthe Open Reaction Databaseから反応データを探索しましょう!

### データセット



### 将来的機能

現在は2つの反応物間の反応データのみを表示しています.
さらに機能を追加していきます。

- [ ] より高速な応答を実現します。

- [ ] ORD データセットからすべての反応データを検索します。

- [ ] 反応による検索が可能になります。

- [ ] 製品予測パフォーマンスを向上させます。

### パッケージ

本アプリでは以下の機能・技術・データベースを使用しました:

・[The Open Reaction Database](https://docs.open-reaction-database.org/en/latest/)

・[OpenAI](https://platform.openai.com/overview)

・[Streamlit](https://streamlit.io/)

・[streamlit-ketcher](https://github-com.translate.goog/mik-laj/streamlit-ketcher?ref=blog.streamlit.io&_x_tr_sl=en&_x_tr_tl=ja&_x_tr_hl=ja&_x_tr_pto=sc)

・[RDKit](https://www.rdkit.org/docs/index.html)

・[PubChem](https://pubchem.ncbi.nlm.nih.gov/)

・[Python](https://www.python.org/)

# 謝辞

[田中　統](https://suguru-tanaka.com/)

このアプリケーションは田中氏のコーチングプログラムに参加して制作しました. 
3ヶ月前まで素人だった私がアプリまで制作できたのは田中さんの親切で解りやすい指導のおかげです.
この場を借りて感謝申し上げます.

2023.07.17
山下　黄

# License
