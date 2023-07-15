# ORD Chemical Reaction Predictor

[日本語版はこちら](https://github.com/KatsumiYamashita/React_ABY/blob/main/README_japanese.md)

You can check out the app on [hugging face](https://huggingface.co/spaces/kumasan681104/React_St) !

![image](https://github.com/KatsumiYamashita/React_ABY/blob/main/img/ord_finder_title_white.jpg?raw=true)

## About

This app can find organic reaction data containing your input compounds from the dataset of over 300,000 entries in [the Open Reaction Database](https://open-reaction-database.org/client/browse). 
Even when the compound is not present in the dataset, app attempts products prediction using GPT-3.5 trained on reaction SMILES.

## Usage

A simple usage of the app is as follows;

1. Let's access [app](https://huggingface.co/spaces/kumasan681104/React_St)!

2. Draw the structures of the two compounds whose reactions you want to search for! 
   Your compound information can be retrieved from PubChem database!

3. Press the "Find" button to explore reaction data from the Open Reaction Database!

## Datasets


## Future Features

For now, only two reactant reactions are shown.  
I will be adding more function.

- [ ] Achieve faster response.

- [ ] Find all reaction data from the ORD dataset.

- [ ] Enables searching by reaction.

- [ ] Improve product prediction performance.

# Packages

This app uses the following packages:

・[The Open Reaction Database](https://docs.open-reaction-database.org/en/latest/)

・[OpenAI](https://platform.openai.com/overview)

・[Streamlit](https://streamlit.io/)

・[streamlit-ketcher](https://github-com.translate.goog/mik-laj/streamlit-ketcher?ref=blog.streamlit.io&_x_tr_sl=en&_x_tr_tl=ja&_x_tr_hl=ja&_x_tr_pto=sc)

・[RDKit](https://www.rdkit.org/docs/index.html)

・[PubChem](https://pubchem.ncbi.nlm.nih.gov/)

・[Python](https://www.python.org/)

# Acknowledgment

[Suguru Tanaka](https://suguru-tanaka.com/)

I created my first St. App as part of his coaching program. 
I am truly grateful for his passionate guidance.

2023.07.17
Katsumi Yamashita

# License
