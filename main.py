import openai
import re

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
    
    response_text = response["choices"][0]["message"]["content"]
    pattern = r":\s(.+)" #テンプレートの作成
    product_Y_candidates = re.findall(pattern, response_text) #テンプレートをもとにSMILESを抜き出す
    return product_Y_candidates