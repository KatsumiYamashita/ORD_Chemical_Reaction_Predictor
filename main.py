import openai
import re

def get_prodY_SMILES(A: str, B: str, training_dataset):

    context = f"A: {A}\
                B: {B}\
                Y1:\
                Y2:\
                Y3:\
                Y4:\
                Y5:\
                "

    response = openai.ChatCompletion.create(
    model="gpt-3.5-turbo-16k",
    messages=[
        #{"role": "system", "content": "Synthesize compound Y from test data compounds A and B."},
        {"role": "assistant", "content": f"{training_dataset}+\n\+{context}"},
        {"role": "user", "content": f"Answer at least five candidates for 'y1 to y5'."}],
        max_tokens=300,
        temperature=0,
        )
    
    response_text = response["choices"][0]["message"]["content"]
    pattern = r":\s(.+)" #テンプレートの作成
    product_Y_candidates = re.findall(pattern, response_text) #テンプレートをもとにSMILESを抜き出す
    return product_Y_candidates