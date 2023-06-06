import openai
import os

openai.api_key_path = "./openai_api_key.txt"

def get_prodY_SMILES(A: str, B: str, context):

    context =f"{context}\
            test data:\
            {A} + {B} = '?'\
             "

    question = "Answer a candidate for '?'"

    response = openai.Completion.create(
        engine="text-davinci-003",
        prompt=f"Question answering:\nContext: {context}\nQuestion: {question}",
        max_tokens=50)
    
    product_Y = response.choices[0].text.strip()
    return product_Y
    
#sub_A = "CC1=NN(C=C1NC2=NC=C(C(=C2)I)C(F)(F)F)C"
#sub_B = "CONC(=O)C1=CC=CC=C1N"

#react_ai(sub_A, sub_B)