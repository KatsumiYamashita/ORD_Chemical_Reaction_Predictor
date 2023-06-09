import openai

def get_prodY_SMILES(A: str, B: str, training_dataset):

    context =f"{training_dataset}\
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