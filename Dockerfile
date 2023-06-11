# React-app/Dockerfile

FROM python:3.9

# set a directorys
WORKDIR /app

# install git and libXrender packages 
RUN apt-get update && apt-get install -y \
    curl \
    git \
    libxrender1

# clone code in a remote repo 
RUN git clone https://github.com/KatsumiYamashita/React_ABY.git

# set a directorys
WORKDIR /app/React_ABY

# install python dependencies from requirements.txt
RUN pip3 install -r requirements.txt

# run app on port ****
EXPOSE 8501


HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health


ENTRYPOINT ["streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]