from fastapi import FastAPI,HTTPException
import uvicorn
import pandas as pd
from modules import *
from pydantic import BaseModel
from fastapi.middleware.cors import CORSMiddleware
from Bio import SeqIO
from Bio.Seq import Seq
import random
import py3Dmol
import requests
import biotite.structure.io as bsio
import base64
from io import BytesIO
import re

zika_sequences = []
with open("./assets/zika_sequences.fasta", "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        zika_sequences.append(str(seq_record.seq))
chikungunya_sequences = []
with open("./assets/chikungunya_sequences.fasta", "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        chikungunya_sequences.append(str(seq_record.seq))
influenza_a_sequences = []
with open("./assets/influenza_a_sequences.fasta", "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        influenza_a_sequences.append(str(seq_record.seq))
influenza_b_sequences = []
with open("./assets/influenza_b_sequences.fasta", "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        influenza_b_sequences.append(str(seq_record.seq))
influenza_c_sequences = []
with open("./assets/influenza_c_sequences.fasta", "r") as fasta_file:
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        influenza_c_sequences.append(str(seq_record.seq))

zika_sequences = random.sample(zika_sequences,100)
chikungunya_sequences = random.sample(chikungunya_sequences,100)
influenza_a_sequences = random.sample(influenza_a_sequences,100)
influenza_b_sequences = random.sample(influenza_b_sequences,100)
influenza_c_sequences = random.sample(influenza_c_sequences,100)

viral_data = {
    "zika":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
    "chikungunya":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
    "influenza_a":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
    "influenza_b":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
    "influenza_c":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
    "unknown":{
        "gc_score":[],
        "at_score":[],
        "molecular_weight":[],
        "hydrophobic_score":[],
        "hydrophilic_score":[],
        "sequence_entropy":[],
    },
}

def create_data():
    for seq in zika_sequences:
        viral_data["zika"]["gc_score"].append(gc_content(seq))
        viral_data["zika"]["at_score"].append(at_content(seq))
        viral_data["zika"]["molecular_weight"].append(molecular_weight(seq))
        viral_data["zika"]["hydrophobic_score"].append(hydrophobicity(seq)[0])
        viral_data["zika"]["hydrophilic_score"].append(hydrophobicity(seq)[1])
        viral_data["zika"]["sequence_entropy"].append(sequence_entropy(seq))
    for seq in chikungunya_sequences:
        viral_data["chikungunya"]["gc_score"].append(gc_content(seq))
        viral_data["chikungunya"]["at_score"].append(at_content(seq))
        viral_data["chikungunya"]["molecular_weight"].append(molecular_weight(seq))
        viral_data["chikungunya"]["hydrophobic_score"].append(hydrophobicity(seq)[0])
        viral_data["chikungunya"]["hydrophilic_score"].append(hydrophobicity(seq)[1])
        viral_data["chikungunya"]["sequence_entropy"].append(sequence_entropy(seq))
    for seq in influenza_a_sequences:
        viral_data["influenza_a"]["gc_score"].append(gc_content(seq))
        viral_data["influenza_a"]["at_score"].append(at_content(seq))
        viral_data["influenza_a"]["molecular_weight"].append(molecular_weight(seq))
        viral_data["influenza_a"]["hydrophobic_score"].append(hydrophobicity(seq)[0])
        viral_data["influenza_a"]["hydrophilic_score"].append(hydrophobicity(seq)[1])
        viral_data["influenza_a"]["sequence_entropy"].append(sequence_entropy(seq))
    for seq in influenza_b_sequences:
        viral_data["influenza_b"]["gc_score"].append(gc_content(seq))
        viral_data["influenza_b"]["at_score"].append(at_content(seq))
        viral_data["influenza_b"]["molecular_weight"].append(molecular_weight(seq))
        viral_data["influenza_b"]["hydrophobic_score"].append(hydrophobicity(seq)[0])
        viral_data["influenza_b"]["hydrophilic_score"].append(hydrophobicity(seq)[1])
        viral_data["influenza_b"]["sequence_entropy"].append(sequence_entropy(seq))
    for seq in influenza_c_sequences:
        viral_data["influenza_c"]["gc_score"].append(gc_content(seq))
        viral_data["influenza_c"]["at_score"].append(at_content(seq))
        viral_data["influenza_c"]["molecular_weight"].append(molecular_weight(seq))
        viral_data["influenza_c"]["hydrophobic_score"].append(hydrophobicity(seq)[0])
        viral_data["influenza_c"]["hydrophilic_score"].append(hydrophobicity(seq)[1])
        viral_data["influenza_c"]["sequence_entropy"].append(sequence_entropy(seq))

def render_mol(pdb):
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.setBackgroundColor('white')
    view.zoomTo()
    view.spin(True)
    return view

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
class SequenceInput(BaseModel):
    sequence: str


@app.get("/")
async def root():
    return {"message": "Welcome to the Viral Gene Prediction API!"}

# Prediction endpoint
@app.post("/predict/")
async def predict(input_data: SequenceInput):
    sequence = input_data.sequence
    # print(sequence)
    df = feature_extraction_pipeline(sequence)
    df = scaler.transform(df)
    prediction = model.predict_proba(df)
    # print(prediction)
    return {
    "prediction": prediction.tolist()
    }

create_data()
@app.post("/getData")
async def getData(input_data: SequenceInput):
    if len(viral_data["unknown"]["gc_score"]) == 0:
        viral_data["unknown"]["gc_score"].append(gc_content(input_data.sequence))
        viral_data["unknown"]["at_score"].append(at_content(input_data.sequence))
        viral_data["unknown"]["molecular_weight"].append(molecular_weight(input_data.sequence))
        viral_data["unknown"]["hydrophobic_score"].append(hydrophobicity(input_data.sequence)[0])
        viral_data["unknown"]["hydrophilic_score"].append(hydrophobicity(input_data.sequence)[1])
        viral_data["unknown"]["sequence_entropy"].append(sequence_entropy(input_data.sequence))
    # print(viral_data)
    return viral_data

@app.post("/getStructure")
async def getStructure(input_data:SequenceInput):
    dna = re.sub(r'[^ATCG]', '', input_data.sequence)
    dna_seq = Seq(dna)
    rna = dna_seq.transcribe()
    rna_sequence = re.sub(r'[^AUCG]', '', str(rna))
    rna = Seq(rna_sequence)
    protein = rna.translate()
    def clean_sequence(seq, valid_chars):
        return ''.join([char for char in seq if char in valid_chars])
    valid_amino_acids = "ARNDCEQGHILKMFPSTWYV"
    protein = clean_sequence(protein, valid_amino_acids)[:400]
    print(protein)
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=protein)

    if response.status_code != 200:
        raise HTTPException(status_code=500, detail="Error fetching protein structure from the API.")

    pdb_string = response.content.decode('utf-8')

    # Return the structure as a response
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_string, 'pdb')
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    view.render()
    return {"pdb": pdb_string}

if __name__ == '__main__':
    uvicorn.run(app, host="127.0.0.1", port=8000)

