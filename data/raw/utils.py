import os
import pandas as pd
from tqdm import tqdm

def show_stats():
    table = pd.read_csv(os.path.join(os.getcwd(), "activities_ic50.csv"), sep=";")
    print(table)
    print(table["Standard Relation"].value_counts())
    print(table["Standard Units"].value_counts())

def download_fasta():
    table = pd.read_csv(os.path.join(os.getcwd(), "targets.csv"), sep=";")
    for index, row in table.iterrows():
        id = row["UniProt Accessions"]
        id = id.split('|')[0]
        os.system(f'wget -P ./fasta/ https://rest.uniprot.org/uniprotkb/{id}.fasta')

def transform():
    table = pd.read_csv(os.path.join(os.getcwd(), "activities_ic50.csv"), sep=";")
    target = pd.read_csv(os.path.join(os.getcwd(), "targets.csv"), sep=";")
    result = open("seq2pot.csv", "w")
    result.write("Molecule ID,Target ID,IC50,Potency,Molecule,Target\n")
    for index, row in tqdm(table.iterrows()):
        if row["Standard Relation"] != "'='" or row["Standard Units"] != "nM":
            continue
        if pd.isnull(row["Standard Value"]):
            potency = 0
        elif row["Standard Value"] < 0.01:
            potency = 100
        else:
            potency = 1 / row["Standard Value"]
        uniprot = target.loc[target["ChEMBL ID"]==row["Target ChEMBL ID"], "UniProt Accessions"].iloc[0]
        uniprot = uniprot.split('|')[0]
        with open(os.path.join(os.getcwd(), "fasta", f"{uniprot}.fasta")) as f:
            lines = f.readlines()
            lines = lines[1:]
            seq = ''.join(lines)
            seq = seq.replace('\n', '')
        mol_id = row["Molecule ChEMBL ID"]
        tar_id = row["Target ChEMBL ID"]
        IC50 = row["Standard Value"]
        pot = potency
        mol = row["Smiles"]
        tar = seq
        result.write(f"{mol_id},{tar_id},{IC50},{pot},{mol},{tar}\n")
    result.close()



if __name__ == '__main__':
    #show_stats()
    #download_fasta()
    transform()