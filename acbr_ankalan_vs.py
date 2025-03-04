########################################
# acbr_anakalan_vs.py

# This code is a part of Acbr_Ankalan webserver
# Github- https://github.com/mchopradu

# This code is owned by Prof. Madhu Chopra
# ACBR, University of Delhi, Delhi, India
# Date- 02/24/2025

#How to run this file 
# python acbr_ankalan_vs.py --target BCR-ABL --file input_smiles.csv --output screening_results.csv

#########################################

import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import argparse
import pandas as pd
import os
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Lipinski
from padelpy import padeldescriptor

# Delete output file if it already exists
if os.path.exists("screening_results.csv"):
    os.remove("screening_results.csv")

################
try:
    from dependencies import *  # Import all dependencies from dependencies.py
    print("✅ All libraries loaded successfully!")
except ImportError as e:
    print(f"❌ Error importing libraries: {e}")

############### Argument Parsing ###############
parser = argparse.ArgumentParser(description="Run virtual screening on a batch of SMILES from a CSV file.")
parser.add_argument("--target", required=True, choices=["BCR-ABL", "HDAC6", "PARP1", "TELOMERASE"],
                    help="Select the target from BCR-ABL, HDAC6, PARP1, TELOMERASE.")
parser.add_argument("--file", required=True, help="Path to input CSV file containing SMILES and ids.")
parser.add_argument("--output", default="screening_results.csv", help="Path to save the output CSV file.")

args = parser.parse_args()

target = args.target
input_file = args.file
output_file = args.output

# Delete the output file if it already exists
if os.path.exists(output_file):
    os.remove(output_file)

if not os.path.exists(input_file):
    print(f"❌ Error: Input file '{input_file}' not found!")
    exit()

print(f"✔️ Running virtual screening for target: {target}")
print(f"📂 Input file: {input_file}")

############## Loading the Model #########################
def load_model(target):
    """Load the trained ML model and preprocessing pipeline."""
    model_path = f'./models/vs/{target}/{target}.joblib'
    pipeline_path = f'./models/vs/{target}/{target}.pkl'
    
    if not os.path.exists(model_path) or not os.path.exists(pipeline_path):
        raise FileNotFoundError(f"❌ Model files for {target} not found!")
    
    model = joblib.load(model_path)
    pipeline = joblib.load(pipeline_path)
    return model, pipeline

######  Calculating Descriptors Using Padelpy  #########
def get_descriptors(smiles, target):
    """Compute molecular descriptors for a given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Skip invalid SMILES

    descriptors = {
        "MolWt": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "HeavyAtoms": Descriptors.HeavyAtomCount(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "AromaticRings": Lipinski.NumAromaticRings(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol)
    }

    with open("molecule.smi", "w") as f:
        f.write(f"{smiles}\tmol_001\n")

    padel_output = 'descriptors.csv'
    padeldescriptor(mol_dir='molecule.smi',
                    d_file=padel_output,
                    descriptortypes=f'./models/vs/{target}/{target}Fingerprinter.xml',
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True)

    if not os.path.exists(padel_output):
        return None, None  # If descriptor file not generated, skip molecule

    descriptors_df = pd.read_csv(padel_output)
    descriptors_df.drop(columns=['Name'], errors='ignore', inplace=True)

    return descriptors, descriptors_df

########## Prediction Function ##########
def predict_activity(smiles, target):
    """Predict the activity of a given molecule."""
    model, pipeline = load_model(target)
    descriptors, descriptors_df = get_descriptors(smiles, target)

    if descriptors_df is None:
        return None  # Skip invalid molecules

    if target == "BCR-ABL":
        X_test = pipeline.named_steps['variance_threshold'].transform(descriptors_df)
        prediction = model.predict_proba(X_test)[0, 1]

    elif target == "HDAC6":
        X_test_hdac_lgb = pipeline.named_steps['variance_threshold'].transform(descriptors_df)
        X_test_reduced_hdac_lgb = pipeline.named_steps['rfe'].transform(X_test_hdac_lgb)
        prediction = model.predict_proba(X_test_reduced_hdac_lgb)[0, 1]

    elif target == "PARP1":
        X_test_parp_xgb = pipeline.named_steps['variance_threshold'].transform(descriptors_df)
        X_test_reduced_parp_xgb = pipeline.named_steps['rfe'].transform(X_test_parp_xgb)
        prediction = model.predict_proba(X_test_reduced_parp_xgb)[0, 1]

    elif target == "TELOMERASE":
        X_test_telo_svc = pipeline.named_steps['variance_threshold'].transform(descriptors_df)
        X_test_reduced_telo_svc = pipeline.named_steps['rfe'].transform(X_test_telo_svc)
        prediction = model.predict_proba(X_test_reduced_telo_svc)[0, 1]

    return round(prediction, 2), descriptors  # Return both probability and descriptors

########## Running Virtual Screening ##########

try:
    df = pd.read_csv(input_file)
    if "id" not in df.columns or "SMILES" not in df.columns:
        raise ValueError("❌ Input CSV must have 'id' and 'SMILES' columns.")

    results = []
    for index, row in df.iterrows():
        smiles = row["SMILES"]
        molecule_id = row["id"]

        prob, descriptors = predict_activity(smiles, target)
        if prob is not None and descriptors is not None:
            prediction_status = "Active" if prob >= 0.5 else "Inactive"
            results.append({
                "id": molecule_id,
                "SMILES": smiles,
                "Prediction Probability": prob,
                "Prediction Status": prediction_status,  # New column for Active/Inactive
                "MolWt": descriptors["MolWt"],
                "LogP": descriptors["LogP"],
                "TPSA": descriptors["TPSA"],
                "HeavyAtoms": descriptors["HeavyAtoms"],
                "HBA": descriptors["HBA"],
                "AromaticRings": descriptors["AromaticRings"],
                "HBD": descriptors["HBD"],
                "RotatableBonds": descriptors["RotatableBonds"]
            })

    results_df = pd.DataFrame(results)

    if results_df.empty:
        print("❌ No valid predictions were made. Check input file!")
    else:
        with open(output_file, "w") as f:
            f.write(f"# Target: {target}\n")  # Add the target name as a header comment
        results_df.to_csv(output_file, index=False, mode="a")  # Append results below the header
        print(f"✔️ Virtual screening completed! Results saved to {output_file}")

except Exception as e:
    print(f"❌ Error: {e}")