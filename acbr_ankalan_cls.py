########################################
#acbr_anakalan_cls.py

# This code is a part of Acbr_Ankalan webserver
# Github- https://github.com/mchopradu

# This code is owned by Prof. Madhu Chopra
# ACBR, University of Delhi, Delhi, India
# Date- 02/24/2025

#########################################


import warnings
warnings.filterwarnings("ignore", category=UserWarning)



################
try:
    from dependencies import *  # Import all dependencies from dependencies.py
    print("✅ All libraries loaded successfully!")
except ImportError as e:
    print(f"❌ Error importing libraries: {e}")

###############
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Run classification on a given SMILES input.")
parser.add_argument("--target", required=True, choices=["BCR-ABL", "HDAC6", "PARP1", "TELOMERASE"],
                    help="Select the target from BCR-ABL, HDAC6, PARP1, TELOMERASE.")
parser.add_argument("--smiles", required=True, help="Enter a valid SMILES string.")

args = parser.parse_args()

# Extract arguments
target = args.target
smiles = args.smiles

# Validate input
if not smiles:
    print("❌ Warning: No SMILES input provided! Please enter a valid SMILES string.")
    exit()

if not target:
    print("❌ Warning: No target selected! Please select one from BCR-ABL, HDAC6, PARP1, TELOMERASE.")
    exit()

print(f"✔️ Running classification for target: {target} with SMILES: {smiles}")

#####################################



############## Loading the target model #########################

def load_model(target):
    """Load the trained ML model and preprocessing pipeline."""
    model_path = f'./models/classification/{target}/{target}.joblib'
    pipeline_path = f'./models/classification/{target}/{target}.pkl'
    
    if not os.path.exists(model_path) or not os.path.exists(pipeline_path):
        raise FileNotFoundError(f"Model files for {target} not found!")
    
    model = joblib.load(model_path)
    pipeline = joblib.load(pipeline_path)
    return model, pipeline


######  Calculating the descriptors using Padelpy #####################

def get_descriptors(smiles, target):
    """Compute molecular descriptors for a given SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string!")

    descriptors = {
        "MolWt": round(Descriptors.MolWt(mol), 2),
        "LogP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        "HeavyAtoms": Descriptors.HeavyAtomCount(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        #"QED": round(QED.qed(mol), 2),  
        "AromaticRings": Lipinski.NumAromaticRings(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "RotatableBonds": Lipinski.NumRotatableBonds(mol)
    }

    
    Draw.MolToFile(mol, 'molecule.png', width=900)

    
    with open("molecule.smi", "w") as f:
        f.write(f"{smiles}\tmol_001\n")

    # Compute PADEL descriptors
    padel_output = 'descriptors.csv'
    padeldescriptor(mol_dir='molecule.smi',
                    d_file=padel_output,
                    descriptortypes=f'./models/classification/{target}/{target}Fingerprinter.xml',
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True)

   
    if not os.path.exists(padel_output):
        raise FileNotFoundError(f"PADEL descriptor file {padel_output} was not generated!")

    descriptors_df = pd.read_csv(padel_output)
    descriptors_df.drop(columns=['Name'], errors='ignore', inplace=True)

    return descriptors, descriptors_df  


####################################### Performing classification prediction based on selected target ##########################

def predict_activity(smiles, target):
    """Predict the activity of a given molecule."""
    model, pipeline = load_model(target)
    descriptors, descriptors_df = get_descriptors(smiles, target)

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
        
    
    elif target=="TELOMERASE":
        X_test_telo_svc = pipeline.named_steps['variance_threshold'].transform(descriptors_df)
        X_test_reduced_telo_svc = pipeline.named_steps['rfe'].transform(X_test_telo_svc)
        prediction = model.predict_proba(X_test_reduced_telo_svc)[0, 1]

    return round(prediction, 2)

############################## Calculating the probability based on selected target ###################################
try:
   
    descriptors, descriptors_df = get_descriptors(smiles, target)  

    print("\n Molecular Descriptors:")
    for key, value in descriptors.items():
        print(f"{key}: {value}")

    prob = predict_activity(smiles, target)

    print("\n🚀 Prediction Result:")
    print(f"✔️ Active (Probability: {prob:.2f})" if prob >= 0.5 else f"❌ Inactive (Probability: {prob:.2f})")

    print("The structure of the compound is:")
    img = Image.open('molecule.png')
    display(img)

except Exception as e:
    print(f"Error: {e}")

############################# Calculating Morgan Fingerprint for Tanimoto similarity calculation

def get_morgan_fingerprint(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    return None


############################ Computing Tanimoto similarity #################################

def compute_tanimoto_similarity(query_fp, reference_smiles_list):
    similarities = []
    for index, row in reference_smiles_list.iterrows():
        smiles = row['Smiles']
        chembl_id = row['ChEMBL_ID_Molecule_Name']
        target_fp = get_morgan_fingerprint(smiles)
        
        if target_fp:
            similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
            similarities.append({
                'ChEMBL_ID_Molecule_Name': chembl_id,
                'Similarity': round(similarity, 2)
            })
    
    similarity_df = pd.DataFrame(similarities)
    similarity_df = similarity_df.sort_values(by='Similarity', ascending=False)
    return similarity_df

###################################

def process_tanimoto_similarity(target, prob):
    """Compute Tanimoto similarity for the selected target only if probability >= 0.5."""
    if prob < 0.5:
        print("\n❌ Probability is below threshold. Skipping Tanimoto similarity computation.")
        return
    
    print("\n✔️ Probability meets threshold. Computing Tanimoto similarity...")
    tanimoto_file_path = f'./models/classification/{target}/{target}_tanimoto_list.csv'
    
    if os.path.exists(tanimoto_file_path):
        tanimoto_df = pd.read_csv(tanimoto_file_path)
        query_file_path = 'molecule.smi'
        
        with open(query_file_path, 'r') as f:
            query_smiles, query_name = f.read().strip().split()
        
        query_fp = get_morgan_fingerprint(query_smiles)
        if query_fp:
            similarity_df = compute_tanimoto_similarity(query_fp, tanimoto_df)
            
            if not similarity_df.empty:
                top_row = similarity_df.iloc[0]
                print(f"\nTop {target} drug by Tanimoto similarity to the predicted active compound:")
                print(top_row)
            else:
                print("No valid matches found.")
        else:
            print("Invalid query SMILES. Please enter a valid SMILES string.")
    else:
        print(f"Tanimoto similarity file for {target} not found!")

##########################

prob = predict_activity(smiles, target)
process_tanimoto_similarity(target, prob)
