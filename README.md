![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)  


# ACBR_Ankalan

Title of the research article-

### ACBR_Ankalan Webserver: A novel Machine Learning-enabled platform for advancing Drug Discovery of BCR-ABL, HDAC6, PARP1, and Telomerase Inhibitors

This repository is a part of the abovementioned research article and contain details, codes and the associated ML models of the webserver deployed on the following link [ACBR Ankalan](https://bic.acbr.du.ac.in/). It is available for any one to use it freely for educational and research purposes.

### DOI of the research article- 
Coming soon

## Description 
ACBR_Aankalan is an open-source webserver designed to help with target-driven early-stage drug discovery. It uses proprietary Machine Learning classifiers to predict the bioactivity class (active or inactive) of the molecules against various cancer targets.  With its user-friendly interface , ACBR_Aankalan enables researchers to predict, conduct virtual screening, and use the tool into larger drug discovery pipelines.

ACBR_Aankalan has two types of modules for each target (BCR-ABL, HDAC6, PARP1 and Telomerase)-

### (1) Classification- For bioactivity class (active or inactive) prediction of a single compound against any of the selected target.
### (2) Virtual Screening- For bioactivity class (active or inactive) prediction of more than one compounds against any of the selected target.

For both modules, the prediction results are supplemented with key ADME properties, including Molecular Weight, Aromaticity, LogP, Heavy Atom Count, Hydrogen Bond Acceptors (HBA), Hydrogen Bond Donors (HBD), Total Polar Surface Area (TPSA), and Rotatable Bonds—all essential characteristics for assessing a compound's drug-likeness.

The classification module also includes a Tanimoto Similarity component, which compares the predicted active compound to known active compounds in the training dataset.  This aids in the identification of the most closely related active compound, revealing further information about its potential bioactivity.

The Acbr_Ankalan repository on github is organized into two primary modules:

Classification/ → Contains classification models and scripts.
Virtual_screening/ → Contains virtual screening models and scripts.
Each module has separate subfolders for different targets: BCR-ABL, HDAC6, PARP1, and TELOMERASE as shown below

```bash

Acbr_Ankalan/                # Main repository folder
│── classification/          # Classification module
│   ├── BCR-ABL/             # Models & data for BCR-ABL classification
│   ├── HDAC6/               # Models & data for HDAC6 classification
│   ├── PARP1/               # Models & data for PARP1 classification
│   ├── TELOMERASE/          # Models & data for Telomerase classification
│── virtual_screening/       # Virtual screening module
│   ├── BCR-ABL/             # Models & data for BCR-ABL virtual screening
│   ├── HDAC6/               # Models & data for HDAC6 virtual screening
│   ├── PARP1/               # Models & data for PARP1 virtual screening
│   ├── TELOMERASE/          # Models & data for Telomerase virtual screening
│── dependencies.py          # Centralized imports for all scripts
│── requirements.txt         # Required dependencies for running the scripts
│── acbr_ankalan_vs.py       # Virtual screening script
│── acbr_ankalan_cls.py      # Classification script
│── README.md                # Project documentation


```

### Dataset

The compound datasets for these targets were obtained from [ChEMBL v33](https://www.ebi.ac.uk/chembl/), resulting in high-quality, bioactivity-annotated data for model training and validation.


## Dataset used for model development

| Target      | Raw Dataset | Unique Dataset | Active (pIC50 > 6.5, IC50 ≤ 300nM) | Inactive (pIC50 < 6.3, IC50 ≥ 500nM) |
|------------|------------|---------------|----------------------------------|----------------------------------|
| BCR-ABL    | 2225       | 1561          | 886                              | 675                              |
| HDAC6      | 4212       | 3053          | 1781                             | 1134                             |
| PARP1      | 2426       | 2013          | 1338                             | 600                              |
| Telomerase | 388        | 281           | 117                              | 164                              |


## Best Model Details for Each Target

| Cancer Target | Fingerprint   | Algorithm | Pipeline Steps                     | No. of Features | 5-CV and Test Accuracy                        |
|--------------|--------------|-----------|----------------------------------|----------------|------------------------------------------|
| BCR-ABL      | CDK          | LightGBM  | VarianceThreshold, SMOTE         | 868            | 5-CV =  87.78% and Test Accuracy = 89.78%  |
| HDAC6        | CDK          | LightGBM  | VarianceThreshold, SMOTE, RFE    | 50             | 5-CV =  86.61% and Test Accuracy = 86.79%  |
| PARP1        | PubChem      | XGBoost   | VarianceThreshold, SMOTE, RFE    | 200            | 5-CV =  90.56% and Test Accuracy = 89.95%  |
| Telomerase   | Klekotaroth  | SVC       | VarianceThreshold, SMOTE, RFE    | 50             | 5-CV =  84.52% and Test Accuracy = 87.27%  |


## How to use it

Since it contains two modules—Classification and Virtual Screening; predictions for each target can be made independently in both the modules.

## 🔹 How to Use This Repository 

To run the classification and virtual screening models, you can either **clone the repository using Git** (recommended) or **manually download the files** and set up the folder structure.

### **✅ Option 1: Clone the Repository Using Git (Recommended)**
This method ensures you have all necessary files in the correct structure.

```bash
# Clone the repository
git clone https://github.com/mchopradu/Acbr_Ankalan.git
cd Acbr_Ankalan

```

### **✅ Option 2: Downloading the files manually**
If downloading the files manually, go to the repository and click "Download ZIP". Extract the ZIP file on your computer.


**To run the prediction module, make folders and sub-folders like this**-

```bash

Acbr_Ankalan/classification/
├── BCR-ABL/
├── HDAC6/
├── PARP1/
├── TELOMERASE/

Place model files (.joblib, .pkl etc ) inside their respective subfolders.

# Create and activate a virtual environment (Optional but recommended)
python -m venv venv
source env/bin/activate   # On Linux/macOS
venv\\Scripts\\activate    # On Windows

# Install required libraries
pip install -r requirements.txt

From the Acbr_Ankalan main folder, run the classification script as shown below:

python acbr_ankalan_cls.py --target BCR-ABL --smiles "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"

Replace BCR-ABL with HDAC6, PARP1, or TELOMERASE and smiles as needed.

```

**Similarly, to perform virtual screening, make folders and sub-folders as follows-**

```bash

Acbr_Ankalan/virtual_screening/
├── BCR-ABL/
├── HDAC6/
├── PARP1/
├── TELOMERASE/

Place model files (.joblib, .pkl etc ) inside their respective subfolders.

# Create and activate a virtual environment (Optional but recommended)
python -m venv env
source env/bin/activate   # On Linux/macOS
env\\Scripts\\activate    # On Windows

# Install dependencies
pip install -r requirements.txt

**For Virtual_Screening module, use a CSV file with following format--- SMILES (1st column) and id (2nd column)**

From the Acbr_Ankalan main folder, run the virtual_screening script as shown below:

**If the smiles file is in the same directory as the script:**
python acbr_ankalan_vs.py --target BCR-ABL --file input_smiles.csv --output screening_results.csv

**If the smiles file is in some other folder:**
python acbr_ankalan_vs.py --target BCR-ABL --file path/to/input_smiles.csv --output screening_results.csv

Replace BCR-ABL with HDAC6, PARP1, or TELOMERASE and smiles as needed.

```

## Citation

If you find this resource helpful for your study or research, please do cite it:

Citation coming soon

## Acknowledgement

We gratefully acknowledge the Bioinformatics Centre (BIC), National Network Project, Department of Biotechnology, Govt. of India, and the Institution of Eminence, University of Delhi for financial support. We also acknowledge Dr. B. R. Ambedkar Center for Biomedical Research, University of Delhi, for providing a BIC facility.


## License  

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0)**.  

You are free to **use, modify, and distribute** this software under the terms of the GPL-3.0 license. However, any modified or derivative work **must also be released under GPL-3.0**.  

### Copyright  
Copyright (C) 2025 Prof. Madhu Chopra  
ACBR, University of Delhi, Delhi, India  

See the [LICENSE](./LICENSE) file for more details.






This work was supported by the Department of Biotechnology (DBT), Govt of India.
