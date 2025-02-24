# ACBR_Ankalan

## Advancing Drug Discovery with Machine Learning Models for BCR-ABL, HDAC6, PARP1 and Telomerase Inhibitors.

This repository contains details, codes and the associated ML models of the ACBR Ankalan web-server deployed on the following link [ACBR Ankalan](https://bic.acbr.du.ac.in/ankalan) for any one to use it freely for educational andd research purpose.

## Description 
ACBR_Aankalan is an open-source web server designed to help with target-driven early-stage drug discovery. It uses proprietary Machine Learning classifiers to predict the bioactivity class (active or inactive) of the molecules against various cancer targets.  With its user-friendly interface , ACBR_Aankalan enables researchers to predict, conduct virtual screening, and use the tool into larger drug discovery pipelines.

ACBR_Aankalan has two types of modules for each target (BCR-ABL, HDAC6, PARP1 and Telomerase)-

### (1) Classification- For bioactivity class (active or inactive) prediction of a single compound against any of the selected target.
### (2) Virtual Screening- For bioactivity class (active or inactive) prediction of more than one compounds against any of the selected target.

For both modules, the prediction results are supplemented with key ADME properties, including Molecular Weight, Aromaticity, LogP, Heavy Atom Count, Hydrogen Bond Acceptors (HBA), Hydrogen Bond Donors (HBD), Total Polar Surface Area (TPSA), and Rotatable Bonds—all essential characteristics for assessing a compound's drug-likeness.

The classification module also includes a Tanimoto Similarity component, which compares the predicted active compound to known active compounds in the training dataset.  This aids in the identification of the most closely related active compound, revealing further information about its potential bioactivity.

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


## Citation

If you find this resource helpful for your study or research, please do cite it:

Citation coming soon

## Acknowledgement
This work was supported by the Department of Biotechnology (DBT), Govt of India.
