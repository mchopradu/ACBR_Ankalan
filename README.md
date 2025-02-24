# ACBR_Ankalan

## Advancing Drug Discovery with Machine Learning Models for BCR-ABL, HDAC6, PARP1 and Telomerase Inhibitors.

## Description 
ACBR_Aankalan is an open-source web server designed to help with target-driven early-stage drug discovery. It uses proprietary Machine Learning classifiers to predict the bioactivity class (active or inactive) of the molecules against various cancer targets.  With its user-friendly interface , ACBR_Aankalan enables researchers to predict, conduct virtual screening, and use the tool into larger drug discovery pipelines.

ACBR_Aankalan has two types of modules for each target (BCR-ABL, HDAC6, PARP1 and Telomerase)-

### (1) Classification- For bioactivity class (active or inactive) prediction of a single compound against any of the selected target.
### (2) Virtual Screening- For bioactivity class (active or inactive) prediction of more than one compounds against any of the selected target.

For both modules, the prediction results are supplemented with key ADME properties, including Molecular Weight, Aromaticity, LogP, Heavy Atom Count, Hydrogen Bond Acceptors (HBA), Hydrogen Bond Donors (HBD), Total Polar Surface Area (TPSA), and Rotatable Bonds—all essential characteristics for assessing a compound's drug-likeness.

The classification module also includes a Tanimoto Similarity component, which compares the predicted active compound to known active compounds in the training dataset.  This aids in the identification of the most closely related active compound, revealing further information about its potential bioactivity.

# Dataset

## Dataset Details for Model Training and Testing

| Target      | Raw Dataset | Unique Dataset | Active (pIC50 > 6.5, IC50 ≤ 300nM) | Inactive (pIC50 < 6.3, IC50 ≥ 500nM) |
|------------|------------|---------------|----------------------------------|----------------------------------|
| BCR-ABL    | 2225       | 1561          | 886                              | 675                              |
| HDAC6      | 4212       | 3053          | 1781                             | 1134                             |
| PARP1      | 2426       | 2013          | 1338                             | 600                              |
| Telomerase | 388        | 281           | 117                              | 164                              |


## Citation

If you find this resource helpful for your study or research, please do cite it:

Citation coming soon

## Acknowledgement
This work was supported by the Department of Biotechnology (DBT), Govt of India.
