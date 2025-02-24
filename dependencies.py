# dependencies.py

""" This file contains centralized imports for `acbr_ankalan_cls.py` and `acbr_ankalan_vs.py`.

Author: Prof. Madhu Chopra 
Organization: ACBR, University of Delhi (DU), Delhi, India
Github: https://github.com/mchopradu
Date: 2/24/2025 
"""


import os, pathlib
import pandas as pd, numpy as np
from pathlib import Path
from shutil import which
from os.path import abspath, dirname, join
import subprocess, base64, time
from padelpy import padeldescriptor

from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
        roc_auc_score, roc_curve, RocCurveDisplay, accuracy_score, precision_score, recall_score,
        precision_recall_curve, average_precision_score, classification_report, confusion_matrix, auc,
        matthews_corrcoef
    )
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from imblearn.over_sampling import SMOTE
from xgboost import XGBClassifier as xgb
from sklearn.feature_selection import VarianceThreshold, RFE
from sklearn.decomposition import PCA
from lightgbm import LGBMClassifier as lgb

from PIL import Image
import joblib
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw
from rdkit.Chem import Descriptors, PandasTools, QED, Lipinski, MolFromSmiles
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.DataStructs import TanimotoSimilarity

    
