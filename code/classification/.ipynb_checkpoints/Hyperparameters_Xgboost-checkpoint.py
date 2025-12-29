import biom
import pandas as pd
import numpy as np

from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, roc_auc_score
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV

metadata = pd.read_csv("../../data/metadata.tsv", sep = "\t", encoding='gbk', index_col=0)
metadata = metadata.loc[metadata.study.values != "PRJEB2165"]
OTUs_table = biom.load_table("../../data/projects/table_6721_2378.biom")


base_model = XGBClassifier(
    random_state=42,
    eval_metric='logloss',
    use_label_encoder=False,
    n_jobs=-1  
)

param_dist = {
    'n_estimators': [100, 300, 500, 800],       
    'max_depth': [3, 4, 6, 8, 10],               
    'learning_rate': [0.01, 0.05, 0.1, 0.2],    
    'subsample': [0.7, 0.8, 0.9, 1.0],           
    'colsample_bytree': [0.7, 0.8, 0.9, 1.0],   
    'reg_alpha': [0.01, 0.1, 1.0],              
    'reg_lambda': [0.1, 1.0, 10.0]               
}

random_search = RandomizedSearchCV(
    estimator=base_model,
    param_distributions=param_dist,
    n_iter=50,  
    cv=5,      
    scoring='roc_auc',
    n_jobs=-1,
    random_state=42,
    verbose=1  
)

i = "PRJNA428736"
df = metadata.loc[metadata.study.values != i]
train_sid_1 = df.loc[df.group.values == 0].drop_duplicates(subset=['subject_id']).index.values
train_sid_2 = df.loc[df.group.values == 1].drop_duplicates(subset=['subject_id']).index.values
train_sid = list(train_sid_1) + list(train_sid_2)
test_sid = metadata.loc[metadata.study.values == i].index.values

X_train = OTUs_table.filter(train_sid, axis="sample", inplace=False)
X_train.remove_empty()
fid = X_train.ids(axis="observation")
X_train.norm(axis='sample')
X_train = X_train.to_dataframe().T

X_test = OTUs_table.filter(test_sid, axis="sample", inplace=False)
X_test = X_test.filter(fid, axis="observation", inplace=False)
X_test.norm(axis='sample')
X_test = X_test.to_dataframe().T

y_train = metadata.loc[train_sid].group.values
y_test = metadata.loc[test_sid].group.values


random_search.fit(X_train, y_train)

print("best Hyperparameters: ")
print(random_search.best_params_)

best_xgb_model = random_search.best_estimator_

y_pred = best_xgb_model.predict(X_test)
y_prob = best_xgb_model.predict_proba(X_test)

auc = roc_auc_score(y_test, y_prob[:, 1])
print(f"在测试集上的最终 AUC: {auc:.4f}")
