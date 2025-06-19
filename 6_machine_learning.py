import os
import re
import pandas as pd
import numpy as np

data_dir = 'mechine_learning' # 10 files
file_1 = pd.read_csv('mechine_learning/fun_ROI-193_degree.csv')
file_1.columns = ['participant_id', 'func_193_degree']
file_2 = pd.read_csv('mechine_learning/fun_ROI-192_dist.csv')
file_2.columns = ['participant_id', 'func_192_dist']
file_3 = pd.read_csv('mechine_learning/fun_ROI-193_dist.csv')
file_3.columns = ['participant_id', 'func_193_dist']
file_4 = pd.read_csv('mechine_learning/fun_ROI-189_Enodal.csv')
file_4.columns = ['participant_id', 'func_189_Enodal']
file_5 = pd.read_csv('mechine_learning/fun_ROI-190_Enodal.csv')
file_5.columns = ['participant_id', 'func_190_Enodal']
file_6 = pd.read_csv('mechine_learning/fun_ROI-193_Enodal.csv')
file_6.columns = ['participant_id', 'func_193_Enodal']
file_7 = pd.read_csv('mechine_learning/fiber_ROI-198_degree.csv')
file_7.columns = ['participant_id', 'fiber_198_degree']
#file_8 = pd.read_csv('mechine_learning/fiber_ROI-115_dist.csv') # 1
#file_8.columns = ['participant_id', 'fiber_115_dist']
file_9 = pd.read_csv('mechine_learning/fiber_ROI-29_Enodalwt.csv')
file_9.columns = ['participant_id', 'fiber_29_Enodalwt']
file_10 = pd.read_csv('mechine_learning/thickness_exclude_Subject_019.csv')

clinic_info = pd.read_csv('all_subs_152_info.csv')

file_12 = pd.merge(file_1, file_2, on='participant_id')
file_123 = pd.merge(file_12, file_3, on='participant_id')
file_1234 = pd.merge(file_123, file_4, on='participant_id')
file_12345 = pd.merge(file_1234, file_5, on='participant_id')
file_123456 = pd.merge(file_12345, file_6, on='participant_id')
file_1234567 = pd.merge(file_123456, file_7, on='participant_id')
#file_12345678 = pd.merge(file_1234567, file_8, on='participant_id')
file_123456789 = pd.merge(file_1234567, file_9, on='participant_id')
file_12345678910 = pd.merge(file_123456789, file_10, on='participant_id')


ml_base = clinic_info.loc[:, ['participant_id', 'subject', 'treatment', 'group']]
ml_base = ml_base[ml_base['subject']!='Subject_019']
#ml_base = ml_base[ml_base['treatment']=='before']
ml_base = ml_base.loc[:, ['participant_id', 'treatment', 'group']]

########### machine learning
ml_org = pd.merge(ml_base, file_12345678910, on='participant_id', how='inner')
ml_org.to_csv('mechine_learning/total_before_after.csv', index=None)
X_org = ml_org.iloc[:, 2:]
y_org = ml_org.loc[:, 'group']

from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_org)
from sklearn.preprocessing import LabelEncoder
encoder = LabelEncoder()
y_encoded = encoder.fit_transform(y_org)

from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y_encoded, test_size=0.5, random_state=42)

model = LogisticRegression()
model.fit(X_train, y_train)

# Step 5: Predict and evaluate the model
y_pred = model.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)

print("Accuracy:", accuracy)
