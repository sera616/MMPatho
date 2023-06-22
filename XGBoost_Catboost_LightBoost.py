#!/usr/bin/env python
# coding: utf-8

def comEvaluation(y_true,y_pred):
    from sklearn.metrics import confusion_matrix
    #####------------  tn，fp, fn, tp  ------------###
    def tn(y_true, y_pred):
        return metrics.confusion_matrix(y_true, y_pred)[0, 0]
    def fp(y_true, y_pred):
        return metrics.confusion_matrix(y_true, y_pred)[0, 1]
    def fn(y_true, y_pred):
        return metrics.confusion_matrix(y_true, y_pred)[1, 0]
    def tp(y_true, y_pred):
        return metrics.confusion_matrix(y_true, y_pred)[1, 1]
    TN = tn(y_true, y_pred)
    FP = fp(y_true, y_pred)
    TP = tp(y_true, y_pred)
    FN = fn(y_true, y_pred)
    #sensitivity, recall, hit rate, true positive rate ：TPR = TP / (TP + FN)
    SN = TP*1.0/(TP + FN)*1.0 ## 
    #specificity, true negative rate:TNR = TN / (TN + FP)
    SP = TN / (TN + FP)*1.0  ## 
    #precision, prositive predictive value:PPV = TP / (TP + FP)
    precision = TP / (TP + FP)*1.0
    #negative predictive value:NPV = TN / (TN + FN)
    NPV = TN / (TN + FN)*1.0
    # F1 score is the harmonic mean of precision and sensitivity
    F1= 2*TP / (2*TP + FP+FN)*1.0
    # recall
    recall= TP/(TP+FN)*1.0
    return SN, SP,precision, NPV,F1,recall,TP, TN, FP, FN
print('define sucess!')


import catboost
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from catboost import *
import shap
from sklearn.metrics import accuracy_score
shap.initjs()

## load data
data = pd.read_excel('../data/3--VEP_output/4--Missense_Mutation_GRch37_for_VEP_info_output_mapped_features.xlsx', sheet_name='individual_outputs')
X =  data.loc[:,'SIFT_pred':'MutPred_Top5features_5']
y =  data['label']
print(X.shape, y.shape)
from sklearn.model_selection import train_test_split, GridSearchCV
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=33)

## --------  CatBoost --------##
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score
from catboost import CatBoostClassifier
param_grid = {
    'iterations': [100, 200, 300],
    'learning_rate': [0.1, 0.01, 0.001],
    'depth': [4, 6, 8],}
# Create an instance of CatBoostClassifier
model = CatBoostClassifier()
# Perform grid search with cross-validation
grid_search = GridSearchCV(model, param_grid, cv=5)
grid_search.fit(X_train, y_train)
# Get the best model from grid search
best_model = grid_search.best_estimator_
# Make predictions on the testing data
y_pred = best_model.predict(X_test)
# Calculate the accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
# Get the best parameters from grid search
best_params = grid_search.best_params_
print('Catboost的最优参数是：\t', best_params)
y_pred_proba_catboost = best_model.predict_proba(X_test)
print(y_pred_proba_catboost.shape, y_pred_proba_catboost[:5])

## 
from sklearn.model_selection import train_test_split
import xgboost
import shap
import numpy as np
import pandas as pd
from numpy import array
import matplotlib.pylab as pl
get_ipython().run_line_magic('matplotlib', 'inline')
shap.initjs()
d_train = xgboost.DMatrix(X_train, label=y_train)
d_test = xgboost.DMatrix(X_test, label=y_test)
print(X_train.shape, y_train.shape)
print(X_test.shape, y_test.shape)
params = {
    "eta": 0.01,
    "objective": "binary:logistic",
    "subsample": 0.5,
    "base_score": np.mean(y_train),
    "eval_metric": "logloss" }
model = xgboost.train(params, d_train, 5000, evals = [(d_test, "test")], verbose_eval=100, early_stopping_rounds=42)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
plt.rcdefaults()
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc('font',family='Arial')
plt.rcParams.update({"font.size":20})
fig = plt.figure(figsize=(50,20))
# plt.legend(loc='upper center',bbox_to_anchor=(0.40,1.01), ncol=2, fontsize=18)   
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25) 
# shap.summary_plot(shap_values, X_train, show = False, max_display = 7 , plot_type="violin")
xgboost.plot_tree(model, num_trees=0, rankdir='UT')
plt.rcParams['savefig.dpi'] = 7200 #图片像素
plt.rcParams['figure.dpi'] = 7200 #分辨率
plt.tight_layout()
plt.savefig('../data/4--figures/5--XGBoost_tree.pdf', dpi=7200) 
plt.show()
## ----------------------------- XGBoost functions ------------------------#
# Get all the functions in the Booster object
functions = [func for func in dir(model) if callable(getattr(model, func))]
# Print the list of functions
for func in functions:
    print(func)

## --------  LightGBM --------##
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score
from lightgbm import LGBMClassifier
# Define the parameter grid for grid search
param_grid = {
    'n_estimators': [100, 200, 300,400],
    'learning_rate': [0.1, 0.01, 0.001],
    'max_depth': [4, 6, 8], }

# Create an instance of LGBMClassifier
model = LGBMClassifier()
# Perform grid search with cross-validation
grid_search = GridSearchCV(model, param_grid, cv=5)
grid_search.fit(X_train, y_train)
# Get the best model from grid search
best_model = grid_search.best_estimator_
# Make predictions on the testing data
y_pred = best_model.predict(X_test)
# Calculate the accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)
# Get the best parameters from grid search
best_params = grid_search.best_params_
print('LightGBM optimal parameters：\t', best_params)
y_pred_proba_lightgbm = best_model.predict_proba(X_test)
print(y_pred_proba_lightgbm.shape, y_pred_proba_lightgbm[:5])
## save y_pred_proba
import pandas as pd
y_test_2 = pd.DataFrame(y_test)
y_pred_proba_catboost_2 = pd.DataFrame(y_pred_proba_catboost)
y_pred_proba_xgboost_2 = pd.DataFrame(y_pred_proba_xgboost)
y_pred_proba_lightgbm_2 = pd.DataFrame(y_pred_proba_lightgbm)
writer = pd.ExcelWriter('../data/3--VEP_output/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_prob.xlsx')
y_test_2.to_excel(writer,'y_test')
y_pred_proba_catboost_2.to_excel(writer,'y_pred_proba_catboost')
y_pred_proba_xgboost_2.to_excel(writer,'y_pred_proba_xgboost')
y_pred_proba_lightgbm_2.to_excel(writer,'y_pred_proba_lightgbm')
writer.close()
print('save finish！')

## 2.draw plot
import sklearn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc('font',family='Arial')
# plt.rc('font',family='Times New Roman')
plt.rcParams.update({"font.size":20})
fig = plt.figure(figsize=(28,12))
y_test2 = pd.read_excel('../data/3--VEP_output/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_prob.xlsx',
                        sheet_name='y_test')
y_pred_proba_catboost2 = pd.read_excel('../data/3--VEP_output/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_prob.xlsx',
                                       sheet_name='y_pred_proba_catboost')
y_pred_proba_xgboost2 = pd.read_excel('../data/3--VEP_output/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_prob.xlsx',
                                sheet_name='y_pred_proba_xgboost')
y_pred_proba_lightgbm2 = pd.read_excel('../data/3--VEP_output/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_prob.xlsx',
                               sheet_name='y_pred_proba_lightgbm')
y_test3 = y_test2.iloc[:,1]
y_test3 = y_test3.tolist()
y_pred_proba_catboost3 = y_pred_proba_catboost2.iloc[:,2]
y_pred_proba_catboost3 = y_pred_proba_catboost3.tolist()##
y_pred_proba_xgboost3 = y_pred_proba_xgboost2.iloc[:,1]
y_pred_proba_xgboost3 = y_pred_proba_xgboost3.tolist()##
y_pred_proba_lightgbm3 = y_pred_proba_lightgbm2.iloc[:,2]
y_pred_proba_lightgbm3 = y_pred_proba_lightgbm3.tolist()##
fpr_catboost, tpr_catboost, thresholds_catboost = sklearn.metrics.roc_curve(y_test3,y_pred_proba_catboost3)
roc_auc_catboost = auc(fpr_catboost, tpr_catboost)
fpr_xgboost, tpr_xgboost, thresholds_xgboost = sklearn.metrics.roc_curve(y_test3,y_pred_proba_xgboost3)
roc_auc_xgboost = auc(fpr_xgboost, tpr_xgboost)
fpr_lightgbm, tpr_lightgbm, thresholds_lightgbm = sklearn.metrics.roc_curve(y_test3,y_pred_proba_lightgbm3)
roc_auc_lightgbm = auc(fpr_lightgbm, tpr_lightgbm)
#### ------------------------------- ########################
ax = fig.add_subplot(121)
lw = 2
####################################
ax.plot( fpr_catboost, tpr_catboost,color="plum",lw=lw, label="CatBoost (AUROC = %0.4f)" % roc_auc_catboost)
ax.plot( fpr_xgboost, tpr_xgboost,color="peru",lw=lw, label="XGBoost (AUROC = %0.4f)" % roc_auc_xgboost)
ax.plot( fpr_lightgbm, tpr_lightgbm,color="darkred",lw=lw, label="LightGBM (AUROC = %0.4f)" % roc_auc_lightgbm)
#####################################
ax.plot([0, 1], [0, 1], color="navy", lw=lw/1.0, linestyle="--")
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.0])
ax.set_xlabel("False Positive Rate", fontsize=20)
ax.set_ylabel("True Positive Rate", fontsize=20)
# ax.set_title("On test data")
# ax.set_title("Receiver operating characteristic curve", fontsize=16)
ax.legend(loc="lower right", ncol=1, fontsize=18)
plt.tick_params(labelsize=20)

##############------------- precision-recall curve -------------##################
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import PrecisionRecallDisplay
from sklearn.metrics import average_precision_score
ax = fig.add_subplot(122)
lw = 2
#### ---- Del-Neu在三种策略下的预测结果 ---- ########################
prec_catboost, recall_catboost, _ = precision_recall_curve(y_test3,y_pred_proba_catboost3, pos_label=1)
AP_catboost = average_precision_score(y_test3,y_pred_proba_catboost3)
prec_xgboost, recall_xgboost, _ = precision_recall_curve(y_test3,y_pred_proba_xgboost3, pos_label=1)
AP_xgboost = average_precision_score(y_test3,y_pred_proba_xgboost3)
prec_lightgbm, recall_lightgbm, _ = precision_recall_curve(y_test3,y_pred_proba_lightgbm3, pos_label=1)
AP_lightgbm = average_precision_score(y_test3,y_pred_proba_lightgbm3)
#######################################################################
############################
ax.plot( prec_catboost, recall_catboost,color="plum",lw=lw, label="CatBoost (AUPR = %0.4f)" % AP_catboost)
ax.plot( prec_xgboost, recall_xgboost,color="peru",lw=lw, label="XGBoost (AUPR = %0.4f)" % AP_xgboost)
ax.plot( prec_lightgbm, recall_lightgbm,color="darkred",lw=lw, label="LightGBM (AUPR = %0.4f)" % AP_lightgbm)
############################
ax.plot([0, 1], [1, 0], color="navy", lw=lw/1.0, linestyle="--")
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.0])
ax.set_xlabel("Recall")
ax.set_ylabel("Precision")
# ax.set_title("On test data")
# ax.set_title("Precision-recall curve", fontsize=16)
ax.legend(loc="lower left", ncol=1, fontsize=18)
# ax.set_xlabel('(B) Precision-recall curve', fontsize=20)
plt.subplots_adjust(wspace =0.1, hspace =0)
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
plt.savefig('../data/4--figures/5--Missense_Mutation_100ind_catboost_xgboost_lightgbm_y_pred_ROC_PR.tif', dpi=300,bbox_inches ='tight') #指定分辨率保存
plt.show()

# encoding=utf-8
from matplotlib import pyplot
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
import numpy as np
#从pyplot导入MultipleLocator类，这个类用于设置刻度间隔
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc('font',family='Arial')
plt.rcParams.update({"font.size":14})

fig = plt.figure(figsize=(16, 4))
ax = plt.axes()
# 取消边框
for key, spine in ax.spines.items():
    # 'left', 'right', 'bottom', 'top'
    if key == 'right' or key == 'top' :
        spine.set_visible(False)
#-------------------- MM ----------------------##
x_data = [0.15, 1.15, 2.15, 3.15, 4.15, 5.15, 6.15]
##------------------- Catboost -------------------------##
Catboost = (0.9304, 0.9024, 0.9085, 0.9256, 0.9167, 0.8334, 0.9193)
##------------------- XGBoost -------------------------##
XGBoost = (0.9334, 0.9076, 0.9136, 0.9287, 0.9208, 0.8416, 0.9234)
##------------------- LightGBM -------------------------##
LightGBM = (0.9305, 0.9056, 0.9119, 0.9255, 0.9184, 0.8368, 0.9211)
bar_width = 0.3
plt.bar(np.arange(len(x_data))+0.05,Catboost, label='Catboost',color='#EEE8AA',width=bar_width,ec='#A05220',
       linewidth=1, ecolor='#800000', align='edge', bottom=0, alpha=0.8, linestyle='-',
       fill=True, joinstyle='bevel', hatch='/')
plt.bar(np.arange(len(x_data))+bar_width+0.05,XGBoost, label='XGBoost',color= '#B4D38B',width=bar_width,ec='#A05220',
       linewidth=1, ecolor='#800000', align='edge', bottom=0, alpha=0.8, linestyle='-',
       fill=True, joinstyle='bevel', hatch='\\\\')
plt.bar(np.arange(len(x_data))+2*bar_width+0.05,LightGBM, label='LightGBM',color='#CCCC99',width=bar_width,ec='#A05220',
       linewidth=1, ecolor='#800000', align='edge', bottom=0, alpha=0.8, linestyle='-',fill=True, joinstyle='bevel', 
        hatch='\\')
for x,y in enumerate(Catboost):
    plt.text(x+0.18, y-0.0001,'%0.4f'%y,ha='center',va='bottom',fontsize=11)    
for x,y in enumerate(XGBoost):
    plt.text(x+bar_width+0.18, y+0.001, '%0.4f'%y,ha='center',va='bottom',fontsize=11)    
for x,y in enumerate(LightGBM):
    plt.text(x+bar_width*2+0.18, y+0.0010, '%0.4f'%y,ha='center',va='bottom',fontsize=11)  

ax.set_ylabel('Values')
ax.set_ylim(0.8,0.95)
index = [0.48, 1.53, 2.53, 3.53, 4.53, 5.53, 6.53]
ax.set_xticks(index )
ax.set_xticklabels(('Pre', 'NPV','Recall', 'Spe','ACC', 'MCC',r'$F_1$'),fontstyle='italic')
# ax1.set_xlabel('(B) HumVar')
# fig.tight_layout()
# plt.subplots_adjust(wspace=0, hspace=0.4)
plt.legend(loc='upper center',bbox_to_anchor=(0.5,1.1), ncol =3, fontsize=13)
# plt.savefig('./8--Images/0--KVIST-FFMSRes-MutP--minor-revision.eps',dpi=300,bbox_inches ='tight')
plt.savefig('../data/4--figures/5--CatBooat_XGBoost_LightGBM_comaprisons_bar.tif',dpi=300,bbox_inches ='tight')
# plt.savefig('./8--Images/0--KVIST-FFMSRes-MutP--minor-revision.pdf',dpi=300,bbox_inches ='tight')
plt.rcParams['savefig.dpi'] = 300 #
plt.rcParams['figure.dpi'] = 300 #
plt.show()