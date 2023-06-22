#!/usr/bin/env python
# coding: utf-8

def comEvaluation(y_true,y_pred):
    from sklearn.metrics import confusion_matrix
    #####------------   定义：tn，fp, fn, tp  ------------###
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
    SN = TP*1.0/(TP + FN)*1.0 ## 也就是：SN
    #specificity, true negative rate:TNR = TN / (TN + FP)
    SP = TN / (TN + FP)*1.0  ## 也就是：SP
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

## train_test_split
from sklearn.model_selection import train_test_split
import xgboost
import shap
import numpy as np
import pandas as pd
from numpy import array
import matplotlib.pylab as pl
get_ipython().run_line_magic('matplotlib', 'inline')
# print the JS visualization code to the notebook
shap.initjs()
## load data
data = pd.read_excel('../data/3--VEP_output/4--Missense_Mutation_GRch37_for_VEP_info_output_mapped_features.xlsx', sheet_name='individual_outputs')
X =  data.loc[:,'SIFT_pred':'MutPred_Top5features_5']
y =  data['label']
print(X.shape, y.shape)
# create a train/test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=33)
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
model = xgboost.train(params, d_train, 5000, evals = [(d_test, "test")], verbose_eval=100, early_stopping_rounds=30)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
plt.rcdefaults()
%matplotlib inline
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
plt.savefig('../data/4--figures/5--XGBoost_tree.pdf', dpi=7200) #指定分辨率保存
plt.show()
print('finish !')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
plt.rcdefaults()
%matplotlib inline
plt.rc('font',family='Arial')
plt.rcParams.update({"font.size":10})
fig = plt.figure(figsize=(200,10))
xgboost.plot_importance(model)
plt.title("xgboost.plot_importance(model)")
plt.show()
xgboost.plot_importance(model, height=0.2, ylim = (0,30), max_num_features = 30, grid =False )
plt.show()
xgboost.plot_importance(model, importance_type="cover", height=0.2, ylim = (0,30), max_num_features = 30, grid =False )
pl.title('xgboost.plot_importance(model, importance_type="cover")')
pl.show()
xgboost.plot_importance(model, importance_type="gain", height=0.2, ylim = (0,30), max_num_features = 30, grid =False )
pl.title('xgboost.plot_importance(model, importance_type="gain")')
pl.show()

# # Explain predictions
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_train)
# summary_plot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
plt.rcdefaults()
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc('font',family='Arial')
plt.rcParams.update({"font.size":20})
fig = plt.figure(figsize=(11,8.5))
# plt.legend(loc='upper center',bbox_to_anchor=(0.40,1.01), ncol=2, fontsize=18)   
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25) 
plt.rcParams['savefig.dpi'] = 300 
plt.rcParams['figure.dpi'] = 300 
plt.tight_layout()
shap.summary_plot(shap_values, X_train, show = False, max_display = 7 , plot_type="violin")
plt.savefig('../data/4--figures/1--XGBoost_SHAP_feature_importance_top7.tif', dpi=300) #指定分辨率保存
print('finish !')
shap.summary_plot(shap_values, plot_type="layered_violin", color='coolwarm')
shap.summary_plot(shap_values, plot_type="layered_violin")
shap.summary_plot(shap_values, X, plot_type="violin")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
plt.rcdefaults()
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rc('font',family='Arial')
plt.rcParams.update({"font.size":10})
fig = plt.figure(figsize=(11,10))
# plt.legend(loc='upper center',bbox_to_anchor=(0.40,1.01), ncol=2, fontsize=18)   
#plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.25) 
plt.rcParams['savefig.dpi'] = 300 #图片像素
plt.rcParams['figure.dpi'] = 300 #分辨率
plt.tight_layout()

shap.force_plot(explainer.expected_value, shap_values[0,:], X_test.iloc[0,:], show = False, matplotlib=True)
plt.savefig('../data/4--figures/2--XGBoost_SHAP_force_plot_0.tif', dpi=300) #指定分辨率保存

shap.force_plot(explainer.expected_value, shap_values[10,:], X_test.iloc[10,:], show = False, matplotlib=True)
plt.savefig('../data/4--figures/2--XGBoost_SHAP_force_plot_10.tif', dpi=300) #指定分辨率保存
print('finish !')
## 多个force_plot_many
xx = shap.force_plot(explainer.expected_value, shap_values[:1000,:], X_test.iloc[:1000,:])
shap.save_html('../data/4--figures/2--XGBoost_SHAP_force_plot_many.html',xx)

for name in X_train.columns:
    shap.dependence_plot(name, shap_values, X_train, display_features=X_train, show = False)
    plt.savefig('../data/4--figures/0--dependence_plot/'+ name + '_dependence_plot.tif', dpi=300) #指定分辨率保存
    