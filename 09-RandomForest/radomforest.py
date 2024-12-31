# -*- coding: UTF-8 -*-
# from sklearn.datasets import load_boston
# from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score  # 交叉检验
from sklearn.model_selection import GridSearchCV #网络搜索和交叉验证
from sklearn.metrics import roc_auc_score,roc_curve,auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
from sklearn.inspection import permutation_importance
import os
from sklearn.datasets import fetch_california_housing
from sklearn.datasets import fetch_openml
import openpyxl

os.getcwd()
all=pd.read_csv("Protein.DEP.top50.csv",index_col=0)
allData=all.replace(["Subgroup1","Subgroup2","Subgroup3"],["1","2","3"])
random.seed(300)
num=random.sample(range(0, 145), 102)
train=allData.iloc[num,]
test=allData.iloc[list(set(range(0,145,1))-set(num)),]
X_train=train.iloc[:,0:-1]
Y_train=train.iloc[:,-1]
X_test=test.iloc[:,0:-1]
Y_test=test.iloc[:,-1]
###------------------------ 1选择最优参数------------------------###
# optimal parameter selection process
# Define the parameter range for the Random Forest
n_estimators=list(range(0,101,10))
n_estimators[0]=1
samples_leaf=list(range(0, 21, 1))
samples_leaf.pop(0)
param_test1 = {"n_estimators":n_estimators,
               'min_samples_leaf':samples_leaf}
# 10-fold cross-validation
gsearch1 = GridSearchCV(estimator=RandomForestClassifier(random_state=50,n_jobs=-1),param_grid=param_test1,cv=10)
gsearch1.fit(X_train,Y_train)
# print result
print(gsearch1.best_params_)
# {'min_samples_leaf': 3, 'n_estimators': 60} # final set: n_estimators=60, min_samples_leaf=3
print("best accuracy:%f" % gsearch1.best_score_)
# best accuracy:0.842727
param1Result=pd.DataFrame(gsearch1.cv_results_["params"])
param1Result["mean_accuracy"]=gsearch1.cv_results_['mean_test_score']
param1Result.to_csv("Protein.RandomForest.parameter.test.csv")
###------------------------ 2建模 获取重要性排序 ------------------------###
# 训练数据
model = RandomForestClassifier(n_estimators=60,min_samples_leaf=3,random_state=50,oob_score=False)
model.fit(X_train,Y_train)
print(model.score(test.iloc[:,0:-1],test.iloc[:,-1])) #准确率
# print(model.oob_score_)
# 每个特征重要性
importances=pd.DataFrame(model.feature_importances_,columns=["Importance"])
importances.insert(0,"GeneName",allData.columns.tolist()[0:-1])
importances.sort_values(by=["Importance"],axis=0,ascending=False,inplace=True)
importances.index=range(0,importances.shape[0])
importances.to_csv("Protein.RandomForest.importance.csv")
########累积auc##############
def getAUC(X_train, Y_train, X_test, Y_test,n_estimators, min_samples_leaf):
    model = RandomForestClassifier(n_estimators=n_estimators, min_samples_leaf=min_samples_leaf, random_state=50, oob_score=False)
    model.fit(X_train, Y_train)
    # Binarize the output
    classes = [1, 2, 3]
    Y_test_binar = label_binarize(Y_test, classes=[1, 2, 3])
    Y_pred_binar = label_binarize(model.predict(X_test), classes=[1, 2, 3])
    # Compute micro-average ROC curve and ROC area
    fpr, tpr, thresholds = roc_curve(Y_test_binar.ravel(), Y_pred_binar.ravel())
    roc_auc = auc(fpr, tpr)
    return (roc_auc)
#################################
# top50
AUCs = []
for top in range(1, 51):
    select_features = importances.iloc[0:top,0]
    top_auc = getAUC(X_train.loc[:,select_features.tolist()],Y_train,
                     X_test.loc[:,select_features.tolist()],Y_test,
                     60,3)
    AUCs.append(top_auc)

cumulate_auc=pd.DataFrame(AUCs,columns=["AUC"])
cumulate_auc.insert(0,"GeneName",importances.iloc[0:49,0])
cumulate_auc.to_csv("Protein.RandomForest.Cumulate_auc_top50.csv")
###################################
# top39
AUCs = []
for top in range(1, 40):
    select_features = importances.iloc[0:top,0]
    top_auc = getAUC(X_train.loc[:,select_features.tolist()],Y_train,
                     X_test.loc[:,select_features.tolist()],Y_test,
                     60,3)
    AUCs.append(top_auc)

cumulate_auc=pd.DataFrame(AUCs,columns=["AUC"])
cumulate_auc.insert(0,"GeneName",importances.iloc[0:39,0])
cumulate_auc.to_csv("Protein.RandomForest.Cumulate_auc_top39.csv")
###------------------------ 6 训练集测试集选择特征后的AUC曲线 ------------------------###
candidate=importances.iloc[0:39,:]
candidateGene=candidate["GeneName"].tolist()
selectData=train.loc[:,candidateGene]
model_reduce = RandomForestClassifier(n_estimators=60,min_samples_leaf=3,random_state=50,oob_score=False)
model_reduce.fit(selectData,Y_train)

classes=[1,2,3]
##train
Y_train_reduce_binar = label_binarize(Y_train, classes=[1,2,3])
Y_train_reduce_pred_binar=label_binarize(model_reduce.predict(selectData),classes=[1,2,3])

fpr_train = dict()
tpr_train = dict()
roc_auc_train = dict()
for i in range(len(classes)):
    fpr_train[i], tpr_train[i], thresholds = roc_curve(Y_train_reduce_binar[:, i],Y_train_reduce_pred_binar[:, i])
    roc_auc_train[i] = auc(fpr_train[i], tpr_train[i])

# micro-average ROC curve（方法一）
fpr_train["micro"], tpr_train["micro"], thresholds = roc_curve(Y_train_reduce_binar.ravel(),Y_train_reduce_pred_binar.ravel())
roc_auc_train["micro"] = auc(fpr_train["micro"], tpr_train["micro"])

##test
Y_test_reduce_binar = label_binarize(Y_test, classes=[1,2,3])
Y_test_reduce_pred_binar=label_binarize(model_reduce.predict(X_test.loc[:,candidateGene]),classes=[1,2,3])

fpr_test = dict()
tpr_test = dict()
roc_auc_test = dict()
for i in range(len(classes)):
    fpr_test[i], tpr_test[i], thresholds = roc_curve(Y_test_reduce_binar[:, i],Y_test_reduce_pred_binar[:, i])
    roc_auc_test[i] = auc(fpr_test[i], tpr_test[i])

# micro-average ROC curve（方法一）
fpr_test["micro"], tpr_test["micro"], thresholds = roc_curve(Y_test_reduce_binar.ravel(),Y_test_reduce_pred_binar.ravel())
roc_auc_test["micro"] = auc(fpr_test["micro"], tpr_test["micro"])

# 输出结果文件
Y_train_data_micro = pd.merge(pd.DataFrame(Y_train_reduce_binar.ravel()), pd.DataFrame(Y_train_reduce_pred_binar.ravel()), left_index=True, right_index=True)
Y_train_data_micro.columns = ["micro_2Category","micro"]
Y_train_data_micro.to_csv("Y_train_reduce_binar_ravel_label.top39.csv")

Y_test_data_micro = pd.merge(pd.DataFrame(Y_test_reduce_binar.ravel()), pd.DataFrame(Y_test_reduce_pred_binar.ravel()), left_index=True, right_index=True)
Y_test_data_micro.columns = ["micro_2Category","micro"]
Y_test_data_micro.to_csv("Y_test_reduce_binar_ravel_label.top39.csv")

#画图部分
plt.figure()
plt.plot(fpr_train["micro"], tpr_train["micro"],color="black",
         label='Reduce train ROC curve (AUC = {0:0.4f})'
               ''.format(roc_auc_train["micro"]),linewidth=3)
plt.plot(fpr_test["micro"], tpr_test["micro"],color="darkgrey",
         label='Reduce test ROC curve (AUC = {0:0.4f})'
               ''.format(roc_auc_test["micro"]),linewidth=3)

plt.plot([0,1],[0,1],'r--')
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc="lower right")
plt.grid(linestyle='-.')
plt.grid(True)
plt.savefig("Protein.RandomForest.testANDtrain_AUC curves.top39.png", dpi=300)
plt.show()
# ###############利用累积auc选择的结果建模画图#####################
# candidateGene_auc=pd.read_csv("top39.list",header = None)
# candidateGene1=candidateGene_auc.loc[:,0].tolist()
candidateGene_auc=pd.read_csv("Protein.RandomForest.Cumulate_auc_top39.csv")
candidate=candidateGene_auc
candidateGene1=candidate["GeneName"].tolist()
selectData_train=X_train.loc[:,candidateGene1]
selectData_test=X_test.loc[:,candidateGene1]
# model
model_reduce1 = RandomForestClassifier(n_estimators=60,min_samples_leaf=3,random_state=50,oob_score=False)
model_reduce1.fit(selectData_train,Y_train)
#print(model_reduce1.best_params_)
classes=[1,2,3]
#train
Y_train_reduce_binar = label_binarize(Y_train, classes=[1,2,3])
Y_train_reduce_pred_binar=model_reduce1.predict_proba(selectData_train) #不需要label_binarize

fpr_train = dict()
tpr_train = dict()
roc_auc_train = dict()
for i in range(len(classes)):
    fpr_train[i], tpr_train[i], thresholds = roc_curve(Y_train_reduce_binar[:, i],Y_train_reduce_pred_binar[:, i])
    roc_auc_train[i] = auc(fpr_train[i], tpr_train[i])

# micro-average ROC curve（方法一）
fpr_train["micro"], tpr_train["micro"], thresholds = roc_curve(Y_train_reduce_binar.ravel(),Y_train_reduce_pred_binar.ravel())
roc_auc_train["micro"] = auc(fpr_train["micro"], tpr_train["micro"])

##test
Y_test_reduce_binar = label_binarize(Y_test, classes=[1,2,3])
Y_test_reduce_pred_binar=model_reduce1.predict_proba(selectData_test)  #不需要label_binarize

fpr_test = dict()
tpr_test = dict()
roc_auc_test = dict()
for i in range(len(classes)):
    fpr_test[i], tpr_test[i], thresholds = roc_curve(Y_test_reduce_binar[:, i],Y_test_reduce_pred_binar[:, i])
    roc_auc_test[i] = auc(fpr_test[i], tpr_test[i])

# micro-average ROC curve（方法一）
fpr_test["micro"], tpr_test["micro"], thresholds = roc_curve(Y_test_reduce_binar.ravel(),Y_test_reduce_pred_binar.ravel())
roc_auc_test["micro"] = auc(fpr_test["micro"], tpr_test["micro"])

# macro-average ROC curve（方法二）
all_fpr = np.unique(np.concatenate([fpr_test[i] for i in range(len(classes))]))
mean_tpr = np.zeros_like(all_fpr)
for i in range(len(classes)):
    mean_tpr += np.interp(all_fpr, fpr_test[i], tpr_test[i])

mean_tpr /= len(classes)
mean_tpr[0] = 0.0

fpr_test["macro"] = all_fpr
tpr_test["macro"] = mean_tpr
roc_auc_test["macro"] = auc(fpr_test["macro"], tpr_test["macro"])

# Save reuslt
auc_df = pd.DataFrame(roc_auc_test, index=["AUCs"]).T
auc_df["Label"] = auc_df.index.values
auc_df.to_excel("AUC.top39.xlsx")

rlt = pd.DataFrame()
for cls in tpr_test.keys():
    df = pd.DataFrame(
        {
            "TPRs": tpr_test[cls],
            "FPRs": fpr_test[cls],
            "Label": np.array([cls] * len(tpr_test[cls]))
        }
    )
    rlt = pd.concat([rlt, df])

rlt.to_excel("TPR_FPR_Label.top39.xlsx")
