---
title: "交叉验证和特征选择"
date: 2020-02-24 04:21:00 +0800
categories: [Machine Learning, Feature Selection]
tags: [Feature Selection, Cross Validation, Prediction Variance]
---

交叉验证通常用于模型选择，那么如何使用交叉验证进行特征选择呢？

# 交叉验证

模型在训练的时候往往是**高方差**估计. 从bias-variance tradeoff的角度看， 在有限数据的情况下训练的结果往往是bias较小单varance很大从而造成模型过拟合，因此降低了模型的泛化能力。为了在有限数据集的条件下使得模型尽可能满足泛化能力，提出了交叉验证 (Cross  Validation, CV)，目的是增加模型训练的bias，降低variance, 从而在bias-variance tradeoff的角度上达到一定的平衡，提高模型的泛化能力。换句话说，交叉验证可以表征不同数据集的variance。

## K-fold 交叉验证

交叉验证本质上是对给定的训练数据集进行切分。比如5-fold(即K=5)交叉验证的处理如下：首先将整个数据集等分为5份， 然后在一次轮迭代中使用第一等分用于模型训练过程中的验证集(validation set)， 剩下的4等份数据用于模型训练；第二轮迭代中使用第二等份作为验证集，剩下的用于模型训练，如此迭代5轮结束。 

![5-fold cross validation]({{ "/assets/img/blogs/5fold-cross-validation.png" | relative_url }}) 图一, 5-fold cross validation, image from [here](http://ethen8181.github.io/machine-learning/model_selection/model_selection.html)

例如用于模型最优参数选择（即模型选择）的伪代码如下：

```python
# optional parameters list for a given type of model, eg, random forest
param_arr = ...

# splitting origianl whole training data
val_index, train_index = split_whole_training_data(whole_training_dataset)

# val_index and train_index are 5 element of array, and each element is also an array containing the index of i-fold validation dataset and training dataset, respectively.

pref_metrix=[]

for i_param in param_arr:

    i_param_perform = []
    for i_val_idx, i_train_idx in zip(val_index, train_index):
        i_val_data = whole_training_dataset[i_val_idx]
        i_train_data = whole_training_dataset[i_train_idx]

        # training model
        i_model = train(i_train_data)

        # evaluate model, like AUC, F1 and so on.
        preform_metrics = evaluate(i_model, i_val_data)
        i_param_perform.append(preform_metrics)
    
    # calculate mean value
    pref_metrix.append(mean(i_param_perform))

# finally, select the parameter according to the max pref_metrix.

```

以上是基于交叉验证的模型选择（也就是选择模型针对给定数据集的最有参数）的通用步(tao)骤(lu)。既然通过数据集切分的交叉验证表征了不同数据集的variance，那么就可以在这些切分后的数据集中找到一个在整体上使得模型预测的validation最小的“最优模型”，从而提高模型泛化能力。

# 使用交叉验证进行特征选择

那么如何通过交叉验证技术进行特征选择呢？它的过程是否与通过交叉验证技术进行模型选择一样呢？如果我们先在整个数据集上进行特征选择，然后再利用交叉验证进行模型选择，那么这样的过程是否是正确呢？是否人为的引入了偏差(bias)? 在回答这些问题前，我们先了解什么是数据泄露(data leakage)。



在 [stackexchange](https://stats.stackexchange.com/questions/27750/feature-selection-and-cross-validation) 有人提问：

> ..., surely if you select a feature set and then evaluate your model using *only* the selected features using cross validation, then you are getting an unbiased estimate of generalized model performance on those features (this assumes the sample under study are representive of the populatation)?

> With this procedure one cannot of course claim an optimal feature set but can one report the performance of the selected feature set on unseen data as valid? ...




# DEMO

https://zhuanlan.zhihu.com/p/74198735