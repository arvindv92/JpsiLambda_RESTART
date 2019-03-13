import numpy
import root_numpy
print 'poop1'
import pandas
from hep_ml import reweight
print 'poop2'
from matplotlib import pyplot as plt

print 'poop'
columns = ['log10(Lb_ConsLb_chi2)', 'p_PIDp', 'L_dm', 'log10(Lb_MINIPCHI2)', 'log10(acos(Lb_DIRA_OWNPV))',
           'log10(Lb_FD_OWNPV)', 'log10(Jpsi_MINIPCHI2)', 'log10(Jpsi_M)', 
           'log10(p_MINIPCHI2)', 'p_PT', 'p_ProbNNp', 'log10(L_FDCHI2_ORIVX)', 
           'log10(L_FD_ORIVX)', 'log10(acos(L_DIRA_OWNPV))', 'log10(L_MINIPCHI2)', 
           'log10(pi_MINIPCHI2)', 'log10(pi_PT)', 'Lb_PT', 'Lb_P', 'Lb_ETA', 'SW']#, 

#columns = ['p_PIDp', 'p_ProbNNp', 'log10(pi_PT)', 'SW']
#           'p_ProbNNp', 'L_FDCHI2_ORIVX']#'L_DIRA_ORIVX', 'L_FD_ORIVX', 'L_DIRA_OWNPV']#,
#           'L_MINIPCHI2', 'pi_MINIPCHI2', 'pi_PT', 'SW']

original = root_numpy.root2array('../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_sanity_LL_withsw.root', 'MyTuple', branches=columns, selection='(Lb_BKGCAT==0||Lb_BKGCAT==50)')
target = root_numpy.root2array('../rootFiles/dataFiles/JpsiLambda/run1/sWeightSanity/jpsilambda_LL_sanity_withsw.root', 'MyTuple', branches=columns)

original = pandas.DataFrame(original, dtype=float)
target = pandas.DataFrame(target, dtype=float)

#original_weights = numpy.ones(len(original))
original_weights = numpy.array(original['SW'])
target_weights = numpy.array(target['SW'])

from sklearn.model_selection import train_test_split
original_train, original_test = train_test_split(original, random_state=42)
target_train, target_test = train_test_split(target, random_state=42)

original_weights_train = original_train['SW']
original_weights_test = original_test['SW']

target_weights_train = target_train['SW']
target_weights_test = target_test['SW']

from hep_ml.metrics_utils import ks_2samp_weighted
hist_settings = {'bins': 50, 'alpha': 0.7, 'density': True}

def draw_distributions(myoriginal, mytarget, new_original_weights, targetwts):
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[0:6], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]), [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim, **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim, **hist_settings)
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(myoriginal[column], mytarget[column],
                                                           weights1=new_original_weights, weights2=targetwts))
    plt.draw()
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[6:12], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]), [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim, **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim, **hist_settings)
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(myoriginal[column], mytarget[column],
                                                           weights1=new_original_weights, weights2=targetwts))
    plt.draw()
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[12:18], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]), [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim, **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim, **hist_settings)
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(myoriginal[column], mytarget[column],
                                                           weights1=new_original_weights, weights2=targetwts))
    plt.draw()
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[18:20], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]), [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim, **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim, **hist_settings)
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(myoriginal[column], mytarget[column],
                                                           weights1=new_original_weights, weights2=targetwts))
    plt.draw()

print 'Original'
draw_distributions(original.iloc[:,:-1],target.iloc[:,:-1], original_weights, target_weights)

#**********Binned Re-weighting****************
# bins_reweighter = reweight.BinsReweighter(n_bins=30, n_neighs=2.)
# bins_reweighter.fit(original_train.iloc[:,:-1], target_train.iloc[:,:-1], original_weights_train, target_weights_train)

# bins_weights_test = bins_reweighter.predict_weights(original_test.iloc[:,:-1])
# # validate reweighting rule on the test part comparing 1d projections

# print 'After binned re-weighting'
# draw_distributions(original_test.iloc[:,:-1], target_test.iloc[:,:-1], bins_weights_test, target_weights_test)
#**********************************************

#*********Gradient Boosted Re-weighting********
reweighter = reweight.GBReweighter(n_estimators=500, learning_rate=0.05, max_depth=3, min_samples_leaf=50, 
                                   gb_args={'subsample': 0.3, 'random_state':42})
reweighter.fit(original_train.iloc[:,:-1], target_train.iloc[:,:-1], original_weights_train, target_weights_train)

gb_weights_test = reweighter.predict_weights(original_test.iloc[:,:-1])
gb_weights = reweighter.predict_weights(original.iloc[:,:-1])
# validate reweighting rule on the test part comparing 1d projections
print 'After GB reweighting'
draw_distributions(original_test.iloc[:,:-1], target_test.iloc[:,:-1], gb_weights_test, target_weights_test)
#**********************************************

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

data = numpy.concatenate([original_test.iloc[:,:-1], target_test.iloc[:,:-1]])
labels = numpy.array([0] * len(original_test.iloc[:,:-1]) + [1] * len(target_test.iloc[:,:-1]))

weights = {}
weights['original'] = original_weights_test
#weights['bins'] = bins_weights_test
weights['gb_weights'] = gb_weights_test

for name, new_weights in weights.items():
    W = numpy.concatenate([new_weights / new_weights.sum() * len(target_test.iloc[:,:-1]), target_weights_test / target_weights_test.sum() * len(target_test.iloc[:,:-1])])
    Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data, labels, W, random_state=42, train_size=0.51)
    clf = GradientBoostingClassifier(subsample=0.3, learning_rate = 0.1, n_estimators=1000).fit(Xtr, Ytr, sample_weight=Wtr)
    
    print(name, roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts))

