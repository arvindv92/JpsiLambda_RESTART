import sys
import numpy
import root_numpy
import pandas
import cPickle as pickle
from hep_ml import reweight
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from hep_ml.metrics_utils import ks_2samp_weighted
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

run = int(sys.argv[1])

columns = ['Lb_P', 'Lb_PT', 'Lb_ETA', 'Jpsi_P', 'Jpsi_PT', 'Jpsi_ETA', 'L_P',
           'L_PT', 'L_ETA', 'p_P', 'p_PT', 'p_ETA', 'pi_P', 'pi_PT', 'pi_ETA',
           'ntracks', 'SW']

mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'.format(run)
dataPath = '../rootFiles/dataFiles/JpsiLambda/run{}/sWeightSanity/'.format(run)


original = root_numpy.root2array(mcPath + 'jpsilambda_withsw.root',
                                 'Lb2JpsiLTree/MyTuple', branches=columns,
                                 selection='(Lb_BKGCAT==0||Lb_BKGCAT==50)')
original_noTM = root_numpy.root2array(mcPath + 'jpsilambda_withsw.root',
                                      'Lb2JpsiLTree/MyTuple', branches=columns)
target = root_numpy.root2array(dataPath + 'jpsilambda_LL_sanity_withsw_noPID.root',
                               'MyTuple', branches=columns)

original = pandas.DataFrame(original, dtype=float)
original_noTM = pandas.DataFrame(original_noTM, dtype=float)
target = pandas.DataFrame(target, dtype=float)

original_weights = numpy.array(original['SW'])
target_weights = numpy.array(target['SW'])

original_train, original_test = train_test_split(original, random_state=42)
target_train, target_test = train_test_split(target, random_state=42)

original_weights_train = original_train['SW']
original_weights_test = original_test['SW']

target_weights_train = target_train['SW']
target_weights_test = target_test['SW']

hist_settings = {'bins': 50, 'alpha': 0.7, 'density': True}


def draw_distributions(myoriginal, mytarget, new_original_weights, targetwts):
    sum_ks = 0
    ctr = 0
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[0:6], 1):
        ctr = ctr + 1
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]),
                                [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim,
                 **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim,
                 **hist_settings)
        plt.title(column)
        myks = ks_2samp_weighted(myoriginal[column], mytarget[column],
                                 weights1=new_original_weights,
                                 weights2=targetwts)
        sum_ks = sum_ks + myks
        # print('KS over ', column, ' = ', myks)
    plt.draw()
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[6:12], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]),
                                [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim,
                 **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim,
                 **hist_settings)
        plt.title(column)
        myks = ks_2samp_weighted(myoriginal[column], mytarget[column],
                                 weights1=new_original_weights,
                                 weights2=targetwts)
        sum_ks = sum_ks + myks
        # print('KS over ', column, ' = ', myks)
    plt.draw()
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(columns[12:16], 1):
        xlim = numpy.percentile(numpy.hstack([mytarget[column]]),
                                [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(myoriginal[column], weights=new_original_weights, range=xlim,
                 **hist_settings)
        plt.hist(mytarget[column], weights=targetwts, range=xlim,
                 **hist_settings)
        plt.title(column)
        myks = ks_2samp_weighted(myoriginal[column], mytarget[column],
                                 weights1=new_original_weights,
                                 weights2=targetwts)
        sum_ks = sum_ks + myks
        # print('KS over ', column, ' = ', myks)
    plt.draw()
    avg_ks = sum_ks / ctr
    print('average of KS distances = ', avg_ks)
    return avg_ks


print 'Original'
avgks_orig = draw_distributions(original.iloc[:, :-1], target.iloc[:, :-1],
                                original_weights, target_weights)

# **********Binned Re-weighting****************
# bins_reweighter = reweight.BinsReweighter(n_bins=30, n_neighs=2.)
# bins_reweighter.fit(original_train.iloc[:,:-1], target_train.iloc[:,:-1],
# original_weights_train, target_weights_train)
#
# bins_weights_test = bins_reweighter.predict_weights(original_test.iloc[:,:-1])
# # validate reweighting rule on the test part comparing 1d projections
#
# print 'After binned re-weighting'
# draw_distributions(original_test.iloc[:,:-1], target_test.iloc[:,:-1],
# bins_weights_test, target_weights_test)
# **********************************************

# *********Gradient Boosted Re-weighting********


# This is currently the best config for Run1
# reweighter = reweight.GBReweighter(n_estimators=200, learning_rate=0.1,
#                                         max_depth=3, min_samples_leaf=50,
#                                         gb_args={'subsample': 0.2,
#                                                  'random_state': 42})

reweighter = reweight.GBReweighter(n_estimators=100, learning_rate=0.1,
                                   max_depth=3, min_samples_leaf=50,
                                   gb_args={'subsample': 0.5,
                                            'random_state': 42})

# reweighter.fit(original_train.iloc[:, :-1], target_train.iloc[:, :-1],
#                original_weights_train, target_weights_train)
reweighter.fit(original.iloc[:, :-1], target.iloc[:, :-1],
               original_weights, target_weights)

# gb_weights_test = reweighter.predict_weights(original_test.iloc[:, :-1])
gb_weights = reweighter.predict_weights(original.iloc[:, :-1])

gb_weights_noTM = reweighter.predict_weights(original_noTM.iloc[:, :-1])

gb_weights_noTM.dtype = [('gb_wts', 'float64')]
# validate reweighting rule on the test part comparing 1d projections
# print 'After GB reweighting on test sample'
# draw_distributions(original_test.iloc[:, :-1], target_test.iloc[:, :-1],
#                    gb_weights_test, target_weights_test)
print 'After GB reweighting on all'
avgks_rw = draw_distributions(original.iloc[:, :-1], target.iloc[:, :-1],
                              gb_weights, target_weights)

root_numpy.array2root(gb_weights_noTM,
                      mcPath + 'jpsilambda_weighted.root',
                      treename='MyTuple', mode='recreate')

# *******Exporting RW formula for re-use********
with open(mcPath + 'gb_wts.pkl', 'w') as f:
    pickle.dump(reweighter, f)
# **********************************************


# **********************************************
#
# data = numpy.concatenate([original_test.iloc[:, :-1],
#                           target_test.iloc[:, :-1]])
# labels = numpy.array([0] * len(original_test.iloc[:, :-1])
#                      + [1] * len(target_test.iloc[:, :-1]))
#
# weights = {}
# weights['original'] = original_weights_test
# # weights['bins'] = bins_weights_test
# weights['gb_weights'] = gb_weights_test
#
# for name, new_weights in weights.items():
#     W = numpy.concatenate([new_weights / new_weights.sum() * len(target_test.iloc[:, :-1]),
#                            target_weights_test / target_weights_test.sum() * len(target_test.iloc[:, :-1])])
#     Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data, labels, W,
#                                                     random_state=42, train_size=0.51)
#     clf = GradientBoostingClassifier(subsample=0.3, learning_rate=0.1,
#                                      n_estimators=1000).fit(Xtr, Ytr,
#                                                             sample_weight=Wtr)
#
#     print(name, roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1],
#                               sample_weight=Wts))
