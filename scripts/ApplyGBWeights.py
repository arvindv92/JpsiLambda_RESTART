import sys
# import numpy
import root_numpy
import pandas
import cPickle as pickle
# from hep_ml import reweight
# from sklearn.model_selection import train_test_split

# mcOpt = 1 for JpsiLambda
# mcOpt = 2 for JpsiSigma
# mcOpt = 3 for JpsiXi (reco'd JpsiLambda)
# mcOpt = 10 for JpsiXi (reco'd JpsiXi)
# mcOpt = 11 for Xib0->J/psi Lambda
run = int(sys.argv[1])
mcOpt = int(sys.argv[2])
isGen = int(sys.argv[3])  # set to 1 to run over generator, 0 to run over rec.

print 'Processing Run', run, 'mcOpt ', mcOpt
if isGen:
    print 'Generator\n'
else:
    print 'Reconstructed\n'

if not isGen:
    columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT', 'L_P',
               'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT', 'p_ProbNNghost', 'pi_ProbNNghost']
else:
    columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT', 'L_P',
               'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT']
# columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT',
#            'L_P', 'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT']


# mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'.format(run)
mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/'
if mcOpt == 1:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'
    part = 'jpsilambda'
elif mcOpt == 2:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run{}/'
    part = 'jpsisigma'
elif mcOpt == 3:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiXi/run{}/'
    part = 'jpsixi'
elif mcOpt == 6:
    filePath = '../rootFiles/mcFiles/JpsiLambda/Lst1405/run{}/'
    part = 'lst1405'
elif mcOpt == 7:
    filePath = '../rootFiles/mcFiles/JpsiLambda/Lst1520/run{}/'
    part = 'lst1520'
elif mcOpt == 8:
    filePath = '../rootFiles/mcFiles/JpsiLambda/Lst1600/run{}/'
    part = 'lst1600'
elif mcOpt == 9:
    filePath = '../rootFiles/mcFiles/JpsiLambda/chiC1/run{}/'
    part = 'chic1'
elif mcOpt == 10:
    filePath = '../rootFiles/mcFiles/JpsiXi/run{}/'
    part = 'jpsixi'
elif mcOpt == 11:
    filePath = '../rootFiles/mcFiles/JpsiLambda/Xib0/run{}/'
    part = 'xib0'
if isGen:
    fileName = part + '.root'
    treeName = 'MCTuple/MCDecayTree'
else:
    treeName = 'MyTuple'
    if(mcOpt != 10):
        fileName = part + '_pidgen.root'
    else:
        fileName = part + '_cut_LL.root'
filePath = filePath.format(run)

if not isGen:
    wtsFileName = 'gb_wts.pkl'
else:
    wtsFileName = 'gb_wts_gen.pkl'
# Get the reweighter
with open(mcPath + wtsFileName) as f:
    reweighter = pickle.load(f)

    # Get the input file
    original = root_numpy.root2array(filePath + fileName,
                                     treeName, branches=columns)
    original = pandas.DataFrame(original, dtype=float)

    # Predict weights
    gb_weights = reweighter.predict_weights(original)
    gb_weights.dtype = [('GB_WT', 'float64')]
    if isGen:
        # Write out weights to separate ROOT file
        root_numpy.array2root(gb_weights,
                              filePath + 'RW/gbWeights_gen.root',
                              treename='MyTuple', mode='recreate')
    else:
        root_numpy.array2root(gb_weights,
                              filePath + fileName,
                              treename='MyTuple', mode='update')
