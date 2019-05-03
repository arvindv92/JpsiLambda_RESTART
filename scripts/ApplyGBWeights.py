import sys
import numpy
import root_numpy
import pandas
import cPickle as pickle
from hep_ml import reweight
from sklearn.model_selection import train_test_split

# mcOpt = 1 for JpsiLambda
# mcOpt = 2 for JpsiSigma
# mcOpt = 3 for JpsiXi (reco'd JpsiLambda)
# mcOpt = 4 for JpsiXi (reco'd JpsiXi)

run = int(sys.argv[1])
mcOpt = int(sys.argv[2])
isGen = int(sys.argv[3])  # set to 1 to run over generator, 0 to run over rec.

print 'Processing Run', run, 'mcOpt ', mcOpt
if isGen:
    print 'Generator\n'
else:
    print 'Reconstructed\n'
columns = ['Lb_P', 'Lb_PT', 'Lb_ETA', 'Jpsi_P', 'Jpsi_PT', 'Jpsi_ETA',
           'L_P', 'L_PT', 'L_ETA', 'p_P', 'p_PT', 'p_ETA', 'pi_P', 'pi_PT',
           'pi_ETA']


mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'.format(run)
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
# elif mcOpt == 4:
#     filePath = '../rootFiles/mcFiles/JpsiXi/run{}/'
#     # fileName = 'jpsixi.root'
#     fileName = 'jpsilambda_aliased.root'
#     # treeName = 'Xib2JpsiXiTree/MyTuple'
#     treeName = 'MCDecayTree'

if isGen:
    fileName = part + '.root'
    treeName = 'MCTuple/MCDecayTree'
else:
    fileName = part + '_pidgen.root'
    treeName = 'MyTuple'
filePath = filePath.format(run)

# Get the reweighter
with open(mcPath + 'gb_wts.pkl') as f:
    reweighter = pickle.load(f)

    # Get the input file
    original = root_numpy.root2array(filePath + fileName,
                                     treeName, branches=columns)
    original = pandas.DataFrame(original, dtype=float)

    # Predict weights
    gb_weights = reweighter.predict_weights(original)
    gb_weights.dtype = [('gb_wts', 'float64')]
    if isGen:
        # Write out weights to separate ROOT file
        root_numpy.array2root(gb_weights,
                              filePath + 'RW/gbWeights_gen.root',
                              treename='MyTuple', mode='recreate')
    else:
        root_numpy.array2root(gb_weights,
                              filePath + fileName,
                              treename='MyTuple', mode='update')
