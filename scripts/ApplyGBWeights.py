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

print 'Processing Run', run, 'mcOpt ', mcOpt
columns = ['Lb_P', 'Lb_PT', 'Lb_ETA', 'Jpsi_P', 'Jpsi_PT', 'Jpsi_ETA',
           'L_P', 'L_PT', 'L_ETA', 'p_P', 'p_PT', 'p_ETA', 'pi_P', 'pi_PT',
           'pi_ETA', 'p_ProbNNp', 'pi_ProbNNpi', 'p_PIDp']


mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'.format(run)
if mcOpt == 1:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'
    fileName = 'jpsilambda.root'
    # fileName = 'jpsilambda_aliased.root'
    treeName = 'Lb2JpsiLTree/MyTuple'
    # treeName = 'MCDecayTree'
elif mcOpt == 2:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run{}/'
    fileName = 'jpsisigma.root'
    # fileName = 'jpsisigma_aliased.root'
    treeName = 'Lb2JpsiLTree/MyTuple'
    # treeName = 'MCDecayTree'
elif mcOpt == 3:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiXi/run{}/'
    fileName = 'jpsixi.root'
    # fileName = 'jpsixi_aliased.root'
    treeName = 'Lb2JpsiLTree/MyTuple'
    # treeName = 'MCDecayTree'
elif mcOpt == 4:
    filePath = '../rootFiles/mcFiles/JpsiXi/run{}/'
    # fileName = 'jpsixi.root'
    fileName = 'jpsilambda_aliased.root'
    # treeName = 'Xib2JpsiXiTree/MyTuple'
    treeName = 'MCDecayTree'

filePath = filePath.format(run)

# Get the reweighter
with open(mcPath + 'gb_wts_pid.pkl') as f:
    reweighter = pickle.load(f)

    # Get the input file
    original = root_numpy.root2array(filePath + fileName,
                                     treeName, branches=columns)
    original = pandas.DataFrame(original, dtype=float)

    # Predict weights
    gb_weights = reweighter.predict_weights(original)
    gb_weights.dtype = [('gb_wts_pid', 'float64')]
    # Write out weights to separate ROOT file
    root_numpy.array2root(gb_weights,
                          filePath + 'RW/gbWeights_pid_rec.root',
                          treename='MyTuple', mode='recreate')
