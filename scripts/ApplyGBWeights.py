# This code applies weights to MC to correct data-MC differences
# The RW is applied to both reconstructed and generated MC
# The RW is applied to the pidgen MC files (i.e. before any selections are applied on it)

import sys
# import numpy
import root_numpy
import pandas
import cPickle as pickle
# from hep_ml import reweight
# from sklearn.model_selection import train_test_split

# mcType = 1 for Lb -> J/psi Lambda MC
# mcType = 2 for Lb -> J/psi Sigma MC        (reco'd JpsiLambda)
# mcType = 3 for Xib -> Jpsi Xi MC           (reco'd JpsiLambda)
# mcType = 8 for Xib0 -> J/psi Lambda MC
# mcType = 9 for Xib- -> J/psi Xi- MC        (reco'd J/psi Xi-)

run = int(sys.argv[1])
mcType = int(sys.argv[2])
isGen = int(sys.argv[3])  # set to 1 to run over generator, 0 to run over rec.

list = [1, 2, 3, 8, 9]
if mcType not in list:
    print 'mcType not in allowed list. Exiting!\n'
    sys.exit()

print 'Processing Run', run, 'mcType ', mcType
if isGen:
    print 'Generator\n'
else:
    print 'Reconstructed\n'

if mcType != 3 and mcType != 9:
    if not isGen:
        columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT', 'L_P',
                   'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT', 'p_ProbNNghost', 'pi_ProbNNghost']
    else:
        columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT', 'L_P',
                   'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT']
else:
    columns = ['Lb_P', 'Lb_PT', '(p_P-pi_P)/(p_P+pi_P)']
# columns = ['Lb_P', 'Lb_PT', 'Jpsi_P', 'Jpsi_PT',
#            'L_P', 'L_PT', 'p_P', 'p_PT', 'pi_P', 'pi_PT']


# mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'.format(run)
mcPath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/'
if mcType == 1:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run{}/'
    part = 'jpsilambda'
elif mcType == 2:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run{}/'
    part = 'jpsisigma'
elif mcType == 3:
    filePath = '../rootFiles/mcFiles/JpsiLambda/JpsiXi/run{}/'
    part = 'jpsixi'
elif mcType == 8:
    filePath = '../rootFiles/mcFiles/JpsiLambda/Xib0/run{}/'
    part = 'xib0'
elif mcType == 9:
    filePath = '../rootFiles/mcFiles/JpsiXi/run{}/'
    part = 'jpsixi'

if isGen:
    fileName = part + '.root'
    treeName = 'MCTuple/MCDecayTree'
else:
    treeName = 'MyTuple'
    if mcType != 9:
        fileName = part + '_pidgen.root'
    else:
        fileName = part + '_cut_LL.root'
filePath = filePath.format(run)

if mcType != 3 and mcType != 9:
    if not isGen:
        wtsFileName = 'gb_wts.pkl'
    else:
        wtsFileName = 'gb_wts_gen.pkl'
else:
    wtsFileName = 'gb_wts_xib.pkl'
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
