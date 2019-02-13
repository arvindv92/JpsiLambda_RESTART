import glob
import os

files = glob.glob("../logs/data/JpsiLambda/run2/OptimizeFinalBDT*.txt")
ctr = 0
myFOM_nonZero = 0.0
myFOM_Zero = 0.0
maxFOM_nonZero = -1.0
maxFOM_Zero = -1.0
bdtCut_Zero = -1.0
bdtCut_nonZero = -1.0
bdtConf_best_nonZero = -1
isoConf_best_nonZero = -1
isoVersion_best_nonZero = -1
newFlag_best_nonZero = False
bdtCut_best_nonZero = -1.0
bdtConf_best_Zero = -1
isoConf_best_Zero = -1
isoVersion_best_Zero = -1
newFlag_best_Zero = False
bdtCut_best_Zero = -1.0

for fName in files:
    bdtConf = int(fName.split('BDT')[1][0])
    isoConf = int(fName.split('iso')[1][0])
    isoVersion = int(fName.split('_v')[1][0])
    if 'new' in fName:
        newFlag = True
    with open(fName) as f:
        ctr = 0
        lines = f.readlines()
        for line in lines:
            if 'MAXIMUM FOM' in line:
                if ctr == 0:
                    myFOM_nonZero = float(line.split()[3])
                    bdtCut_nonZero = float(line.split()[7])
                if ctr == 1:
                    myFOM_Zero = float(line.split()[3])
                    bdtCut_Zero = float(line.split()[7])
                ctr = ctr + 1
    if myFOM_nonZero > maxFOM_nonZero:
        maxFOM_nonZero = myFOM_nonZero
        bdtConf_best_nonZero = bdtConf
        isoConf_best_nonZero = isoConf
        isoVersion_best_nonZero = isoVersion
        newFlag_best_nonZero = newFlag
        bdtCut_best_nonZero = bdtCut_nonZero
    if myFOM_Zero > maxFOM_Zero:
        maxFOM_Zero = myFOM_Zero
        bdtConf_best_Zero = bdtConf
        newFlag_best_Zero = newFlag
        bdtCut_best_Zero = bdtCut_Zero
print ('Best Config for nonZeroTracks is isoVersion', isoVersion_best_nonZero,
       'isoConf', isoConf_best_nonZero, 'FinalBDTconf', bdtConf_best_nonZero,
       'newFlag', newFlag_best_nonZero, 'with Final BDT > ',
       bdtCut_best_nonZero)
print ('Best Config for ZeroTracks is FinalBDTconf', bdtConf_best_Zero,
       'newFlag', newFlag_best_Zero, 'with Final BDT > ', bdtCut_best_Zero)
