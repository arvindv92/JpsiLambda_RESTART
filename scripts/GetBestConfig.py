import glob
import argparse

parser = argparse.ArgumentParser(description='Get Best Config')

parser.add_argument('-run', '--RUN', type=int, nargs=1,
                    help='Which run do you want best config for?')
parser.add_argument('-iso', '--ISO', type=bool, nargs=1,
                    help='IsoFlag?')

args = parser.parse_args()
run = args.RUN
isoFlag = args.ISO

# isoFlag = True
files = glob.glob("../logs/data/JpsiLambda/run{}/OptimizeFinalBDT*_Punzi_Sigma.txt".format(run))
ctr = 0
myFOM_nonZero = 0.0
myFOM_Zero = 0.0
mySigEff_nonZero = 0.0
mySigEff_Zero = 0.0
myBkgEff_nonZero = 0.0
myBkgEff_Zero = 0.0
maxFOM_nonZero = -1.0
maxFOM_Zero = -1.0
bdtCut_Zero = -1.0
maxSigEff_nonZero = 0.0
maxSigEff_Zero = 0.0
maxBkgEff_nonZero = 0.0
maxBkgEff_Zero = 0.0
bdtCut_nonZero = -1.0
bdtConf_best_nonZero = -1
isoConf_best_nonZero = -1
isoVersion_best_nonZero = -1
# newFlag_best_nonZero = False
bdtCut_best_nonZero = -1.0
bdtConf_best_Zero = -1
isoConf_best_Zero = -1
isoVersion_best_Zero = -1
# newFlag_best_Zero = False
bdtCut_best_Zero = -1.0

for fName in files:
    print fName
    bdtConf = int(fName.split('BDT')[1][0])
    if isoFlag:
        isoConf = int(fName.split('iso')[1][0])
    if isoFlag:
        isoVersion = int(fName.split('_v')[1][0])
    with open(fName) as f:
        ctr = 0
        lines = f.readlines()
        for line in lines:
            line = line.translate(None, '%')
            if 'MAXIMUM FOM' in line:
                if ctr == 0:
                    myFOM_nonZero = float(line.split()[3])
                    bdtCut_nonZero = float(line.split()[7])
                    mySigEff_nonZero = float(line.split()[14])  # this gives weighted signal efficiency
                    myBkgEff_nonZero = float(line.split()[18])
                if ctr == 1:
                    myFOM_Zero = float(line.split()[3])
                    bdtCut_Zero = float(line.split()[7])
                    mySigEff_Zero = float(line.split()[14])  # this gives weighted signal efficiency
                    myBkgEff_Zero = float(line.split()[18])
                ctr = ctr + 1
    print ('finalBDTconf', bdtConf, 'isoConf', isoConf,
           'isoVersion', isoVersion, 'FOM_nonZero', myFOM_nonZero,
           'eff_nonZero', mySigEff_nonZero, 'eff_Zero', mySigEff_Zero)
    if myFOM_nonZero > maxFOM_nonZero:
        maxFOM_nonZero = myFOM_nonZero
        bdtConf_best_nonZero = bdtConf
        isoConf_best_nonZero = isoConf
        isoVersion_best_nonZero = isoVersion
        # newFlag_best_nonZero = newFlag
        bdtCut_best_nonZero = bdtCut_nonZero
        maxSigEff_nonZero = mySigEff_nonZero
        maxBkgEff_nonZero = myBkgEff_nonZero
    if myFOM_Zero > maxFOM_Zero:
        maxFOM_Zero = myFOM_Zero
        bdtConf_best_Zero = bdtConf
        # newFlag_best_Zero = newFlag
        bdtCut_best_Zero = bdtCut_Zero
        maxSigEff_Zero = mySigEff_Zero
        maxBkgEff_Zero = myBkgEff_Zero
    # if myFOM_nonZero > maxFOM_nonZero:
    #     maxFOM_nonZero = myFOM_nonZero
    #     bdtConf_best_nonZero = bdtConf
    #     if isoFlag:
    #         isoConf_best_nonZero = isoConf
    #         isoVersion_best_nonZero = isoVersion
    #     # newFlag_best_nonZero = newFlag
    #     bdtCut_best_nonZero = bdtCut_nonZero
    #     maxSigEff_nonZero = mySigEff_nonZero
    #     maxBkgEff_nonZero = myBkgEff_nonZero
    # if myFOM_Zero > maxFOM_Zero:
    #     maxFOM_Zero = myFOM_Zero
    #     bdtConf_best_Zero = bdtConf
    #     # newFlag_best_Zero = newFlag
    #     bdtCut_best_Zero = bdtCut_Zero
    #     maxSigEff_Zero = mySigEff_Zero
    #     maxBkgEff_Zero = myBkgEff_Zero
print ('Best Config for nonZeroTracks is isoVersion', isoVersion_best_nonZero,
       'isoConf', isoConf_best_nonZero, 'FinalBDTconf', bdtConf_best_nonZero,
       'with Final BDT > ', bdtCut_best_nonZero, 'FOM = ', maxFOM_nonZero,
       'Sig Eff = ', maxSigEff_nonZero, 'Bkg Eff = ', maxBkgEff_nonZero)
print ('Best Config for ZeroTracks is FinalBDTconf', bdtConf_best_Zero,
       'with Final BDT > ', bdtCut_best_Zero, 'FOM = ', maxFOM_Zero,
       'Sig Eff = ', maxSigEff_Zero, 'Bkg Eff = ', maxBkgEff_Zero)
# print ('Best Config for noIso is FinalBDTconf', bdtConf_best_Zero,
#        'with Final BDT > ', bdtCut_best_Zero, 'FOM = ', maxFOM_Zero,
#        'Sig Eff = ', maxSigEff_Zero, 'Bkg Eff = ', maxBkgEff_Zero)
