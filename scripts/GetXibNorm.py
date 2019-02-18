import math


def GetNorm(run=1, isoVersion="v0", isoConf=1, finalBDTConf=1, newFlag=True):
    if newFlag:
        myNew = "_new"
    else:
        myNew = ""
    xibDataLog = open("../logs/data/JpsiXi/run" + str(run)
                      + "/Fit_JpsiXi_log.txt")
    lines = xibDataLog.readlines()
    for line in lines:
        if 'RooRealVar::nsig' in line:
            xibYield = float(line.split()[2])
            xibYieldErr = float(line.split()[4])
    relErr_xibYield = xibYieldErr / xibYield
    xibMcLog = open("../logs/mc/JpsiXi/run" + str(run)
                    + "/cuts_JpsiXi_log.txt")
    lines = xibMcLog.readlines()
    for line in lines:
        if 'inclusive' in line:
            xibEff = float(line.split()[7]) / 100
            xibEffErr = float(line.split()[10]) / 100
    relErr_xibEff = xibEffErr / xibEff
    xibMcLog_JpsiLambda = open("../logs/mc/JpsiLambda/JpsiXi/run"
                               + str(run) + "/CutFinalBDT" + str(finalBDTConf)
                               + "_LL_iso" + str(isoConf) + "_" + isoVersion
                               + myNew + ".txt")
    lines = xibMcLog_JpsiLambda.readlines()
    for line in lines:
        if 'inclusive' in line:
            line = line.translate(None, '%')
            xibEff_JpsiLambda = float(line.split()[8]) / 100
            xibEffErr_JpsiLambda = float(line.split()[10]) / 100
    relErr_xibEff_JpsiLambda = xibEffErr_JpsiLambda / xibEff_JpsiLambda
    xibNorm = xibYield * xibEff_JpsiLambda / xibEff
    xibNormErr = xibNorm * math.sqrt(pow(relErr_xibYield, 2)
                                     + pow(relErr_xibEff, 2)
                                     + pow(relErr_xibEff_JpsiLambda, 2))
    return xibNorm, xibNormErr
