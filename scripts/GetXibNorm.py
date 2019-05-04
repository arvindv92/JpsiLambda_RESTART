from __future__ import division
import math
from ROOT import TFile, TGraph


def GetNorm(run=1, isoVersion="v0", isoConf=1, finalBDTConf_nonZero=1,
            finalBDTConf_Zero=1, bdtCut_nonZero=-1.0, bdtCut_Zero=-1.0):
    genLog = open("../logs/mc/JpsiXi/run" + str(run) + "/gen_log.txt")
    line = genLog.readline()
    genYield = int(line.rstrip('\n'))

    xibDataLog = open("../logs/data/JpsiXi/run" + str(run)
                      + "/Fit_JpsiXi_log.txt")
    lines = xibDataLog.readlines()
    for line in lines:
        if 'RooRealVar::nsig' in line:
            xibYield = float(line.split()[2])
            xibYieldErr = float(line.split()[4])
    relErr_xibYield = xibYieldErr / xibYield

    # print 'xibYield = ' + str(xibYield) + '+/-' + str(xibYieldErr)
    ###############################
    xibMcLog = open("../logs/mc/JpsiXi/run" + str(run)
                    + "/cuts_JpsiXi_log.txt")
    lines = xibMcLog.readlines()
    for line in lines:
        if 'inclusive' in line:
            xibEff = float(line.split()[7]) / 100
            xibEffErr = float(line.split()[10]) / 100
    relErr_xibEff = xibEffErr / xibEff

    # print 'xibEff = ' + str(xibEff) + '+/-' + str(xibEffErr)
    ###############################
    path = "../rootFiles/mcFiles/JpsiLambda/JpsiXi/run" + str(run) + "/"
    nonZeroTracksFile = TFile(path
                              + "jpsixi_cutoutks_LL_nonZeroTracks_noPID.root")
    ZeroTracksFile = TFile(path
                           + "jpsixi_cutoutks_LL_ZeroTracks_noPID.root")
    nonZeroTracksTree = nonZeroTracksFile.MyTuple
    ZeroTracksTree = ZeroTracksFile.MyTuple
    ###############################
    nonZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_LL_FinalBDT"
                                + str(finalBDTConf_nonZero) + "_iso" + str(isoConf)
                                + "_" + isoVersion + "_noPID.root")
    ZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_zeroTracksLL_FinalBDT"
                             + str(finalBDTConf_Zero) + "_noPID.root")
    ###############################

    num = nonZeroTracksTree.GetEntries("Lb_BKGCAT==40 && BDT"
                                       + str(finalBDTConf_nonZero) + ">"
                                       + str(bdtCut_nonZero))
    + ZeroTracksTree.GetEntries("Lb_BKGCAT==40 && BDT" + str(finalBDTConf_Zero)
                                + ">" + str(bdtCut_Zero))
    # print "num =" + str(num)

    xibEff_JpsiLambda = (num / genYield)

    if num == 0:
        xibEffErr_JpsiLambda = 0.0
    else:
        xibEffErr_JpsiLambda = xibEff_JpsiLambda * math.sqrt((1.0 / num)
                                                             + (1.0 / genYield))
    # print xibEff_JpsiLambda

    # xibMcLog_JpsiLambda = open("../logs/mc/JpsiLambda/JpsiXi/run"
    #                            + str(run) + "/CutFinalBDT"
    #                            + str(finalBDTConf)
    #                            + "_LL_iso" + str(isoConf) + "_" + isoVersion
    #                            + ".txt")
    # lines = xibMcLog_JpsiLambda.readlines()
    # for line in lines:
    #     if 'inclusive' in line:
    #         line = line.translate(None, '%')
    #         xibEff_JpsiLambda = float(line.split()[8]) / 100
    #         xibEffErr_JpsiLambda = float(line.split()[10]) / 100

    if xibEff_JpsiLambda > 0:
        relErr_xibEff_JpsiLambda = xibEffErr_JpsiLambda / xibEff_JpsiLambda
        xibNorm = xibYield * xibEff_JpsiLambda / xibEff
        xibNormErr = xibNorm * math.sqrt(pow(relErr_xibYield, 2)
                                         + pow(relErr_xibEff, 2)
                                         + pow(relErr_xibEff_JpsiLambda, 2))
    else:
        xibNorm = 0.0
        xibNormErr = 0.0
    xibNormLog = open("../logs/mc/JpsiXi/run" + str(run) + "/xibNorm_log.txt",
                      "w")
    xibNormLog.write(str(xibNorm) + "\n")
    xibNormLog.write(str(xibNormErr) + "\n")
    return xibNorm, xibNormErr
