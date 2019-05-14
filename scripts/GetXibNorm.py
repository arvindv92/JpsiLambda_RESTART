from __future__ import division
import math
from ROOT import TFile, TGraph, gDirectory


def GetNorm(run=1, isoVersion="v0", isoConf=1, finalBDTConf_nonZero=1,
            finalBDTConf_Zero=1, bdtCut_nonZero=-1.0, bdtCut_Zero=-1.0):
    # Get generated yield of Xib -> J/psi Xi MC ###############
    genLog = open("../logs/mc/JpsiXi/run" + str(run) + "/gen_log.txt")
    line = genLog.readline()
    genYield = int(line.rstrip('\n'))
    ###############################

    # Get fitted yield from fully reco'd Xib -> J/psi Xi in data
    xibDataLog = open("../logs/data/JpsiXi/run" + str(run)
                      + "/Fit_JpsiXi_log.txt")
    lines = xibDataLog.readlines()
    for line in lines:
        if 'RooRealVar::nsig' in line:
            xibYield = float(line.split()[2])
            xibYieldErr = float(line.split()[4])
    relErr_xibYield = xibYieldErr / xibYield

    print 'xibYield = ' + str(xibYield) + '+/-' + str(xibYieldErr)
    ###############################

    # Get reconstruction eff of Xib -> J/psi Xi from MC
    xibMcLog = open("../logs/mc/JpsiXi/run" + str(run)
                    + "/Cuts_JpsiXi_log.txt")
    lines = xibMcLog.readlines()
    for line in lines:
        if 'inclusive' in line:
            xibEff = float(line.split()[7]) / 100
            xibEffErr = float(line.split()[10]) / 100
    relErr_xibEff = xibEffErr / xibEff

    file_Xi_rec = TFile("../rootFiles/mcFiles/JpsiXi/run{}/"
                        "jpsixi_cut_LL.root".format(run))
    tree_Xi_rec = file_Xi_rec.MyTuple

    file_Xi_gen = TFile("../rootFiles/mcFiles/JpsiXi/run{}/RW/"
                        "gbWeights_gen.root".format(run))
    tree_Xi_gen = file_Xi_gen.MyTuple

    tree_Xi_rec.Draw("gb_wts>>hxi_rec", "", "goff")
    tree_Xi_gen.Draw("gb_wts>>hxi_gen", "", "goff")

    hxi_rec = gDirectory.Get("hxi_rec")
    hxi_gen = gDirectory.Get("hxi_gen")

    xibnum_wt = hxi_rec.GetEntries() * hxi_rec.GetMean()
    xibden_wt = hxi_gen.GetEntries() * hxi_gen.GetMean()
    xibEff_wt = xibnum_wt / xibden_wt
    xibEffErr_wt = math.sqrt(xibEff_wt * (1 - xibEff_wt) / xibden_wt)
    relErr_xibEff_wt = xibEffErr_wt / xibEff_wt

    print 'xibEff = ' + str('%.4f' % (xibEff_wt * 100)) + ' % +/- ' + str('%.4f' % (xibEffErr_wt * 100)) + ' %'
    ###############################

    # Get efficiency for reco'ing Xib -> J/psi Lambda
    # NB Truth Matching used to get eff.
    path = "../rootFiles/mcFiles/JpsiLambda/JpsiXi/run" + str(run) + "/"
    nonZeroTracksFile = TFile(path
                              + "jpsixi_cutoutks_LL_nonZeroTracks_noPID.root")
    ZeroTracksFile = TFile(path
                           + "jpsixi_cutoutks_LL_ZeroTracks_noPID.root")
    nonZeroTracksTree = nonZeroTracksFile.MyTuple
    ZeroTracksTree = ZeroTracksFile.MyTuple

    nonZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_LL_FinalBDT"
                                + str(finalBDTConf_nonZero) + "_iso"
                                + str(isoConf) + "_" + isoVersion
                                + "_noPID.root")
    ZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_zeroTracksLL_FinalBDT"
                             + str(finalBDTConf_Zero) + "_noPID.root")

    num = nonZeroTracksTree.GetEntries("Lb_BKGCAT==40 && BDT"
                                       + str(finalBDTConf_nonZero) + ">"
                                       + str(bdtCut_nonZero))
    + ZeroTracksTree.GetEntries("Lb_BKGCAT==40 && BDT" + str(finalBDTConf_Zero)
                                + ">" + str(bdtCut_Zero))

    nonZeroTracksTree.Draw("gb_wts>>h0", "Lb_BKGCAT==40 && BDT"
                           + str(finalBDTConf_nonZero) + ">"
                           + str(bdtCut_nonZero), "goff")
    ZeroTracksTree.Draw("gb_wts>>h1", "Lb_BKGCAT==40 && BDT"
                        + str(finalBDTConf_Zero) + ">"
                        + str(bdtCut_Zero), "goff")
    h0 = gDirectory.Get("h0")
    h1 = gDirectory.Get("h1")

    num_wt = h0.GetEntries() * h0.GetMean() + h1.GetEntries() * h1.GetMean()

    genWtsFile = TFile(path + "RW/gbWeights_gen.root")
    genWtsTree = genWtsFile.MyTuple
    genWtsTree.Draw("gb_wts>>hgen", "", "goff")

    hgen = gDirectory.Get("hgen")
    den_wt = hgen.GetEntries() * hgen.GetMean()

    xibEff_JpsiLambda_wt = (num_wt / den_wt)
    # print "num =" + str(num)

    xibEff_JpsiLambda = (num / genYield)

    if num == 0:
        xibEffErr_JpsiLambda = 0.0
    else:
        xibEffErr_JpsiLambda = xibEff_JpsiLambda * math.sqrt((1.0 / num)
                                                             + (1.0 / genYield))
    if num_wt == 0:
        xibEffErr_JpsiLambda_wt = 0.0
    else:
        xibEffErr_JpsiLambda_wt = xibEff_JpsiLambda_wt * math.sqrt((1.0 / num_wt)
                                                                   + (1.0 / den_wt))
    # print xibEff_JpsiLambda
    print 'xibEff_JpsiLambda = ' + str('%.4f' % (xibEff_JpsiLambda_wt * 100)) + ' % +/-' + str('%.4f' % (xibEffErr_JpsiLambda_wt * 100)) + ' %'

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
    ###############################
    # Determine normalization
    if xibEff_JpsiLambda > 0:
        relErr_xibEff_JpsiLambda = xibEffErr_JpsiLambda / xibEff_JpsiLambda
        xibNorm = xibYield * xibEff_JpsiLambda / xibEff
        xibNormErr = xibNorm * math.sqrt(pow(relErr_xibYield, 2)
                                         + pow(relErr_xibEff, 2)
                                         + pow(relErr_xibEff_JpsiLambda, 2))
    else:
        xibNorm = 0.0
        xibNormErr = 0.0
    if xibEff_JpsiLambda_wt > 0:
        relErr_xibEff_JpsiLambda_wt = xibEffErr_JpsiLambda_wt / xibEff_JpsiLambda_wt
        xibNorm_wt = xibYield * xibEff_JpsiLambda_wt / xibEff_wt
        xibNormErr_wt = xibNorm_wt * math.sqrt(pow(relErr_xibYield, 2)
                                               + pow(relErr_xibEff_wt, 2)
                                               + pow(relErr_xibEff_JpsiLambda_wt, 2))
    else:
        xibNorm_wt = 0.0
        xibNormErr_wt = 0.0
    print 'Norm = ' + str(xibNorm_wt) + '+/-' + str(xibNormErr_wt)
    xibNormLog = open("../logs/mc/JpsiXi/run" + str(run) + "/xibNorm_log.txt",
                      "w")
    xibNormLog.write(str(xibNorm) + "\n")
    xibNormLog.write(str(xibNormErr) + "\n")
    xibNormLog.write(str(xibNorm_wt) + "\n")
    xibNormLog.write(str(xibNormErr_wt) + "\n")
    return xibNorm, xibNormErr
