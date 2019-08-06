from __future__ import division
import math
from ROOT import TFile, gDirectory


def GetNorm(run=1, isoVersion="v0", isoConf=1, finalBDTConf_nonZero=1,
            finalBDTConf_Zero=1, bdtCut_nonZero=-1.0, bdtCut_Zero=-1.0,
            shift_trEff=0.0):
    # (tau Xib-/tau Xib0)
    tauXibmin_tauXib0 = 1.087 # despite the name of the variable, this is (tau(Xib-)*B(Xi- -> Lambda pi-))/(tau(Xib0)*B(Xi0 -> Lambda pi0)
    err_tauXibmin_tauXib0 = 0.036

    scale_factor = 1.0 + (1.0 / tauXibmin_tauXib0)
    err_scale_factor = scale_factor * (err_tauXibmin_tauXib0 / tauXibmin_tauXib0)

    relErr_scale_factor = (err_scale_factor / scale_factor)

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

    if run == 1:
        xibYield_relSyst = 0.09  # relative systematic on Xib yield. Comes from choice of fit model
    elif run == 2:
        xibYield_relSyst = 0.06

    print 'xibYield = ' + str(xibYield) + '+/-' + str(xibYieldErr) + '+/-' + str(xibYield * xibYield_relSyst)
    ###############################

    # Get reconstruction eff of Xib -> J/psi Xi from MC
    xibMcLog = open("../logs/mc/JpsiXi/run" + str(run)
                    + "/Cuts_JpsiXi_log.txt")
    lines = xibMcLog.readlines()
    for line in lines:
        if 'inclusive' in line:
            xibEff = float(line.split()[7]) / 100  # Unweighted efficiency for Xib -> J/psi Xi
            xibEffErr = float(line.split()[10]) / 100
    relErr_xibEff = xibEffErr / xibEff

    # Now get the weighted efficiency
    file_Xi_rec = TFile("../rootFiles/mcFiles/JpsiXi/run{}/"
                        "jpsixi_cut_LL.root".format(run))
    tree_Xi_rec = file_Xi_rec.MyTuple

    file_Xi_gen = TFile("../rootFiles/mcFiles/JpsiXi/run{}/RW/"
                        "gbWeights_gen.root".format(run))
    tree_Xi_gen = file_Xi_gen.MyTuple

    tree_Xi_rec.Draw("GB_WT*(wt_tracking+{}*wtErr_tracking)>>hxi_rec".format(shift_trEff), "", "goff")
    tree_Xi_gen.Draw("GB_WT>>hxi_gen", "", "goff")

    hxi_rec = gDirectory.Get("hxi_rec")
    hxi_gen = gDirectory.Get("hxi_gen")

    xibnum_wt = hxi_rec.GetEntries() * hxi_rec.GetMean()
    xibden_wt = hxi_gen.GetEntries() * hxi_gen.GetMean()

    xibEff_wt = xibnum_wt / xibden_wt  # Weighted efficiency for Xib->JpsiXi
    xibEffErr_wt = math.sqrt(xibEff_wt * (1 - xibEff_wt) / xibden_wt)  # Binomial uncertainty

    relErr_xibEff_wt = xibEffErr_wt / xibEff_wt

    print 'UNWEIGHTED Xib->JpsiXi Eff = ' + str('%.4f' % (xibEff * 100)) + ' % +/- ' + str('%.4f' % (xibEffErr * 100)) + ' %'
    print 'WEIGHTED Xib->JpsiXi Eff = ' + str('%.4f' % (xibEff_wt * 100)) + ' % +/- ' + str('%.4f' % (xibEffErr_wt * 100)) + ' %'
    ###############################

    # Get efficiency for reco'ing Xib -> J/psi Lambda
    # NB Not using TM anymore.
    path = "../rootFiles/mcFiles/JpsiLambda/JpsiXi/run" + str(run) + "/"
    nonZeroTracksFile = TFile(path
                              + "jpsixi_cutoutks_LL_nonZeroTracks.root")
    ZeroTracksFile = TFile(path
                           + "jpsixi_cutoutks_LL_ZeroTracks.root")
    nonZeroTracksTree = nonZeroTracksFile.MyTuple
    ZeroTracksTree = ZeroTracksFile.MyTuple

    nonZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_LL_FinalBDT"
                                + str(finalBDTConf_nonZero) + "_iso"
                                + str(isoConf) + "_" + isoVersion
                                + ".root")
    ZeroTracksTree.AddFriend("MyTuple", path + "jpsixi_zeroTracksLL_FinalBDT"
                             + str(finalBDTConf_Zero) + ".root")

    num = nonZeroTracksTree.GetEntries("BDT" + str(finalBDTConf_nonZero) + ">"
                                       + str(bdtCut_nonZero))
    + ZeroTracksTree.GetEntries("BDT" + str(finalBDTConf_Zero)
                                + ">" + str(bdtCut_Zero))  # counting no. of entries passing BDT cut.

    nonZeroTracksTree.Draw("GB_WT>>h0", "BDT"
                           + str(finalBDTConf_nonZero) + ">"
                           + str(bdtCut_nonZero), "goff")
    ZeroTracksTree.Draw("GB_WT>>h1", "BDT"
                        + str(finalBDTConf_Zero) + ">"
                        + str(bdtCut_Zero), "goff")

    h0 = gDirectory.Get("h0")
    h1 = gDirectory.Get("h1")

    num_wt = h0.GetEntries() * h0.GetMean() + h1.GetEntries() * h1.GetMean()  # sum of weights passing BDT cut

    genWtsFile = TFile(path + "RW/gbWeights_gen.root")  # Weighted generator MC
    genWtsTree = genWtsFile.MyTuple

    genWtsTree.Draw("GB_WT>>hgen", "", "goff")
    hgen = gDirectory.Get("hgen")
    den_wt = hgen.GetEntries() * hgen.GetMean()  # sum of wts. for generated events

    xibEff_JpsiLambda_wt = (num_wt / den_wt)  # Weighted efficiency for Xib->JpsiLambda

    xibEff_JpsiLambda = (num / genYield)  # Unweighted efficiency for Xib->JpsiLambda

    xibEffErr_JpsiLambda = math.sqrt((xibEff_JpsiLambda * (1 - xibEff_JpsiLambda)) / genYield)
    xibEffErr_JpsiLambda_wt = math.sqrt((xibEff_JpsiLambda_wt * (1 - xibEff_JpsiLambda_wt)) / den_wt)

    relErr_xibEff_JpsiLambda = xibEffErr_JpsiLambda / xibEff_JpsiLambda
    relErr_xibEff_JpsiLambda_wt = xibEffErr_JpsiLambda_wt / xibEff_JpsiLambda_wt

    # print xibEff_JpsiLambda
    print 'UNWEIGHTED Xib->JpsiLambda Eff = ' + str('%.4f' % (xibEff_JpsiLambda * 100)) + ' % +/-' + str('%.4f' % (xibEffErr_JpsiLambda * 100)) + ' %'
    print 'WEIGHTED Xib->JpsiLambda Eff = ' + str('%.4f' % (xibEff_JpsiLambda_wt * 100)) + ' % +/-' + str('%.4f' % (xibEffErr_JpsiLambda_wt * 100)) + ' %'

    ###############################
    # Determine normalization

    trackingErr = 0.05  # rel error. includes 4.5% tracking unc, and 2% uncertainty for material effects from hadronic interactions
    xiVtxUnc = 0.014  # rel error. Unc from vertexing the Xi-. Comes from Steves ANA

    # Unweighted Norm
    if xibEff_JpsiLambda > 0:
        xibNorm = scale_factor * xibYield * xibEff_JpsiLambda / xibEff  # NB: Accounting for Xib0 right here
        xibNormErr_stat = scale_factor * xibYieldErr * xibEff_JpsiLambda / xibEff
        xibNormErr_syst = xibNorm * math.sqrt(pow(trackingErr, 2) + pow(xibYield_relSyst, 2)
                                              + pow(relErr_xibEff_JpsiLambda, 2) + pow(relErr_xibEff, 2)
                                              + pow(xiVtxUnc, 2) + pow(relErr_scale_factor, 2))
        xibNormErr = math.sqrt(pow(xibNormErr_stat, 2) + pow(xibNormErr_syst, 2))
    else:
        xibNorm = 0.0
        xibNormErr_stat = 0.0
        xibNormErr_syst = 0.0
        xibNormErr = 0.0

    # Weighted Norm
    if xibEff_JpsiLambda_wt > 0:
        xibNorm_wt = scale_factor * xibYield * xibEff_JpsiLambda_wt / xibEff_wt
        xibNormErr_wt_stat = scale_factor * xibYieldErr * xibEff_JpsiLambda_wt / xibEff_wt
        xibNormErr_wt_syst = xibNorm_wt * math.sqrt(pow(trackingErr, 2) + pow(xibYield_relSyst, 2)
                                                    + pow(relErr_xibEff_JpsiLambda_wt, 2) + pow(relErr_xibEff_wt, 2)
                                                    + pow(xiVtxUnc, 2) + pow(relErr_scale_factor, 2))
        xibNormErr_wt = math.sqrt(pow(xibNormErr_wt_stat, 2) + pow(xibNormErr_wt_syst, 2))
    else:
        xibNorm_wt = 0.0
        xibNormErr_wt_stat = 0.0
        xibNormErr_wt_syst = 0.0
        xibNormErr_wt = 0.0

    print 'UNWEIGTED Norm = ' + str(xibNorm) + '+/-' + str(xibNormErr_stat) + '+/-' + str(xibNormErr_syst)
    print 'WEIGHTED Norm = ' + str(xibNorm_wt) + '+/-' + str(xibNormErr_wt_stat) + '+/-' + str(xibNormErr_wt_syst)
    print 'BREAKUP OF SYSTEMATICS:'
    print 'Overall Syst     = ' + str((xibNormErr_wt_syst / xibNorm_wt) * 100) + '%'
    print '*Tracking        = ' + str(trackingErr * 100) + '%'
    print '*Xi Vtx          = ' + str(xiVtxUnc * 100) + '%'
    print '*Xib Yield       = ' + str(xibYield_relSyst * 100) + '%'
    print '*Xib->JpsiL eff  = ' + str(relErr_xibEff_JpsiLambda_wt * 100) + '%'
    print '*Xib->JpsiXi eff = ' + str(relErr_xibEff_wt * 100) + '%'
    print '*tauXib-/tauXib0 = ' + str(relErr_scale_factor * 100) + '%'
    # Write results out to log file
    xibNormLog = open("../logs/mc/JpsiXi/run" + str(run) + "/xibNorm_log.txt",
                      "w")
    xibNormLog.write(str(xibNorm) + "\n")
    xibNormLog.write(str(xibNormErr_stat) + "\n")
    xibNormLog.write(str(xibNormErr_syst) + "\n")
    xibNormLog.write(str(xibNorm_wt) + "\n")
    xibNormLog.write(str(xibNormErr_wt_stat) + "\n")
    xibNormLog.write(str(xibNormErr_wt_syst) + "\n")
    return xibNorm, xibNormErr
