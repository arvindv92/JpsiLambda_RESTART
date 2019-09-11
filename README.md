CollateFiles	 - When processing data, this script will form a TChain out of all the raw ROOT files that come from DaVinci. When processing MC, it will hadd the PID corrected ROOT files to form a total ROOT file.
Trigger      	  - This script applies the trigger cut and writes out an output ROOT file. For Run 1, all years are processed simultaneously. For Run 2, each year is processed separately.
Sanity       	   - This script applies the "sanity" cuts, which are loose offline preselections, to the output of Trigger and writes out an output ROOT file. For Run 1, all years are processed simultaneously. For Run 2, each year is processed separately.
DoSWeight_Sanity - This script applies sPlot to the data after Sanity (on events b/w 5.2 GeV and 6 GeV). It writes an output ROOT file for the events in the range with all the original branches + an SW branch.
CutOutKs    	  - This script applies the PID veto to the output of Sanity, to get rid of Ks -> Lambda mis-ID. It writes out separate output ROOT files for the nonZeroTracks and ZeroTracks subsets. For Run 1, all years are processed simultaneously. For Run 2, each year is processed separately. 
DoSWeight   	   - This	script applies sPlot to	the data after CutOutKs (on events b/w 5.2 GeV and 6 GeV). It writes an output ROOT file for the events in the range with all the original branches + an SW branch. nonZeroTracks and ZeroTracks are sWeighted separately
TrainIsolation   - This trains the isolation BDT. Signal training sample: sWeighted nonZeroTracks data. Background training sample: sWeighted B+ -> J/psi K+ data. 
ApplyIsolation   - This applies the isolation BDT to all nonZeroTracks data and MC that has survived CutOutKs. Writes out a ROOT file with a single branch containing the isolation BDT output. This is meant to be attached as a "friend" tree to the data tree.
TrainFinalBDT    - This trains the 2 versions of the final BDT, one with isolation among the training variables (IsoBDT), one without (noIsoBDT). IsoBDT is meant to be applied to nonZeroTracks data. noIsoBDT is meant to be applied to the ZeroTracks data. 
		      For IsoBDT  : Signal Training sample: sWeighted nonZeroTracks data with isolation BDT applied to it. Background training sample: Upper sideband nonZeroTracks events with isolation BDT applied to it. 
		      	     For NoIsoBDT: Signal Training sample: All sWeighted data surviving CutOutKs. Background training sample: All upper sideband events surviving CutOutKs.
ApplyFinalBDT    - This applies the final BDT to all data. For nonZeroTracks data, IsoBDT is applied. For ZeroTracks data, NoIsoBDT is applied. Writes out a ROOT file with a single branch containing the final BDT output. This is meant to be attached as a "friend" tree to the data tree.
OptimizeFinalBDT - This script chooses the final BDT cut(s) to maximize the Punzi FoM.
Fitscript_simul  - Main fit script. Performs a simultaneous binned fit to Run 1 and Run 2 data post all selections. Output of the fit is written in a ModelConfig in a ROOT file.
StdHypoTestInvDemo - This script takes the ModelConfig from the fit as input and performs the UL calculation.
