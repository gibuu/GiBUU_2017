//
//
//
//
//  


#define Analyze_cxx
#include "Analyze.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <THStack.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>

#include "ConstantsGiBUU.h"

#include <iostream>
#include <fstream>
#include <vector>


using namespace std;


#include "Loop.C"


#include "LoopMuonKinematics.C"


#include "LoopEnergyReconstruction.C"


#include "LoopPionsInCLAS.C"


// #include "LoopInclusiveElectronDdiff.C"


// #include "LoopNumberOutgoing.C"



