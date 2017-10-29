//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 26 17:18:29 2012 by ROOT version 5.28/00c
// from TTree T/events from the ascii file
// found on file: GiBUUEvents.root
//////////////////////////////////////////////////////////

#ifndef Analyze_h
#define Analyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


#include <vector>



class Analyze {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   static const Int_t     kMaxNumberOutParticles=115;

   // Declaration of leaf types
   Int_t           eventNumber;
   Int_t           origin;
   Float_t         perweight;
   Float_t         enu;
   Int_t           numberOutParticles;
   Int_t           particleID[kMaxNumberOutParticles];   //[numberOutParticles]
   Int_t           charge[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         positionX[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         positionY[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         positionZ[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         energy[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         momentumX[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         momentumY[kMaxNumberOutParticles];   //[numberOutParticles]
   Float_t         momentumZ[kMaxNumberOutParticles];   //[numberOutParticles]
   Long64_t        history[kMaxNumberOutParticles];   //[numberOutParticles]

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_origin;   //!
   TBranch        *b_perweight;   //!
   TBranch        *b_enu;   //!
   TBranch        *b_numberOutParticles;   //!
   TBranch        *b_particleID;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_positionX;   //!
   TBranch        *b_positionY;   //!
   TBranch        *b_positionZ;   //!
   TBranch        *b_energy;   //!
   TBranch        *b_momentumX;   //!
   TBranch        *b_momentumY;   //!
   TBranch        *b_momentumZ;   //!
   TBranch        *b_history;   //!

   Analyze(TTree *tree=0);
   virtual ~Analyze();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   virtual void     Loop();
   virtual void     LoopMuonKinematics(TString reaction);
   virtual void     LoopEnergyReconstruction(TString reaction);
   virtual void     LoopPionsInCLAS();
//    virtual void     LoopInclusiveElectronDdiff();
//    virtual void     LoopNumberOutgoing();

};

#endif

#ifdef Analyze_cxx
Analyze::Analyze(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("GiBUUEvents.root");
      if (!f) {
         f = new TFile("GiBUUEvents.root");
      }
      tree = (TTree*)gDirectory->Get("T");

   }
   Init(tree);
}

Analyze::~Analyze()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Analyze::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("origin", &origin, &b_origin);
   fChain->SetBranchAddress("perweight", &perweight, &b_perweight);
   fChain->SetBranchAddress("enu", &enu, &b_enu);
   fChain->SetBranchAddress("numberOutParticles", &numberOutParticles, &b_numberOutParticles);
   fChain->SetBranchAddress("particleID", particleID, &b_particleID);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("positionX", positionX, &b_positionX);
   fChain->SetBranchAddress("positionY", positionY, &b_positionY);
   fChain->SetBranchAddress("positionZ", positionZ, &b_positionZ);
   fChain->SetBranchAddress("energy", energy, &b_energy);
   fChain->SetBranchAddress("momentumX", momentumX, &b_momentumX);
   fChain->SetBranchAddress("momentumY", momentumY, &b_momentumY);
   fChain->SetBranchAddress("momentumZ", momentumZ, &b_momentumZ);
   fChain->SetBranchAddress("history", history, &b_history);
   Notify();
}

Bool_t Analyze::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyze::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyze_cxx
