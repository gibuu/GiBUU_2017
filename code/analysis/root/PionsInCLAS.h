
#ifndef PionsInCLAS_h
#define PionsInCLAS_h

#include <TROOT.h>
#include <TFile.h>


#include <TVector3.h>
#include <TLorentzVector.h>


#include <vector>
#include <cmath>

#include "ConstantsGiBUU.h"
#include "MuonKinematics.h"


class PionsInCLAS : public MuonKinematics{

   private:
//   TLorentzVector   muon4Momentum, boson4Momentum, hadron4Momentum;
//   Float_t          Enu;
     TLorentzVector   hadron4Momentum;


   public:
   // contructor
   PionsInCLAS(Float_t enu=1.0, Float_t emu=1.0, Float_t pxmu=0, Float_t pymu=0, Float_t pzmu=0, Float_t eH=1.0, Float_t pxH=0, Float_t pyH=0, Float_t pzH=0);
//   ~PionsInCLAS();


//   Float_t Q2();
//   Float_t Wfree();
//   Float_t thetamu();
   Float_t ppi();
   Float_t thetapi();


   Bool_t SelectOutPionAngle(const Float_t angleLowerThreshold=(7./180.*TMath::Pi()), const Float_t angleUpperThreshold=(55./180.*TMath::Pi()) );

   Bool_t SelectOutPionMomentum(const Float_t momentumLowerThreshold=0.1, const Float_t momentumUpperThreshold=4.0);



};


//constructor
PionsInCLAS::PionsInCLAS(Float_t enu, Float_t emu, Float_t pxmu, Float_t pymu, Float_t pzmu, Float_t eH, Float_t pxH, Float_t pyH, Float_t pzH) : MuonKinematics(enu,emu,pxmu,pymu,pzmu)
{
hadron4Momentum.SetPxPyPzE(pxH,pyH,pzH,eH);
}


//destructor
//PionsInCLAS::~PionsInCLAS()  : ~MuonKinematics
//{
//hadron4Momentum.SetPxPyPzE(0.0,0.0,0.0,1.0);
//}




//Float_t PionsInCLAS::Q2(){ return -boson4Momentum.Mag2(); }
//Float_t PionsInCLAS::Wfree(){Float_t W2=m_N2 + 2*m_N*boson4Momentum.E() + boson4Momentum.Mag2();   return sqrt(W2); }
//Float_t PionsInCLAS::thetamu(){ return muon4Momentum.Theta();  } 

Float_t PionsInCLAS::ppi(){ return hadron4Momentum.Rho(); }

Float_t PionsInCLAS::thetapi(){ return hadron4Momentum.Theta();  } 




Bool_t PionsInCLAS::SelectOutPionAngle(const Float_t angleLowerThreshold, const Float_t angleUpperThreshold)
{ 
return ( thetapi()>angleLowerThreshold && thetapi()<angleUpperThreshold );
}


Bool_t PionsInCLAS::SelectOutPionMomentum(const Float_t momentumLowerThreshold, const Float_t momentumUpperThreshold)
{ 
return ( ppi()>momentumLowerThreshold && ppi()<momentumUpperThreshold );
}




#endif

