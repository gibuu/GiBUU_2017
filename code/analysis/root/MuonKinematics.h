
#ifndef MuonKinematics_h
#define MuonKinematics_h

#include <TROOT.h>
#include <TFile.h>

#include <TLorentzVector.h>

#include <cmath>

#include "ConstantsGiBUU.h"

class MuonKinematics{

   private:
   TLorentzVector   muon4Momentum, boson4Momentum;
   Float_t          Enu;


   public:
   // contructor
   MuonKinematics(Float_t enu=1.0, Float_t emu=1.0, Float_t pxmu=0, Float_t pymu=0, Float_t pzmu=0);
//   ~MuonKinematics();


   Float_t Nu();
   Float_t Q2();
   Float_t Wfree();
   Float_t thetamu();
   Float_t pmu();
   
   
   Bool_t SelectInElectronEnergy(const Float_t energyLowerThreshold=1.1,                    const Float_t energyUpperThreshold=5.0);
   Bool_t SelectQ2(const Float_t Q2LowerThreshold=1.1,                                    const Float_t Q2UpperThreshold=5.0);
   Bool_t SelectOutElectronY(const Float_t yLowerThreshold=0.0,                           const Float_t yUpperThreshold=0.872);
   Bool_t SelectOutElectronAngle(const Float_t angleLowerThreshold=(7./180.*TMath::Pi()), const Float_t angleUpperThreshold=(55./180.*TMath::Pi()) );

   

};


//constructor
MuonKinematics::MuonKinematics(Float_t enu, Float_t emu, Float_t pxmu, Float_t pymu, Float_t pzmu)
{
Enu=enu;
muon4Momentum.SetPxPyPzE(pxmu,pymu,pzmu,emu);
boson4Momentum.SetPxPyPzE(-pxmu,-pymu,enu-pzmu,enu-emu);
}



//MuonKinematics::~MuonKinematics()
//{
//Enu=1.0;
//muon4Momentum.SetPxPyPzE(0.0,0.0,0.0,1.0);
//boson4Momentum.SetPxPyPzE(0.0,0.0,0.0,1.0);
//}


Float_t MuonKinematics::Nu(){ return boson4Momentum.E(); }

Float_t MuonKinematics::Q2(){ return -boson4Momentum.Mag2(); }

Float_t MuonKinematics::Wfree(){ Float_t W2=m_N2 + 2*m_N*boson4Momentum.E() + boson4Momentum.Mag2();   return sqrt(W2); }

Float_t MuonKinematics::thetamu(){ return muon4Momentum.Theta();  } 

Float_t MuonKinematics::pmu(){ return muon4Momentum.Rho();  } 




Bool_t MuonKinematics::SelectQ2(const Float_t Q2LowerThreshold, const Float_t Q2UpperThreshold)
{
return ( Q2() > Q2LowerThreshold  &&  Q2() < Q2UpperThreshold );
}



Bool_t MuonKinematics::SelectOutElectronY(const Float_t yLowerThreshold, const Float_t yUpperThreshold)
{
return ( boson4Momentum.E()/Enu > yLowerThreshold  &&  boson4Momentum.E()/Enu < yUpperThreshold );
}


Bool_t MuonKinematics::SelectOutElectronAngle(const Float_t angleLowerThreshold, const Float_t angleUpperThreshold)
{ 
return ( thetamu()>angleLowerThreshold && thetamu()<angleUpperThreshold );
}


Bool_t MuonKinematics::SelectInElectronEnergy(const Float_t energyLowerThreshold,  const Float_t energyUpperThreshold)
{ 
return ( Enu>energyLowerThreshold && Enu<energyUpperThreshold );
}




#endif

