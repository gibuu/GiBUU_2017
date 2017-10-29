
#ifndef EnergyReconstruction_h
#define EnergyReconstruction_h

#include <TROOT.h>
#include <TFile.h>


#include <TVector3.h>
#include <TLorentzVector.h>


#include <vector>

class EnergyReconstruction{

   private:
//   TVector3       neutrinoDirection(0.0,0.0,1.0);
//   TVector3       muon3Momentum, hadron3Momentum;
//   Float_t        muonAbsMomentum, hadronAbsMomentum;
   TLorentzVector   muon4Momentum, hadron4Momentum;
//   Float_t        cosNuMuonAngle, cosNuHadronAngle;



//   static const  Float_t  m_mu2=0.105*0.105, m_pi2=0.140*0.140, m_N=0.938, m_N2=0.938*0.938, M_Delta2=1.232*1.232; 

   public:
   // contructor
   EnergyReconstruction(Float_t emu=1.0, Float_t pxmu=0, Float_t pymu=0, Float_t pzmu=0, Float_t eH=1.0, Float_t pxH=0, Float_t pyH=0, Float_t pzH=0);
   ~EnergyReconstruction();

   Float_t OnePionEnergyReconstruction();
   Float_t DeltaEnergyReconstruction();


};


//constructor
EnergyReconstruction::EnergyReconstruction(Float_t emu, Float_t pxmu, Float_t pymu, Float_t pzmu, Float_t eH, Float_t pxH, Float_t pyH, Float_t pzH)
{
//muon3Momentum.SetXYZ(pxmu,pymu,pzmu);
//hadron3Momentum.SetXYZ(pxH,pyH,pzH);
muon4Momentum.SetPxPyPzE(pxmu,pymu,pzmu,emu);
hadron4Momentum.SetPxPyPzE(pxH,pyH,pzH,eH);
//muonAbsMomentum
//hadronAbsMomentum
}



EnergyReconstruction::~EnergyReconstruction()
{
muon4Momentum.SetPxPyPzE(0.0,0.0,0.0,1.0);
hadron4Momentum.SetPxPyPzE(0.0,0.0,0.0,1.0);
}





// assumes hadron is pion
// used for 1pi+ events by MiniBooNE
// formular is taken from "Measurement of the neutrino-induced charged-current charged pion production ..." PRD83(2011)052007
// the same formula but with other notations is used for 1pi0 events, see Nelson "CCpi0 event reconstruction at MinibooNE" Proc. NuInt09
Float_t EnergyReconstruction::OnePionEnergyReconstruction()
{
Float_t Erec;

//muon4Momentum.Print();
//hadron4Momentum.Print();

Erec=( m_mu2 + m_pi2 -2.*m_N*(muon4Momentum.E()+hadron4Momentum.E()) + 2.*muon4Momentum*hadron4Momentum )/2./( muon4Momentum.E() + hadron4Momentum.E() - muon4Momentum.Pz() - hadron4Momentum.Pz() - m_N);

return Erec;
}








// assumes that Delta is produced on-shell
// earlier used for 1pi+ events by MiniBooNE 
// formular is taken from taken from PRL103 (2009) 081801  "Measurement of the ratio ..."
Float_t EnergyReconstruction::DeltaEnergyReconstruction()
{
Float_t Erec;

//muon4Momentum.Print();
//hadron4Momentum.Print();

Erec=( 2.*m_N*muon4Momentum.E() + M_Delta2 - m_N2 -m_mu2 )/2./( m_N - muon4Momentum.E()  + muon4Momentum.Pz() );

return Erec;
}









#endif

