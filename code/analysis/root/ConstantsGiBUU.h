#ifndef ConstantsGiBUU_h
#define ConstantsGiBUU_h

#include <vector>
#include <TString.h>

using namespace std;

const  Float_t  m_mu2=0.105*0.105, m_pi2=0.140*0.140, m_N=0.938, m_N2=0.938*0.938, M_Delta2=1.232*1.232; 

// 1=QE, 2-Delta, 3=highRES, 4=1piBG 5=DIS 6=2p2hQE, 7=2p2hDelta, 8=2pi 
const Int_t myints[] = {99,1,2,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,3,3,3,3,3,3,3,3,3, 3,4,4,5,6,7,8};
const vector<Int_t> KToHist (myints, myints + sizeof(myints) / sizeof(Int_t) );

const TString HistToName[9] =  {"all","QE","Delta","highRES","1-pi bgr","DIS","2p2h-NN","2p2h-N#Delta","2-pi bgr"};
//vector<TString> HistToName;
//for (Int_t i=0; i<10; i++){HistToName.push_back(mystrings[i]);}

const TString CToCharge[4] =  {"all charges","charge -1","charge 0","charge +1"};



#endif

