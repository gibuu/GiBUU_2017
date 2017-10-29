
void Analyze::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Analyze.C
//      Root > Analyze t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;


// histograms versus neutrino energy
Int_t kk=40;
Double_t Enu_min=0.0, Enu_max=4.0;
TH1 * eventsQE    = new TH1D ("eventsQE",    "event distribution for QE scattering vs neutrino energy",    kk, Enu_min, Enu_max);  // 40 bins for energy from 0 to 4
Double_t EnergyBinWidth = eventsQE->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs enu
TH1 * eventsDelta = new TH1D ("eventsDelta", "event distribution for Delta production vs neutrino energy", kk, Enu_min, Enu_max);  // 40 bins for energy from 0 to 4




//histograms vs kinetic energy
Int_t nT=40;
Double_t T_min=0, T_max=2.0;
TH1 * events1Piplus_vs_kineticE  = new TH1D ("events1Piplus_vs_kineticE", "event distribution for 1pi+ events vs pion kinetic energy", nT, T_min, T_max);
Double_t KinetEnergyBinWidth = events1Piplus_vs_kineticE->GetBinWidth(0);     // bins width to be used in determining xsec as perweight/bin vs kinetic energy
TH1 * events1Proton_vs_kineticE  = new TH1D ("events1Proton_vs_kineticE",  "event distribution for 1proton events vs proton kinetic energy",   nT, T_min, T_max);
TH1 * events1Neutron_vs_kineticE = new TH1D ("events1Neutron_vs_kineticE", "event distribution for 1neutron events vs neutron kinetic energy", nT, T_min, T_max);


















// initializing number of particles for counting
std::vector<Int_t> numberOutPions(4);       //[0] is total number of pions,    [1] of pi-, [2] of pi0, [3] of pi+
std::vector<Int_t> numberOutNucleons(4);    //[0] is total number of nucleons, [1] of antiprotons, [2] of neutrons, [3] of protons




   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

//#######################
//####   LOOP OVER ENTRIES   BEGINS
//#######################

   // loop-for-entries
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



// set number of each particle species to 0 
for(Int_t i=0; i<4; i++){numberOutPions[i]=0; numberOutNucleons[i]=0;}





//########################################################################  distributions BEFORE particle counting
//event distribution for QE processes 
if(origin==1){    eventsQE->Fill( enu, perweight/EnergyBinWidth ); } 
if(origin==2){ eventsDelta->Fill( enu, perweight/EnergyBinWidth ); } 
      







//###########################################################################   counting particles in each event



//cout<<"    origin="<<origin<<"    enu= "<<enu<<"     perweight="<<perweight<<endl;
//cout<<"numberOutParticles= "<<numberOutParticles<<endl;

	    // particle-loop-1  for counting particles in each event 
	    for (Int_t  j=0; j<numberOutParticles; j++){
	               //cout<<"\t\t\t j= "<<j<<"   particleId="<<particleID[j]<<endl;
	        //pions
	        if(particleID[j]==101){ (numberOutPions[0])++;    (numberOutPions[(charge[j]+2)])++;      }
	        //nucleons
	        if(particleID[j]==1){   (numberOutNucleons[0])++; (numberOutNucleons[(charge[j]+2)])++;  }




	    } // end particle-loop-1





//if(numberOutPions[0]>1){ cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
//check that the pion number is consistent
if(numberOutPions[0]!=(numberOutPions[1]+numberOutPions[2]+numberOutPions[3]))
                        { cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
if(numberOutNucleons[0]!=(numberOutNucleons[1]+numberOutNucleons[2]+numberOutNucleons[3]))
                        { cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl; }








//########################################################################  distributions AFTER particle counting
//event distribution for QE processes 


	    // ////////////  particle-loop-2  for distributions of the outgoin particles 
	    for (Int_t  j=0; j<numberOutParticles; j++){

	    //1-pi+ events 
            if(numberOutPions[0]==1 && numberOutPions[3]==1){
                 if(particleID[j]==101 && charge[j]==1){ events1Piplus_vs_kineticE->Fill(energy[j]-0.140, perweight/KinetEnergyBinWidth ); }   //distribution for 1pi+ 
                                                            }
	    //1-proton events
            if(numberOutNucleons[0]==1 && numberOutNucleons[3]==1){
                 if(particleID[j]==1 && charge[j]==1){ events1Proton_vs_kineticE->Fill(energy[j]-0.938, perweight/KinetEnergyBinWidth ); }   //distribution for 1proton
                                                                  }
	    //1-neutron events
            if(numberOutNucleons[0]==1 && numberOutNucleons[2]==1){
                 if(particleID[j]==1 && charge[j]==0){  events1Neutron_vs_kineticE->Fill(energy[j]-0.938, perweight/KinetEnergyBinWidth ); }   //distribution for 1neutron
								  }





	    } // end particle-loop-2





   } // end loop-for-entries

//#######################
//####   LOOP OVER ENTRIES   IS OVER
//#######################






// draw QE event distribution 
eventsQE->Draw();
eventsDelta->Draw();



//canvas with 1-particle distributions
TCanvas * dEkin = new TCanvas("dEkin", "1-particle distributions versus kinetic energy");
dEkin->Divide(2,2);
dEkin->cd(1);
events1Piplus_vs_kineticE->Draw();
dEkin->cd(2);
events1Proton_vs_kineticE->Draw();
dEkin->cd(3);
events1Neutron_vs_kineticE->Draw();
dEkin->cd(4);





}  // end Analyze::Loop()




