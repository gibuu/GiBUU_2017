#include "MuonKinematics.h"

void Analyze::LoopMuonKinematics(TString reaction)
{

const Short_t numOrigins=6;

// histograms versus Q2
Int_t kk_Q2=25;
Float_t Q2_min=0.0, Q2_max=2.0;
vector <TH1F *> eventsAll_vs_Q2(numOrigins+1);
for (int i=0; i<=numOrigins; i++){ eventsAll_vs_Q2[i] = new TH1F ("eventsAll_vs_Q2"+HistToName[i], "d#sigma/dQ^2. Origin: "+HistToName[i], kk_Q2, Q2_min, Q2_max); }
Float_t Q2BinWidth = eventsAll_vs_Q2[0]->GetBinWidth(0);        // bins width to be used in determining xsec as perweight/bin vs Q2

vector <TH1F *> events0pions_vs_Q2(numOrigins+1);
for (int i=0; i<=numOrigins; i++){ events0pions_vs_Q2[i] = new TH1F ("0pions_events_vs_Q2"+HistToName[i], "0#pi: d#sigma/dQ^2. Origin: "+HistToName[i], kk_Q2, Q2_min, Q2_max); }


// histograms versus W
Int_t kk_W=28;
Float_t W_min=1.0, W_max=2.9;
vector <TH1F *> eventsAll_vs_Wfree(numOrigins+1);
for (int i=0; i<=numOrigins; i++){eventsAll_vs_Wfree[i] = new TH1F ("eventsAll_vs_Wfree"+HistToName[i], "d#sigma/dW_{free}. Origin: "+HistToName[i], kk_W, W_min, W_max); }
Float_t WBinWidth = eventsAll_vs_Wfree[0]->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs W

vector <TH1F *> events0pions_vs_Wfree(numOrigins+1);
for (int i=0; i<=numOrigins; i++){events0pions_vs_Wfree[i] = new TH1F ("0pions_events_vs_Wfree"+HistToName[i], "0#pi: d#sigma/dW_{free}. Origin: "+HistToName[i], kk_W, W_min, W_max); }



// histograms versus muon angle
Int_t kk_theta=20;
Float_t theta_min=0.0, theta_max=40.00; //degrees
vector <TH1F *> eventsAll_vs_thetaMu(numOrigins+1);
for (int i=0; i<=numOrigins; i++){eventsAll_vs_thetaMu[i] = new TH1F ("eventsAll_vs_thetaMu"+HistToName[i], "d#sigma/d#theta_{#mu}. Origin: "+HistToName[i], kk_theta, theta_min, theta_max); }
Float_t thetaBinWidth = eventsAll_vs_thetaMu[0]->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs W

vector <TH1F *> events0pions_vs_thetaMu(numOrigins+1);
for (int i=0; i<=numOrigins; i++){events0pions_vs_thetaMu[i] = new TH1F ("0pions_events_vs_thetaMu"+HistToName[i], "0#pi: d#sigma/d#theta_{#mu}. Origin: "+HistToName[i], kk_theta, theta_min, theta_max); }


// histograms versus muon momentum
Int_t kk_p=25;
Float_t p_min=0.0, p_max=20.00; //GeV
vector <TH1F *> eventsAll_vs_pMu(numOrigins+1);
for (int i=0; i<=numOrigins; i++){eventsAll_vs_pMu[i] = new TH1F ("eventsAll_vs_pMu"+HistToName[i], "d#sigma/dp_{#mu}. Origin: "+HistToName[i], kk_p, p_min, p_max); }
Float_t pBinWidth = eventsAll_vs_pMu[0]->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs W

vector <TH1F *> events0pions_vs_pMu(numOrigins+1);
for (int i=0; i<=numOrigins; i++){events0pions_vs_pMu[i] = new TH1F ("events0pions_vs_pMu"+HistToName[i], "0#pi: d#sigma/dp_{#mu}. Origin: "+HistToName[i], kk_p, p_min, p_max); }



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


        	Float_t Q2, Wfree, thetamu, pmu;
        	Float_t emu, pxmu, pymu, pzmu; 


	        // ////////////  particle-loop-2  for recording muon
		for (Int_t  j=0; j<numberOutParticles; j++){
                if(particleID[j]==903){ emu=energy[j]; pxmu=momentumX[j]; pymu=momentumY[j]; pzmu=momentumZ[j];  }              // muon energy and momentum
                // if(particleID[j]==101 && charge[j]==1 && energy[j]>eH_max){eH_max=energy[j]; pxH_max=momentumX[j]; pyH_max=momentumY[j]; pzH_max=momentumZ[j]; }   //pion of max energy, its energy and momentum
                } // end particle-loop-2



                MuonKinematics b(enu, emu,pxmu,pymu,pzmu);
            	       Q2=b.Q2();
            	       eventsAll_vs_Q2[0]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsAll_vs_Q2[KToHist[origin]]->Fill(Q2, perweight/Q2BinWidth);
            	       Wfree=b.Wfree();
            	       eventsAll_vs_Wfree[0]->Fill(Wfree, perweight/WBinWidth);
            	       eventsAll_vs_Wfree[KToHist[origin]]->Fill(Wfree, perweight/WBinWidth);
            	       thetamu=b.thetamu()/TMath::Pi()*180.;
            	       eventsAll_vs_thetaMu[0]->Fill(thetamu, perweight/thetaBinWidth);
            	       eventsAll_vs_thetaMu[KToHist[origin]]->Fill(thetamu, perweight/thetaBinWidth);
            	       pmu=b.pmu();
            	       eventsAll_vs_pMu[0]->Fill(pmu, perweight/pBinWidth);
            	       eventsAll_vs_pMu[KToHist[origin]]->Fill(pmu, perweight/pBinWidth);


//###########################################################################   counting particles in each event



//cout<<"    origin="<<origin<<"    enu= "<<enu<<"     perweight="<<perweight<<endl;
//cout<<"numberOutParticles= "<<numberOutParticles<<endl;

	// particle-loop-1  for counting particles in each event 
	for (Int_t  j=0; j<numberOutParticles; j++){
	           //cout<<"\t\t\t j= "<<j<<"   particleId="<<particleID[j]<<endl;
	    //pions
	    if(particleID[j]==101 && perweight>0){ (numberOutPions[0])++;    (numberOutPions[(charge[j]+2)])++;     }
	    //nucleons
	    if(particleID[j]==1   && perweight>0){   (numberOutNucleons[0])++; (numberOutNucleons[(charge[j]+2)])++;  }
	} // end particle-loop-1



        //check that the pion and nucleon number is consistent
        if(numberOutPions[0]!=(numberOutPions[1]+numberOutPions[2]+numberOutPions[3]))
                        { cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
        if(numberOutNucleons[0]!=(numberOutNucleons[1]+numberOutNucleons[2]+numberOutNucleons[3]))
                        { cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl; }
        
       /*    
       if(jentry%10000==0) cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl;
       if(jentry%10000==0) cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl;
       */

//########################################################################  distributions AFTER particle counting


		if(numberOutPions[0]==0){
            	       events0pions_vs_Q2[0]->Fill(Q2, perweight/Q2BinWidth);
            	       events0pions_vs_Q2[KToHist[origin]]->Fill(Q2, perweight/Q2BinWidth);
            	       events0pions_vs_Wfree[0]->Fill(Wfree, perweight/WBinWidth);
            	       events0pions_vs_Wfree[KToHist[origin]]->Fill(Wfree, perweight/WBinWidth);
            	       events0pions_vs_thetaMu[0]->Fill(thetamu, perweight/thetaBinWidth);
            	       events0pions_vs_thetaMu[KToHist[origin]]->Fill(thetamu, perweight/thetaBinWidth);
            	       events0pions_vs_pMu[0]->Fill(pmu, perweight/pBinWidth);
            	       events0pions_vs_pMu[KToHist[origin]]->Fill(pmu, perweight/pBinWidth);
		}






   } // end loop-for-entries

//#######################
//####   LOOP OVER ENTRIES   IS OVER
//#######################

// file for saving all histograms
TFile * fileHists = new TFile("MuonKinematics.root","RECREATE");




///////////////////  stacked   Q2, Wfree,  pMu and thetaMu   plots 

THStack * eventsAll_vsQ2_origin = new THStack("Q2_distr_all","Q2-distr all");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_Q2[i]->SetFillColor(i+2);
                        eventsAll_vs_Q2[i]->SetMarkerStyle(21);
                        eventsAll_vs_Q2[i]->SetMarkerColor(i+2);
                        eventsAll_vs_Q2[i]->Write();
                        eventsAll_vsQ2_origin->Add(eventsAll_vs_Q2[i]); 
                       }


THStack * eventsAll_vsWfree_origin = new THStack("Wfree_distr_all","Wfree-distr all");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_Wfree[i]->SetFillColor(i+2);
                        eventsAll_vs_Wfree[i]->SetMarkerStyle(21);
                        eventsAll_vs_Wfree[i]->SetMarkerColor(i+2);
                        eventsAll_vs_Wfree[i]->Write();
                        eventsAll_vsWfree_origin->Add(eventsAll_vs_Wfree[i]); 
                       }


THStack * eventsAll_vspMu_origin = new THStack("muon_momentum_distr_all","muon-momentum distr all");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_pMu[i]->SetFillColor(i+2);
                        eventsAll_vs_pMu[i]->SetMarkerStyle(21);
                        eventsAll_vs_pMu[i]->SetMarkerColor(i+2);
                        eventsAll_vs_pMu[i]->Write();
                        eventsAll_vspMu_origin->Add(eventsAll_vs_pMu[i]); 
                       }


THStack * eventsAll_vsthetaMu_origin = new THStack("muon_angle_distr_all","muon-angle distr all");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_thetaMu[i]->SetFillColor(i+2);
                        eventsAll_vs_thetaMu[i]->SetMarkerStyle(21);
                        eventsAll_vs_thetaMu[i]->SetMarkerColor(i+2);
                        eventsAll_vs_thetaMu[i]->Write();
                        eventsAll_vsthetaMu_origin->Add(eventsAll_vs_thetaMu[i]); 
                       }






///////////////////////////////   reaction title and final state
TLatex * reac = new TLatex();
reac->SetTextAlign(33);
reac->SetNDC();
TLatex * fs = new TLatex();
fs->SetTextAlign(33);
fs->SetNDC();
TString fs_0pions = "0#pi";
TString fs_1piplus = "1#pi^{+} 0#pi^{0} 0#pi^{-} ";





/////////////  legends 

TLegend *legendOrigin = new TLegend(0.7, 0.4, 0.92, 0.7);
legendOrigin->SetBorderSize(0);
for(int i=1; i<=6; i++){ legendOrigin->AddEntry(eventsAll_vs_Q2[i], HistToName[i], "F");   }





///////////////////////   drawing muon-kinematic  xsec for all events 


TCanvas * xsec_origin = new TCanvas("xsec_origin", "xsec: various origins");

xsec_origin->Divide(2,2);


xsec_origin->cd(1);
eventsAll_vsQ2_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
//fs->DrawLatex(0.91,0.8,fs_0pions);
xsec_origin->Update();



xsec_origin->cd(2);
eventsAll_vspMu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
//fs->DrawLatex(0.91,0.8,fs_0pions);
xsec_origin->Update();
// exercises on setttin user range 
//THStack * eventsAll_vsQ2_origin_copy = (THStack *)eventsAll_vsQ2_origin->Clone();
//eventsAll_vsQ2_origin_copy->SetName("all copy");
//eventsAll_vsQ2_origin_copy->Draw();
//eventsAll_vsQ2_origin_copy->GetXaxis()->SetRangeUser(0.2,2.0);
//eventsAll_vsQ2_origin_copy->GetYaxis()->SetRangeUser(0.0,10000.0);


xsec_origin->cd(3);
eventsAll_vsWfree_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
//fs->DrawLatex(0.91,0.8,fs_0pions);


xsec_origin->cd(4);
eventsAll_vsthetaMu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
//fs->DrawLatex(0.91,0.8,fs_0pions);


xsec_origin->Update();
//xsec_eC_piplus_origin->Write();
xsec_origin->Print("xsec_origin.eps");
xsec_origin->Write();

/////////////////////
// Adding ArgoNeut data for dsigma/dthetamu and dsigma/dpMu plots
Int_t A=40;

//TCanvas * ArgoNeuT_xsec_origin = new TCanvas("ArgoNeut_xsec_origin", "ArgoNeut xsec: various origins");
//ArgoNeuT_xsec_origin = (TCanvas*)xsec_origin->DrawClone();

TCanvas * ArgoNeuT_xsec_origin = (TCanvas*)xsec_origin->DrawClone();
ArgoNeuT_xsec_origin->SetName("ArgoNeut_xsec_origin");
ArgoNeuT_xsec_origin->SetTitle("ArgoNeut xsec: various origins");
ArgoNeuT_xsec_origin->Update();


ArgoNeuT_xsec_origin->cd(2);
ArgoNeuT_xsec_origin->SetGrid();
TGraphErrors * ArgoNeut_dpMu = new TGraphErrors("events/ArgoNeuT-numu-numode-dsidpMu.exp","%lg %lg %lg %lg");
// divide data points be A=40
for (int i=0;i<ArgoNeut_dpMu->GetN();i++){ ArgoNeut_dpMu->GetY()[i] *= 1./A; ArgoNeut_dpMu->GetEY()[i] *= 1./A; }
//
ArgoNeut_dpMu->SetTitle("ArgoNeut muon momentum ditribution");
ArgoNeut_dpMu->SetMarkerStyle(kFullCircle);
ArgoNeut_dpMu->SetMarkerColor(kBlack);
ArgoNeut_dpMu->Draw("PSame");


ArgoNeuT_xsec_origin->Modified();
ArgoNeuT_xsec_origin->Update();




ArgoNeuT_xsec_origin->cd(4);
ArgoNeuT_xsec_origin->SetGrid();
TGraphErrors * ArgoNeuT_dthetaMu = new TGraphErrors("events/ArgoNeuT-numu-numode-dsidthetaMu.exp","%lg %lg %lg %lg");
// divide data points by A=40
for (int i=0;i<ArgoNeuT_dthetaMu->GetN();i++){ ArgoNeuT_dthetaMu->GetY()[i] *= 1./A; ArgoNeuT_dthetaMu->GetEY()[i] *= 1./A; }
//
ArgoNeuT_dthetaMu->SetTitle("ArgoNeut muon angle distribution");
ArgoNeuT_dthetaMu->SetMarkerStyle(kFullCircle);
ArgoNeuT_dthetaMu->SetMarkerColor(kBlack);
ArgoNeuT_dthetaMu->SetMaximum(0.18);
ArgoNeuT_dthetaMu->Draw("PSame");
ArgoNeuT_xsec_origin->Update();

//double ymax = TMath::MinElement(ArgoNeuT_dthetaMu->GetN(),ArgoNeuT_dthetaMu->GetY()); 
//ArgoNeuT_dthetaMu->GetHistogram()->SetMaximum(ymax);

ArgoNeuT_xsec_origin->Update();
ArgoNeuT_xsec_origin->Print("ArgoNeut_xsec_origin.eps");
ArgoNeuT_xsec_origin->Write();















////////////////////////////////////////////////////////////////////////////   0pions  events  
///////////////////  stacked   Q2, Wfree,  pMu and thetaMu   plots   

THStack * events0pions_vsQ2_origin = new THStack("Q2_distr_all","Q2-distr all");
for(int i=1; i<=6; i++){events0pions_vs_Q2[i]->SetFillColor(i+2);
                        events0pions_vs_Q2[i]->SetMarkerStyle(21);
                        events0pions_vs_Q2[i]->SetMarkerColor(i+2);
                        events0pions_vs_Q2[i]->Write();
                        events0pions_vsQ2_origin->Add(events0pions_vs_Q2[i]); 
                       }


THStack * events0pions_vsWfree_origin = new THStack("Wfree_distr_all","Wfree-distr all");
for(int i=1; i<=6; i++){events0pions_vs_Wfree[i]->SetFillColor(i+2);
                        events0pions_vs_Wfree[i]->SetMarkerStyle(21);
                        events0pions_vs_Wfree[i]->SetMarkerColor(i+2);
                        events0pions_vs_Wfree[i]->Write();
                        events0pions_vsWfree_origin->Add(events0pions_vs_Wfree[i]); 
                       }


THStack * events0pions_vspMu_origin = new THStack("muon_momentum_distr_all","muon-momentum distr all");
for(int i=1; i<=6; i++){events0pions_vs_pMu[i]->SetFillColor(i+2);
                        events0pions_vs_pMu[i]->SetMarkerStyle(21);
                        events0pions_vs_pMu[i]->SetMarkerColor(i+2);
                        events0pions_vs_pMu[i]->Write();
                        events0pions_vspMu_origin->Add(events0pions_vs_pMu[i]); 
                       }


THStack * events0pions_vsthetaMu_origin = new THStack("muon_angle_distr_all","muon-angle distr all");
for(int i=1; i<=6; i++){events0pions_vs_thetaMu[i]->SetFillColor(i+2);
                        events0pions_vs_thetaMu[i]->SetMarkerStyle(21);
                        events0pions_vs_thetaMu[i]->SetMarkerColor(i+2);
                        events0pions_vs_thetaMu[i]->Write();
                        events0pions_vsthetaMu_origin->Add(events0pions_vs_thetaMu[i]); 
                       }





///////////////////////   drawing muon-kinematic  xsec for events with 0 pions


TCanvas * xsec_origin_0pions = new TCanvas("0pions_xsec_origin", "0#pi: xsec: various origins");

xsec_origin_0pions->Divide(2,2);


xsec_origin_0pions->cd(1);
events0pions_vsQ2_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
xsec_origin_0pions->Update();



xsec_origin_0pions->cd(2);
events0pions_vspMu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
xsec_origin_0pions->Update();


xsec_origin_0pions->cd(3);
events0pions_vsWfree_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


xsec_origin_0pions->cd(4);
events0pions_vsthetaMu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


xsec_origin_0pions->Update();
//xsec_eC_piplus_origin->Write();
xsec_origin_0pions->Print("0pions_xsec_origin.eps");
xsec_origin_0pions->Write();





////////////////
// draw the same canvas with muon kinematics for 0 pion events, but now with experimental points from ArgoNeuT barnumu barnumode

TCanvas * ArgoNeuT_xsec_origin_0pions = new TCanvas("0pions_ArgoNeuT_xsec_origin", "0#pi: ArgoNeuT xsec: various origins");

ArgoNeuT_xsec_origin_0pions->Divide(2,2);


ArgoNeuT_xsec_origin_0pions->cd(1);
events0pions_vsQ2_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
ArgoNeuT_xsec_origin_0pions->Update();



ArgoNeuT_xsec_origin_0pions->cd(2);

TGraphErrors * ArgoNeut_dpMu_barnumu_barnumode = new TGraphErrors("events/ArgoNeuT-barnumu-barnumode-0pions-dsidpMu.exp","%lg %lg %lg %lg");
// divide data points be A=40
Float_t scale=1./35.; //scale from events to xsec
for (int i=0;i<ArgoNeut_dpMu_barnumu_barnumode->GetN();i++){ ArgoNeut_dpMu_barnumu_barnumode->GetY()[i] *= scale/A; 
                                                             ArgoNeut_dpMu_barnumu_barnumode->GetEY()[i] *= scale/A; }
//
events0pions_vspMu_origin->Draw();
//
ArgoNeut_dpMu_barnumu_barnumode->SetTitle("ArgoNeut muon momentum ditribution");
ArgoNeut_dpMu_barnumu_barnumode->SetMarkerStyle(kFullCircle);
ArgoNeut_dpMu_barnumu_barnumode->SetMarkerColor(kBlack);
ArgoNeut_dpMu_barnumu_barnumode->Draw("PSame");
//
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

ArgoNeuT_xsec_origin_0pions->Update();






ArgoNeuT_xsec_origin_0pions->cd(3);
events0pions_vsWfree_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_xsec_origin_0pions->cd(4);
events0pions_vsthetaMu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_xsec_origin_0pions->Update();
//xsec_eC_piplus_origin->Write();
ArgoNeuT_xsec_origin_0pions->Print("0pions_ArgoNeuT_xsec_origin.eps");
ArgoNeuT_xsec_origin_0pions->Write();





// saving histograms to a file
// does NOT work for unknown reason
//fileHists->Write();
fileHists->Close();



}   // end Analyze::LoopMuonKinematics


// change range of the graph
//  mygraph->GetHistogram()->SetMaximum(mymaximum);
//  mygraph->GetHistogram()->SetMinimum(myMinimum);

