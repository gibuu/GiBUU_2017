#include "PionsInCLAS.h"

void Analyze::LoopPionsInCLAS()
{

// histograms versus Q2
Int_t kk_Q2=40;
Double_t Q2_min=1.0, Q2_max=5.0;
vector <TH1F *> eventsMultiPiplus_vs_Q2(7);
for (int i=0; i<=6; i++){ eventsMultiPiplus_vs_Q2[i] = new TH1F ("eventsMultiPiplus_vs_Q2"+HistToName[i], "dsi/dQ2 for multi-pi+ production in CLAS. Origin: "+HistToName[i], kk_Q2, Q2_min, Q2_max); }
Double_t Q2BinWidth = eventsMultiPiplus_vs_Q2[0]->GetBinWidth(0);        // bins width to be used in determining xsec as perweight/bin vs Q2
vector <TH1F *> eventsMultiPiminus_vs_Q2(7);
for (int i=0; i<=6; i++){ eventsMultiPiminus_vs_Q2[i] = new TH1F ("eventsMultiPiminus_vs_Q2"+HistToName[i], "dsi/dQ2 for multi-pi- production in CLAS. Origin: "+HistToName[i], kk_Q2, Q2_min, Q2_max); }


// histograms versus W
Int_t kk_W=28;
Double_t W_min=1.0, W_max=2.9;
vector <TH1F *> eventsMultiPiplus_vs_Wfree(7);
for (int i=0; i<=6; i++){eventsMultiPiplus_vs_Wfree[i] = new TH1F ("eventsMultiPiplus_vs_Wfree"+HistToName[i], "dsi/dWfree for multi-pi+ production in CLAS. Origin: "+HistToName[i], kk_W, W_min, W_max); }
Double_t WBinWidth = eventsMultiPiplus_vs_Wfree[0]->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs W
vector <TH1F *> eventsMultiPiminus_vs_Wfree(7);
for (int i=0; i<=6; i++){eventsMultiPiminus_vs_Wfree[i] = new TH1F ("eventsMultiPiminus_vs_Wfree"+HistToName[i], "dsi/dWfree for multi-pi- production in CLAS. Origin: "+HistToName[i], kk_W, W_min, W_max); }
//vector <TH1F *> eventsAll_vs_Wfree(7);
//for (int i=0; i<=6; i++){eventsAll_vs_Wfree[i] = new TH1F ("eventsAll_vs_Wfree", "dsi/dWfree. Origin: "+HistToName[i], kk_W, W_min, W_max); }





// histograms versus pionMomentum
Int_t kk_ppi=38;
Double_t ppi_min=0.2, ppi_max=4.0;
TH1 * eventsMultiPiplus_vs_ppi_max = new TH1D ("eventsMultiPiplus_vs_ppi_max",   "dsi/dppi for multi-pi+ (pion of max energy) production in CLAS", kk_ppi, ppi_min, ppi_max);  
Double_t ppiBinWidth = eventsMultiPiplus_vs_ppi_max->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs ppi
TH1 * eventsMultiPiminus_vs_ppi_max = new TH1D ("eventsMultiPiminus_vs_ppi_max", "dsi/dppi for multi-pi- (pion of max energy) production in CLAS", kk_ppi, ppi_min, ppi_max);  


// histograms versus pionAngle
Int_t kk_thetapi=10;
Double_t thetapi_min=5.0, thetapi_max=55.0;
TH1 * eventsMultiPiplus_vs_thetapi_max = new TH1D ("eventsMultiPiplus_vs_thetapi_max",   "dsi/dthetapi for multi-pi+ (pion of max energy) production in CLAS", kk_thetapi, thetapi_min, thetapi_max);  
Double_t thetapiBinWidth = eventsMultiPiplus_vs_thetapi_max->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs ppi
TH1 * eventsMultiPiminus_vs_thetapi_max = new TH1D ("eventsMultiPiminus_vs_thetapi_max", "dsi/dthetapi for multi-pi- (pion of max energy) production in CLAS", kk_thetapi, thetapi_min, thetapi_max);  







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



        //check that the pion and nucleon number is consistent
        if(numberOutPions[0]!=(numberOutPions[1]+numberOutPions[2]+numberOutPions[3]))
                        { cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
        if(numberOutNucleons[0]!=(numberOutNucleons[1]+numberOutNucleons[2]+numberOutNucleons[3]))
                        { cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl; }





//########################################################################  distributions AFTER particle counting


        // multi-pi+ and  multi-pi- events as selected in CLAS analysis
        if(  numberOutPions[3]>=1 || numberOutPions[1]>=1){
        	    
        	Float_t Q2, Wfree, ppi, thetapi;
        	Float_t emu, pxmu, pymu, pzmu; 
        	Float_t eH_max=0., pxH_max=0., pyH_max=0., pzH_max=0.;
        	Int_t	charge_max=0;


	        // ////////////  particle-loop-2  for recording muon and recording pion of maximum energy
		for (Int_t  j=0; j<numberOutParticles; j++){
                if(particleID[j]==903){ emu=energy[j]; pxmu=momentumX[j]; pymu=momentumY[j]; pzmu=momentumZ[j];  }              // muon energy and momentum
                if(particleID[j]==101 && ( charge[j]==1 || charge[j]==-1)  && energy[j]>eH_max){
                          charge_max=charge[j];
                          eH_max=energy[j];     pxH_max=momentumX[j]; pyH_max=momentumY[j]; pzH_max=momentumZ[j]; 
                                                                                               }   //pion of max energy, its charge, energy and momentum
                } // end particle-loop-2



                PionsInCLAS b(enu, emu,pxmu,pymu,pzmu,eH_max,pxH_max,pyH_max,pzH_max);
                if (b.SelectOutElectronY() && b.SelectQ2()  && b.SelectOutElectronAngle()   && b.SelectOutPionAngle() && b.SelectOutPionMomentum() ){
            	       Q2=b.Q2();
            	       Wfree=b.Wfree();
            	       ppi=b.ppi();
            	       thetapi=b.thetapi()*180./TMath::Pi();
            	       
            	       if (charge_max==1){
            	       eventsMultiPiplus_vs_Q2[0]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsMultiPiplus_vs_Q2[KToHist[origin]]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsMultiPiplus_vs_Wfree[0]->Fill(Wfree, perweight/WBinWidth);
            	       eventsMultiPiplus_vs_Wfree[KToHist[origin]]->Fill(Wfree, perweight/WBinWidth);
            	       eventsMultiPiplus_vs_ppi_max->Fill(ppi, perweight/ppiBinWidth);
            	       eventsMultiPiplus_vs_thetapi_max->Fill(thetapi, perweight/thetapiBinWidth);
                                         }
                       else if (charge_max==-1){
            	       eventsMultiPiminus_vs_Q2[0]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsMultiPiminus_vs_Q2[KToHist[origin]]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsMultiPiminus_vs_Wfree[0]->Fill(Wfree, perweight/WBinWidth);
            	       eventsMultiPiminus_vs_Wfree[KToHist[origin]]->Fill(Wfree, perweight/WBinWidth);
            	       eventsMultiPiminus_vs_ppi_max->Fill(ppi, perweight/ppiBinWidth);
            	       eventsMultiPiminus_vs_thetapi_max->Fill(thetapi, perweight/thetapiBinWidth);
					       }
            	}
               // cout<<"in Analyze.C:  jentry= "<<jentry<< "        enuRec1Piplus = "<< enuRec1Piplus<<endl;
        }  // end multi-pi+  and multi-pi- events as selected in CLAS analysis








   } // end loop-for-entries

//#######################
//####   LOOP OVER ENTRIES   IS OVER
//#######################


TCanvas * xsec_eC_piplus = new TCanvas("xsec_eC_piplus", "C(12)  multi-pi+ CLAS");
xsec_eC_piplus->Divide(2,2);

xsec_eC_piplus->cd(1);

eventsMultiPiplus_vs_Q2[0]->SetLineColor(2); //2=red
eventsMultiPiplus_vs_Q2[0]->SetLineWidth(2);
eventsMultiPiplus_vs_Q2[0]->GetXaxis()->SetTitle("Q2");
eventsMultiPiplus_vs_Q2[0]->GetYaxis()->SetTitle("multi-pi+, pion of highest energy");
eventsMultiPiplus_vs_Q2[0]->Draw();

/*
events1Piplus_vs_EnuRec->SetLineColor(3); //3=green
events1Piplus_vs_EnuRec->SetLineWidth(2);
events1Piplus_vs_EnuRec->Draw("same");

events1Piplus_vs_EnuRecD->SetLineColor(4); //4=blue 5=yellow 6=magenta
events1Piplus_vs_EnuRecD->SetLineWidth(2);
events1Piplus_vs_EnuRecD->Draw("same");



TLegend *legend4 = new TLegend(0.5, 0.6, 0.85, 0.8);
legend4->SetBorderSize(0);
legend4->AddEntry(events1Piplus_vs_EnuTrue, "true", "L");
legend4->AddEntry(events1Piplus_vs_EnuRec,  "reconstructed via pion and muon", "L");
legend4->AddEntry(events1Piplus_vs_EnuRecD,  "reconstructed via Delta", "L");
legend4->Draw();
*/

xsec_eC_piplus->cd(2);

eventsMultiPiplus_vs_Wfree[0]->SetLineColor(2); //red
eventsMultiPiplus_vs_Wfree[0]->SetLineWidth(2);
eventsMultiPiplus_vs_Wfree[0]->GetXaxis()->SetTitle("invariant mass");
eventsMultiPiplus_vs_Wfree[0]->GetYaxis()->SetTitle("multi-pi+, pion of highest energy");
eventsMultiPiplus_vs_Wfree[0]->Draw();

/*
events1Pi0_vs_EnuRec->SetLineColor(3); //green
events1Pi0_vs_EnuRec->SetLineWidth(2);
events1Pi0_vs_EnuRec->Draw("same");

events1Pi0_vs_EnuRecD->SetLineColor(4); //blue
events1Pi0_vs_EnuRecD->SetLineWidth(2);
events1Pi0_vs_EnuRecD->Draw("same");



TLegend *legend5 = new TLegend(0.5, 0.6, 0.85, 0.8);
legend5->SetBorderSize(0);
legend5->AddEntry(events1Piplus_vs_EnuTrue, "true", "L");
legend5->AddEntry(events1Piplus_vs_EnuRec,  "reconstructed via pion and muon", "L");
legend5->AddEntry(events1Piplus_vs_EnuRecD,  "reconstructed via Delta", "L");
legend5->Draw();
*/



xsec_eC_piplus->cd(3);

eventsMultiPiplus_vs_ppi_max->SetLineColor(2); //red
eventsMultiPiplus_vs_ppi_max->SetLineWidth(2);
eventsMultiPiplus_vs_ppi_max->GetXaxis()->SetTitle("pion momentum");
eventsMultiPiplus_vs_ppi_max->GetYaxis()->SetTitle("multi-pi+, pion of highest energy");
eventsMultiPiplus_vs_ppi_max->Draw();



xsec_eC_piplus->cd(4);

eventsMultiPiplus_vs_thetapi_max->SetLineColor(2); //red
eventsMultiPiplus_vs_thetapi_max->SetLineWidth(2);
eventsMultiPiplus_vs_thetapi_max->GetXaxis()->SetTitle("pion angle");
eventsMultiPiplus_vs_thetapi_max->GetYaxis()->SetTitle("multi-pi+, pion of highest energy");
eventsMultiPiplus_vs_thetapi_max->Draw();


xsec_eC_piplus->Update();
xsec_eC_piplus->Print("xsec_eC_piplus.eps");










///////////////////  stacked pi+

THStack * eventsMultiPiplus_vsQ2_origin = new THStack("multi-pi+","origin of events in Q2-distribution of multi-pi+ events");
for(int i=1; i<=5; i++){eventsMultiPiplus_vs_Q2[i]->SetFillColor(i+2);
                        eventsMultiPiplus_vs_Q2[i]->SetMarkerStyle(21);
                        eventsMultiPiplus_vs_Q2[i]->SetMarkerColor(i+2);
                        eventsMultiPiplus_vsQ2_origin->Add(eventsMultiPiplus_vs_Q2[i]); 
                       }


THStack * eventsMultiPiplus_vsWfree_origin = new THStack("multi-pi+","origin of events in Wfree-distribution of multi-pi+ events");
for(int i=1; i<=5; i++){eventsMultiPiplus_vs_Wfree[i]->SetFillColor(i+2);
                        eventsMultiPiplus_vs_Wfree[i]->SetMarkerStyle(21);
                        eventsMultiPiplus_vs_Wfree[i]->SetMarkerColor(i+2);
                        eventsMultiPiplus_vsWfree_origin->Add(eventsMultiPiplus_vs_Wfree[i]); 
                       }


TCanvas * xsec_eC_piplus_origin = new TCanvas("xsec_eC_piplus_origin", "C(12)  multi-pi+ CLAS");

xsec_eC_piplus_origin->Divide(2,2);

xsec_eC_piplus_origin->cd(1);
eventsMultiPiplus_vs_Q2[0]->Draw();

xsec_eC_piplus_origin->cd(2);
eventsMultiPiplus_vsQ2_origin->Draw();
TLegend *legendOrigin = new TLegend(0.5, 0.6, 0.75, 0.8);
legendOrigin->SetBorderSize(0);
for(int i=2; i<=5; i++){ legendOrigin->AddEntry(eventsMultiPiplus_vs_Q2[i], HistToName[i], "F");   }
legendOrigin->Draw();
xsec_eC_piplus_origin->Update();



xsec_eC_piplus_origin->cd(3);
eventsMultiPiplus_vs_Wfree[0]->Draw();

xsec_eC_piplus_origin->cd(4);
eventsMultiPiplus_vsWfree_origin->Draw();
legendOrigin->Draw();

xsec_eC_piplus_origin->Update();
//xsec_eC_piplus_origin->Write();
xsec_eC_piplus_origin->Print("xsec_eC_piplus_origin.eps");













///////////////////  stacked pi-

THStack * eventsMultiPiminus_vsQ2_origin = new THStack("multi-pi-","origin of events in Q2-distribution of multi-pi- events");
for(int i=1; i<=5; i++){eventsMultiPiminus_vs_Q2[i]->SetFillColor(i+2);
                        eventsMultiPiminus_vs_Q2[i]->SetMarkerStyle(21);
                        eventsMultiPiminus_vs_Q2[i]->SetMarkerColor(i+2);
                        eventsMultiPiminus_vsQ2_origin->Add(eventsMultiPiminus_vs_Q2[i]); 
                       }


THStack * eventsMultiPiminus_vsWfree_origin = new THStack("multi-pi-","origin of events in Wfree-distribution of multi-pi- events");
for(int i=1; i<=5; i++){eventsMultiPiminus_vs_Wfree[i]->SetFillColor(i+2);
                        eventsMultiPiminus_vs_Wfree[i]->SetMarkerStyle(21);
                        eventsMultiPiminus_vs_Wfree[i]->SetMarkerColor(i+2);
                        eventsMultiPiminus_vsWfree_origin->Add(eventsMultiPiminus_vs_Wfree[i]); 
                       }


TCanvas * xsec_eC_piminus_origin = new TCanvas("xsec_eC_piminus_origin", "C(12)  multi-pi- CLAS");

xsec_eC_piminus_origin->Divide(2,2);

xsec_eC_piminus_origin->cd(1);
eventsMultiPiminus_vs_Q2[0]->Draw();

xsec_eC_piminus_origin->cd(2);
eventsMultiPiminus_vsQ2_origin->Draw();

legendOrigin->Clear();
for(int i=2; i<=5; i++){ legendOrigin->AddEntry(eventsMultiPiminus_vs_Q2[i], HistToName[i], "F");   }

legendOrigin->Draw();
xsec_eC_piminus_origin->Update();



xsec_eC_piminus_origin->cd(3);
eventsMultiPiminus_vs_Wfree[0]->Draw();

xsec_eC_piminus_origin->cd(4);
eventsMultiPiminus_vsWfree_origin->Draw();
legendOrigin->Draw();

xsec_eC_piminus_origin->Update();
//xsec_eC_piminus_origin->Write();
xsec_eC_piminus_origin->Print("xsec_eC_piminus_origin.eps");
















///////////////////   pi- 

TCanvas * xsec_eC_piminus = new TCanvas("xsec_eC_piminus", "C(12)  multi-pi- CLAS");
xsec_eC_piminus->Divide(2,2);

xsec_eC_piminus->cd(1);

eventsMultiPiminus_vs_Q2[0]->SetLineColor(2); //2=red
eventsMultiPiminus_vs_Q2[0]->SetLineWidth(2);
eventsMultiPiminus_vs_Q2[0]->GetXaxis()->SetTitle("Q2");
eventsMultiPiminus_vs_Q2[0]->GetYaxis()->SetTitle("multi-pi-, pion of highest energy");
eventsMultiPiminus_vs_Q2[0]->Draw();



xsec_eC_piminus->cd(2);

eventsMultiPiminus_vs_Wfree[0]->SetLineColor(2); //red
eventsMultiPiminus_vs_Wfree[0]->SetLineWidth(2);
eventsMultiPiminus_vs_Wfree[0]->GetXaxis()->SetTitle("invariant mass");
eventsMultiPiminus_vs_Wfree[0]->GetYaxis()->SetTitle("multi-pi-, pion of highest energy");
eventsMultiPiminus_vs_Wfree[0]->Draw();



xsec_eC_piminus->cd(3);

eventsMultiPiminus_vs_ppi_max->SetLineColor(2); //red
eventsMultiPiminus_vs_ppi_max->SetLineWidth(2);
eventsMultiPiminus_vs_ppi_max->GetXaxis()->SetTitle("pion momentum");
eventsMultiPiminus_vs_ppi_max->GetYaxis()->SetTitle("multi-pi-, pion of highest energy");
eventsMultiPiminus_vs_ppi_max->Draw();



xsec_eC_piminus->cd(4);

eventsMultiPiminus_vs_thetapi_max->SetLineColor(2); //red
eventsMultiPiminus_vs_thetapi_max->SetLineWidth(2);
eventsMultiPiminus_vs_thetapi_max->GetXaxis()->SetTitle("pion angle");
eventsMultiPiminus_vs_thetapi_max->GetYaxis()->SetTitle("multi-pi-, pion of highest energy");
eventsMultiPiminus_vs_thetapi_max->Draw();


xsec_eC_piminus->Update();
xsec_eC_piminus->Print("xsec_eC_piminus.eps");





}   // end Analyze::LoopPionsInCLAS)












