#include "EnergyReconstruction.h"
#include "ConstantsGiBUU.h"

void Analyze::LoopEnergyReconstruction(TString reaction)
{


const Float_t thresEkinN=0.021;//0.21 for LAr experiments
const Short_t numOrigins=6;

// histograms versus neutrino energy
Int_t kk_Enu=20;
Double_t Enu_min=0.0, Enu_max=16.0;

TH1F * events1Piplus_vs_EnuTrue = new TH1F ("events1Piplus_vs_EnuTrue", "event distribution 1pi+ production vs true neutrino energy", kk_Enu, Enu_min, Enu_max);  
Double_t EnergyBinWidth = events1Piplus_vs_EnuTrue->GetBinWidth(0);                      // bins width to be used in determining xsec as perweight/bin vs enu
TH1F * events1Piplus_vs_EnuRec  = new TH1F ("events1Piplus_vs_EnuRec",  "event distribution 1pi+ production vs reconstructed neutrino energy", kk_Enu, Enu_min, Enu_max);
TH1F * events1Piplus_vs_EnuRecD  = new TH1F ("events1Piplus_vs_EnuRecD",  "event distribution 1pi+ production vs reconstructed neutrino energy assuming Delta on-shell", kk_Enu, Enu_min, Enu_max);


TH1F * events1Pi0_vs_EnuTrue = new TH1F ("events1Pi0_vs_EnuTrue", "event distribution 1pi0 production vs true neutrino energy", kk_Enu, Enu_min, Enu_max); 
TH1F * events1Pi0_vs_EnuRec  = new TH1F ("events1Pi0_vs_EnuRec",  "event distribution 1pi0 production vs reconstructed neutrino energy", kk_Enu, Enu_min, Enu_max);
TH1F * events1Pi0_vs_EnuRecD  = new TH1F ("events1Pi0_vs_EnuRecD",  "event distribution 1pi0 production vs reconstructed neutrino energy assuming Delta on-shell", kk_Enu, Enu_min, Enu_max);


Int_t numMaxProtons=3;
char name[30], title[60];
TH1F * events0pions_vs_EnuTrue[numMaxProtons+1][numOrigins+1];
for (int outP=0; outP<=numMaxProtons; outP++){
for (int i=0; i<=numOrigins; i++){
    sprintf(name,"0pions_%dprotons_vs_Enu_%s",outP,HistToName[i].Data());
    sprintf(title,"0pions %dprotons: events vs E_{#nu} %s",outP,HistToName[i].Data());
    events0pions_vs_EnuTrue[outP][i] = new TH1F (name,title,kk_Enu, Enu_min, Enu_max); }
}
TH1F * events0pions_Allp_vs_EnuTrue[numOrigins+1];
for (int i=0; i<=numOrigins; i++){
    sprintf(name,"0pions_Xprotons_vs_Enu_%s",HistToName[i].Data());
    sprintf(title,"0pions Xprotons: events vs E_{#nu} %s",HistToName[i].Data());
    events0pions_Allp_vs_EnuTrue[i] = new TH1F (name,title, kk_Enu, Enu_min, Enu_max); 
}


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
	    if(particleID[j]==1 && (energy[j]-m_N)>thresEkinN){   (numberOutNucleons[0])++; (numberOutNucleons[(charge[j]+2)])++;  }
	} // end particle-loop-1



        //check that the pion and nucleon number is consistent
        if(numberOutPions[0]!=(numberOutPions[1]+numberOutPions[2]+numberOutPions[3]))
                        { cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
        if(numberOutNucleons[0]!=(numberOutNucleons[1]+numberOutNucleons[2]+numberOutNucleons[3]))
                        { cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl; }





//########################################################################  distributions AFTER particle counting


        // 1pi+ events
        if(numberOutPions[0]==1 && numberOutPions[3]==1){
        	events1Piplus_vs_EnuTrue->Fill(enu, perweight/EnergyBinWidth );
        	    
        	Float_t enuRec1Piplus, enuRec1PiplusD;
        	Float_t emu, pxmu, pymu, pzmu, eH, pxH, pyH, pzH;


	        // ////////////  particle-loop-2  for energy reconstruction 
		for (Int_t  j=0; j<numberOutParticles; j++){
                if(particleID[j]==903){ emu=energy[j]; pxmu=momentumX[j]; pymu=momentumY[j]; pzmu=momentumZ[j];  }              // muon energy and momentum
                if(particleID[j]==101 && charge[j]==1){eH=energy[j]; pxH=momentumX[j]; pyH=momentumY[j]; pzH=momentumZ[j]; }   //pion energy and momentum
                } // end particle-loop-2


                EnergyReconstruction b(emu,pxmu,pymu,pzmu,eH,pxH,pyH,pzH);
                enuRec1Piplus=b.OnePionEnergyReconstruction();
                enuRec1PiplusD=b.DeltaEnergyReconstruction();
//                cout<<"in Analyze.C:  jentry= "<<jentry<< "        enuRec1Piplus = "<< enuRec1Piplus<<endl;
                events1Piplus_vs_EnuRec->Fill(enuRec1Piplus, perweight/EnergyBinWidth );
                events1Piplus_vs_EnuRecD->Fill(enuRec1PiplusD, perweight/EnergyBinWidth );
        }  // end 1pi+ events









        // 1pi0 events
        if(numberOutPions[0]==1 && numberOutPions[2]==1){
        	events1Pi0_vs_EnuTrue->Fill(enu, perweight/EnergyBinWidth );
        	    
        	Float_t enuRec1Pi0, enuRec1Pi0D;
        	Float_t emu, pxmu, pymu, pzmu, eH, pxH, pyH, pzH;


	        // ////////////  particle-loop-2  for energy reconstruction 
		for (Int_t  j=0; j<numberOutParticles; j++){
                if(particleID[j]==903){ emu=energy[j]; pxmu=momentumX[j]; pymu=momentumY[j]; pzmu=momentumZ[j];  }              // muon energy and momentum
                if(particleID[j]==101 && charge[j]==0){eH=energy[j]; pxH=momentumX[j]; pyH=momentumY[j]; pzH=momentumZ[j]; }   //pion energy and momentum
                } // end particle-loop-2


                EnergyReconstruction b(emu,pxmu,pymu,pzmu,eH,pxH,pyH,pzH);
                enuRec1Pi0=b.OnePionEnergyReconstruction();
                enuRec1Pi0D=b.DeltaEnergyReconstruction();
//                cout<<"in Analyze.C:  jentry= "<<jentry<< "        enuRec1Pi0 = "<< enuRec1Pi0<<endl;
                events1Pi0_vs_EnuRec->Fill(enuRec1Pi0, perweight/EnergyBinWidth );
                events1Pi0_vs_EnuRecD->Fill(enuRec1Pi0D, perweight/EnergyBinWidth );
        }  // end 1pi0 events




        // 0pion events
        if(numberOutPions[0]==0 && numberOutNucleons[3]<=numMaxProtons){
        	events0pions_vs_EnuTrue[numberOutNucleons[3]][0]->Fill(enu, perweight/EnergyBinWidth );
        	events0pions_vs_EnuTrue[numberOutNucleons[3]][KToHist[origin]]->Fill(enu, perweight/EnergyBinWidth );
        	events0pions_Allp_vs_EnuTrue[0]->Fill(enu, perweight/EnergyBinWidth );
        	events0pions_Allp_vs_EnuTrue[KToHist[origin]]->Fill(enu, perweight/EnergyBinWidth );
        	    
        }  // end 0pion events









   } // end loop-for-entries

//#######################
//####   LOOP OVER ENTRIES   IS OVER
//#######################


// file for saving all histograms
TFile * fileHists = new TFile("EnergyReconstruction.root","RECREATE");


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
for(int i=1; i<=6; i++){ legendOrigin->AddEntry(events0pions_Allp_vs_EnuTrue[i], HistToName[i], "F");   }










TCanvas * energyReconstruction_1piplus = new TCanvas("energyReconstruction_1piplus", "1#pi^{+} event distributions");
energyReconstruction_1piplus->Divide(1,2);

energyReconstruction_1piplus->cd(1);

events1Piplus_vs_EnuTrue->SetLineColor(2); //2=red
events1Piplus_vs_EnuTrue->SetLineWidth(2);
events1Piplus_vs_EnuTrue->GetXaxis()->SetTitle("neutrino energy");
events1Piplus_vs_EnuTrue->GetYaxis()->SetTitle("1pi+ event distribution (flux folded xsec)");
events1Piplus_vs_EnuTrue->Draw();

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

energyReconstruction_1piplus->cd(2);

events1Pi0_vs_EnuTrue->SetLineColor(2); //red
events1Pi0_vs_EnuTrue->SetLineWidth(2);
events1Pi0_vs_EnuTrue->GetXaxis()->SetTitle("neutrino energy");
events1Pi0_vs_EnuTrue->GetYaxis()->SetTitle("1pi0 event distribution (flux folded xsec)");
events1Pi0_vs_EnuTrue->Draw();

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

reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

energyReconstruction_1piplus->Update();
energyReconstruction_1piplus->Print("1piplus_energyReconstruction.eps");












///////////////////////   event distribution versus neutrino energy for 0 pions and fixed number of nucleons

THStack * events0pions_Allp_vsEnu_origin = new THStack("events0pions_Allp_vsEnu","0pions Xprotons: event distr vs E_{#nu}");
for(int i=1; i<=numOrigins; i++){events0pions_Allp_vs_EnuTrue[i]->SetFillColor(i+2);
                                 events0pions_Allp_vs_EnuTrue[i]->SetMarkerStyle(21);
                                 events0pions_Allp_vs_EnuTrue[i]->SetMarkerColor(i+2);
                                 events0pions_Allp_vs_EnuTrue[i]->Write();
                                 events0pions_Allp_vsEnu_origin->Add(events0pions_Allp_vs_EnuTrue[i]); 
                       }


THStack * events0pions_0p_vsEnu_origin = new THStack("events0pions_Allp_vsEnu","0pions 0protons: event distr vs E_{#nu}");
for(int i=1; i<=numOrigins; i++){events0pions_vs_EnuTrue[0][i]->SetFillColor(i+2);
                                 events0pions_vs_EnuTrue[0][i]->SetMarkerStyle(21);
                                 events0pions_vs_EnuTrue[0][i]->SetMarkerColor(i+2);
                                 events0pions_vs_EnuTrue[0][i]->Write();
                                 events0pions_0p_vsEnu_origin->Add(events0pions_vs_EnuTrue[0][i]); 
                       }


THStack * events0pions_1p_vsEnu_origin = new THStack("events0pions_1p_vsEnu","0pions 1proton: event distr vs E_{#nu}");
for(int i=1; i<=numOrigins; i++){events0pions_vs_EnuTrue[1][i]->SetFillColor(i+2);
                                 events0pions_vs_EnuTrue[1][i]->SetMarkerStyle(21);
                                 events0pions_vs_EnuTrue[1][i]->SetMarkerColor(i+2);
                                 events0pions_vs_EnuTrue[1][i]->Write();
                                 events0pions_1p_vsEnu_origin->Add(events0pions_vs_EnuTrue[1][i]); 
                       }



THStack * events0pions_2p_vsEnu_origin = new THStack("events0pions_2p_vsEnu","0pions 2protons: event distr vs E_{#nu}");
for(int i=1; i<=numOrigins; i++){events0pions_vs_EnuTrue[2][i]->SetFillColor(i+2);
                                 events0pions_vs_EnuTrue[2][i]->SetMarkerStyle(21);
                                 events0pions_vs_EnuTrue[2][i]->SetMarkerColor(i+2);
                                 events0pions_vs_EnuTrue[2][i]->Write();
                                 events0pions_2p_vsEnu_origin->Add(events0pions_vs_EnuTrue[2][i]); 
                       }


///////////////////////   drawing 0-pion X,0,1,2-proton event disstributions versus Enu


TCanvas * eventDistr_origin_0pions = new TCanvas("0pions_eventDsitr_origin", "0#pi: Event distribution vs E_{#nu}");

eventDistr_origin_0pions->Divide(2,2);


eventDistr_origin_0pions->cd(1);

events0pions_Allp_vsEnu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);



eventDistr_origin_0pions->cd(2);
events0pions_0p_vsEnu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


eventDistr_origin_0pions->cd(3);
events0pions_1p_vsEnu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


eventDistr_origin_0pions->cd(4);
events0pions_2p_vsEnu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


eventDistr_origin_0pions->Update();
eventDistr_origin_0pions->Print("0pions_eventDistr_Enu.eps");
eventDistr_origin_0pions->Write();





//  drawing the same distributions with ArgoNeut barnumu barnu-mode data points 
TCanvas * ArgoNeuT_eventDistr_origin_0pions = new TCanvas("0pions_ArgoNeuT_eventDsitr_origin", "0#pi: ArgoNeuT event distribution vs E_{#nu}");

ArgoNeuT_eventDistr_origin_0pions->Divide(2,2);


ArgoNeuT_eventDistr_origin_0pions->cd(1);
TGraphErrors * ArgoNeut_0pi_Xp_Enu_barnumu_barnumode = new TGraphErrors("events/ArgoNeuT-barnumu-barnumode-0pions-Xp-Enu.exp","%lg %lg %lg %lg");
// divide data points be A=40
Int_t A=40;
Float_t scale=1./35.; //scale from events to xsec
for (int i=0;i<ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->GetN();i++){ ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->GetY()[i] *= scale/A; 
                                                                   ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->GetEY()[i] *= scale/A; }

events0pions_Allp_vsEnu_origin->Draw();
ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->SetTitle("0pions: ArgoNeuT events with any number of protons");
ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->SetMarkerStyle(kFullCircle);
ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->SetMarkerColor(kBlack);
ArgoNeut_0pi_Xp_Enu_barnumu_barnumode->Draw("PSame");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);



ArgoNeuT_eventDistr_origin_0pions->cd(2);
TGraphErrors * ArgoNeut_0pi_0p_Enu_barnumu_barnumode = new TGraphErrors("events/ArgoNeuT-barnumu-barnumode-0pions-0p-Enu.exp","%lg %lg %lg %lg");
for (int i=0;i<ArgoNeut_0pi_0p_Enu_barnumu_barnumode->GetN();i++){ ArgoNeut_0pi_0p_Enu_barnumu_barnumode->GetY()[i] *= scale/A; 
                                                                   ArgoNeut_0pi_0p_Enu_barnumu_barnumode->GetEY()[i] *= scale/A; }
events0pions_0p_vsEnu_origin->Draw();
ArgoNeut_0pi_0p_Enu_barnumu_barnumode->SetTitle("0pions: ArgoNeuT events with 0 protons");
ArgoNeut_0pi_0p_Enu_barnumu_barnumode->SetMarkerStyle(kFullCircle);
ArgoNeut_0pi_0p_Enu_barnumu_barnumode->SetMarkerColor(kBlack);
ArgoNeut_0pi_0p_Enu_barnumu_barnumode->Draw("PSame");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_eventDistr_origin_0pions->cd(3);
TGraphErrors * ArgoNeut_0pi_1p_Enu_barnumu_barnumode = new TGraphErrors("events/ArgoNeuT-barnumu-barnumode-0pions-1p-Enu.exp","%lg %lg %lg %lg");
for (int i=0;i<ArgoNeut_0pi_1p_Enu_barnumu_barnumode->GetN();i++){ ArgoNeut_0pi_1p_Enu_barnumu_barnumode->GetY()[i] *= scale/A; 
                                                                   ArgoNeut_0pi_1p_Enu_barnumu_barnumode->GetEY()[i] *= scale/A; }
events0pions_1p_vsEnu_origin->Draw();
ArgoNeut_0pi_1p_Enu_barnumu_barnumode->SetTitle("0pions: ArgoNeuT events with 1 proton");
ArgoNeut_0pi_1p_Enu_barnumu_barnumode->SetMarkerStyle(kFullCircle);
ArgoNeut_0pi_1p_Enu_barnumu_barnumode->SetMarkerColor(kBlack);
ArgoNeut_0pi_1p_Enu_barnumu_barnumode->Draw("PSame");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_eventDistr_origin_0pions->cd(4);
events0pions_2p_vsEnu_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_eventDistr_origin_0pions->Update();
ArgoNeuT_eventDistr_origin_0pions->Print("0pions_ArgoNeuT_eventDistr_Enu.eps");
ArgoNeuT_eventDistr_origin_0pions->Write();







fileHists->Close();

}   // end Analyze::LoopEnergyReconstruction()









