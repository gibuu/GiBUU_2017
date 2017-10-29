// to change the threshold for nucleon detection change thresEkinN


# include "Particle.h"
using std::ofstream;

void Analyze::LoopNumberOutgoing(TString reaction)
{

const Float_t thresEkinN=0.021;		//0.21 for LAr experiments	
const Short_t numOrigins=6;


// ************************
// histograms versus number of particles
// ************************
Int_t numMaxPion=7;
TH1F *  eventsAll_vs_numberOutPions[4][numOrigins+1];  //4 is charge ALL,-1,0,1 // 7 is origin all, QE, Delta, highRES, 1pi-bgr, DIS, 2p2hQE
for (int ch=0; ch<=3; ch++){
for (int i=0; i<=numOrigins; i++){eventsAll_vs_numberOutPions[ch][i] = new TH1F ("eventsAll_vs_numberOutPions"+CToCharge[ch]+HistToName[i],    
									"Number of outgoing pions  "+CToCharge[ch]+".   Origin: "+HistToName[i], numMaxPion+1,-0.5,numMaxPion+0.5); }
}
Int_t numMaxNucleon=12;
TH1F * eventsAll_vs_numberOutNucleons[4][numOrigins+1];
for (int ch=0; ch<=3; ch++){
for (int i=0; i<=numOrigins; i++){eventsAll_vs_numberOutNucleons[ch][i]    = new TH1F ("eventsAll_vs_numberOutNucleons "+CToCharge[ch]+HistToName[i],    
									      "Number of outgoing nucleons "+CToCharge[ch]+".   Origin: "+HistToName[i], numMaxNucleon+1,-0.5,numMaxNucleon+0.5); }
}
Float_t invArea[4];
// as requested by LAr experiment
TH1F * events0pions_vs_numberOutNucleons[4][numOrigins+1];
for (int ch=0; ch<=3; ch++){
for (int i=0; i<=numOrigins; i++){events0pions_vs_numberOutNucleons[ch][i] = new TH1F ("events0pions_vs_numberOutNucleons "+CToCharge[ch]+HistToName[i],    
									      "0 pions: Number of outgoing nucleons "+CToCharge[ch]+".   Origin: "+HistToName[i], numMaxNucleon+1,-0.5,numMaxNucleon+0.5); }
}
// as requsted by MiniBooNE experiment 
TH1F * events1Piplus_vs_numberOutNucleons[4][numOrigins+1];
for (int ch=0; ch<=3; ch++){
for (int i=0; i<=numOrigins; i++){events1Piplus_vs_numberOutNucleons[ch][i] = new TH1F ("events1Piplus_vs_numberOutNucleons "+CToCharge[ch]+HistToName[i],    
									      "1#pi^{+}: Number of outgoing nucleons "+CToCharge[ch]+".   Origin: "+HistToName[i], numMaxNucleon+1,-0.5,numMaxNucleon+0.5); }
}


// **************************
// histograms versus kinetic energy of the outgoing nucleons (sum over nucleons)
// **************************
char name[30], title[60];
Float_t Ekin_min=0, Ekin_max=2.0;
Int_t kk_Ekin=100;


TH1F *  nucleonsEkinSum[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"nucleonsEkinSum_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"sum Ekin for events with %d "+NucleonToCharge[ch],i);
	    nucleonsEkinSum[i][ch][ori] = new TH1F(name, title,kk_Ekin,Ekin_min,Ekin_max);
	}
    }
}
Float_t EkinBinWidth = nucleonsEkinSum[0][0][0]->GetBinWidth(0);        // bins width to be used in determining xsec as perweight/bin vs Ekin



TH1F *  nucleonsEkinSum_0pions[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"0pions_nucleonsEkinSum_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"0pions: sum Ekin for events with %d "+NucleonToCharge[ch],i);
	    nucleonsEkinSum_0pions[i][ch][ori] = new TH1F(name, title,kk_Ekin,Ekin_min,Ekin_max);
	}
    }
}


TH1F *  nucleonsEkinSum_1Piplus[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"1Piplus_nucleonsEkinSum_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"1#pi^{+}: sum Ekin for events with %d "+NucleonToCharge[ch],i);
	    nucleonsEkinSum_1Piplus[i][ch][ori] = new TH1F(name, title,kk_Ekin,Ekin_min,Ekin_max);
	}
    }
}






// **************************
// histograms versus kinetic energy of the outgoing nucleons (each nucleon)
// **************************
TH1F *  nucleonsEkin[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"nucleonsEkin_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"Ekin for events with %d "+NucleonToCharge[ch],i);
	    nucleonsEkin[i][ch][ori] = new TH1F(name, title,kk_Ekin,Ekin_min,Ekin_max);
	}
    }
}


TH1F *  nucleonsEkin_0pions[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"0pions_nucleonsEkin_%d_"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"0pions: Ekin for events with %d "+NucleonToCharge[ch],i);
	    nucleonsEkin_0pions[i][ch][ori] = new TH1F(name, title,kk_Ekin,Ekin_min,Ekin_max);
	}
    }
}






// **************************
// histograms versus angle of the outgoing protons with respect to neutrino beam (which is z-direction in GiBUU)
// **************************
Float_t theta_min=0, theta_max=180;
Int_t kk_theta=36;

TH1F *  nucleonsTheta[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"nucleonsTheta_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title, NucleonToCharge[ch]+" angle for events with %d "+NucleonToCharge[ch],i);
	    nucleonsTheta[i][ch][ori] = new TH1F(name, title,kk_theta,theta_min,theta_max);
	}
    }
}
Float_t thetaBinWidth = nucleonsTheta[0][0][0]->GetBinWidth(0);        // bins width to be used in determining xsec as perweight/bin vs Ekin


TH1F *  nucleonsTheta_0pions[numMaxNucleon+1][4][numOrigins+1];
for (int i=0; i<=numMaxNucleon; i++){
    for (int ch=0; ch<=3; ch++){ 
	for (int ori=0; ori<=numOrigins; ori++){
	    sprintf(name,"0pions_nucleonsTheta_%d"+NucleonToCharge[ch]+HistToName[ori],i);
	    sprintf(title,"0pions: "+NucleonToCharge[ch]+" angle for events with %d "+NucleonToCharge[ch],i);
	    nucleonsTheta_0pions[i][ch][ori] = new TH1F(name, title,kk_theta,theta_min,theta_max);
	}
    }
}









// initializing number of particles for counting
std::vector<Int_t> numberOutPions(4);       //[0] is total number of pions,    [1] of pi-, [2] of pi0, [3] of pi+
std::vector<Int_t> numberOutNucleons(4);    //[0] is total number of nucleons, [1] of antiprotons, [2] of neutrons, [3] of protons
std::vector<Float_t> sumEkinN(4);

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


	// set number of each particle species to 0   and sum of kinetic energies to 0
	for(Int_t i=0; i<4; i++){numberOutPions[i]=0; numberOutNucleons[i]=0;  sumEkinN[i]=0;} 



std::vector<Float_t> angleRadian(numberOutParticles);    // polar angle of the outgoing particle

        // set angle and kinetic energy  of outgoing particles to zero
        for (int i=0; i<numberOutParticles; i++) {angleRadian[i]=0;}





//########################################################################  distributions BEFORE particle counting

                /*   
        	Float_t Q2, Wfree;
        	Float_t emu, pxmu, pymu, pzmu; 
                

	        // ////////////  particle-loop-2  for recording muon
		for (Int_t  j=0; j<numberOutParticles; j++){
                if(particleID[j]==903){ emu=energy[j]; pxmu=momentumX[j]; pymu=momentumY[j]; pzmu=momentumZ[j];  }              // muon energy and momentum
                } // end particle-loop-2


                
                MuonKinematics b(enu, emu,pxmu,pymu,pzmu);
            	       Q2=b.Q2();
            	       Wfree=b.Wfree();
            	       eventsAll_vs_Q2[0]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsAll_vs_Q2[KToHist[origin]]->Fill(Q2, perweight/Q2BinWidth);
            	       eventsAll_vs_Wfree[0]->Fill(Wfree, perweight/WBinWidth);
            	       eventsAll_vs_Wfree[KToHist[origin]]->Fill(Wfree, perweight/WBinWidth);
                */ 

//###########################################################################   counting particles in each event



//cout<<"    origin="<<origin<<"    enu= "<<enu<<"     perweight="<<perweight<<endl;
//cout<<"numberOutParticles= "<<numberOutParticles<<endl;

	// particle-loop-1  for counting particles in each event 
	for (Int_t  j=0; j<numberOutParticles; j++){
	                            //cout<<"\t\t\t j= "<<j<<"   particleId="<<particleID[j]<<endl;
	    // angle of each particle with respect of neutrino direction       
	    Particle part(energy[j],momentumX[j],momentumY[j],momentumZ[j]); 
	    angleRadian[j]=part.thetaRadian();

	    //pions
	    if(particleID[j]==101){ (numberOutPions[0])++;    (numberOutPions[(charge[j]+2)])++;     }
	    //nucleons // only those above threshold are counted
	    if( particleID[j]==1 && ((energy[j]-m_N) > thresEkinN) ){  numberOutNucleons[0]++;          numberOutNucleons[(charge[j]+2)]++;  
	                                                               sumEkinN[0] += (energy[j]-m_N);  sumEkinN[charge[j]+2] += (energy[j]-m_N);
	                                                            }
	    
	} // end particle-loop-1



        //check that the pion and nucleon number is consistent
        if(numberOutPions[0]!=(numberOutPions[1]+numberOutPions[2]+numberOutPions[3]))
                        { cout<<"entry="<<jentry<<"   numberOutPions= "<<numberOutPions[0]<<"    pi- "<<numberOutPions[1]<<"     pi0 "<<numberOutPions[2]<<"    pi+ "<<numberOutPions[3]<<endl; }
        if(numberOutNucleons[0]!=(numberOutNucleons[1]+numberOutNucleons[2]+numberOutNucleons[3]))
                        { cout<<"entry="<<jentry<<"   numberOutNucleons= "<<numberOutNucleons[0]<<"    pbar "<<numberOutNucleons[1]<<"     n "<<numberOutNucleons[2]<<"    p "<<numberOutNucleons[3]<<endl; }
        

//########################################################################  distributions AFTER particle counting

		// outgoing pion and nucleons for all events 
		for(int k=0; k<=3; k++){ //cout<<"origin is"<<origin<<"  K is "<<(KToHist[origin])<<endl;
            	        eventsAll_vs_numberOutPions[k][0]->Fill(numberOutPions[k], perweight);
            	        eventsAll_vs_numberOutPions[k][KToHist[origin]]->Fill(numberOutPions[k], perweight);
            	        eventsAll_vs_numberOutNucleons[k][0]->Fill(numberOutNucleons[k], perweight);
            	        eventsAll_vs_numberOutNucleons[k][KToHist[origin]]->Fill(numberOutNucleons[k], perweight);
				       }


		// outgoing nucleons for events with 0 pions
		if (numberOutPions[0]==0){
		for(int k=0; k<=3; k++){ 
            	        events0pions_vs_numberOutNucleons[k][0]->Fill(numberOutNucleons[k], perweight);
            	        events0pions_vs_numberOutNucleons[k][KToHist[origin]]->Fill(numberOutNucleons[k], perweight);
				       }
		}


		// outgoing nucleons for events with 1 pilus and no other pions
		if (numberOutPions[3]==1 && numberOutPions[0]==1){
		for(int k=0; k<=3; k++){ 
            	        events1Piplus_vs_numberOutNucleons[k][0]->Fill(numberOutNucleons[k], perweight);
            	        events1Piplus_vs_numberOutNucleons[k][KToHist[origin]]->Fill(numberOutNucleons[k], perweight);
				       }
		}



		// kinetic energy distrubtions (sum ober all nucleons)
//	        cout<<"Before Fill    numberoutNucleons "<<numberOutNucleons[0]<<"   Ekin "<<sumEkinN<<"   perweight "<<perweight<<endl;
		// of nucleons
		if (numberOutNucleons[0]>0 &&  numberOutNucleons[0]<numMaxNucleon){
		    // distributions for each nucleon charge for any number of outgoing nucleons 
		    for(int ch=0; ch<=3; ch++){
			nucleonsEkinSum[0][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			nucleonsEkinSum[0][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			
			if (numberOutPions[0]==0){ 
			    nucleonsEkinSum_0pions[0][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			    nucleonsEkinSum_0pions[0][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			}

			if (numberOutPions[3]==1){ 
			    nucleonsEkinSum_1Piplus[0][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			    nucleonsEkinSum_1Piplus[0][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			}

		    // distributions for each nucleon charge for fixed number of this-charge-outgoing-nucleons
		        if (numberOutNucleons[ch]>0){   nucleonsEkinSum[numberOutNucleons[ch]][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
							nucleonsEkinSum[numberOutNucleons[ch]][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
		        
		        
			    if (numberOutPions[0]==0){ 
				nucleonsEkinSum_0pions[numberOutNucleons[ch]][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
				nucleonsEkinSum_0pions[numberOutNucleons[ch]][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			    }

			    if (numberOutPions[3]==1){
				nucleonsEkinSum_1Piplus[numberOutNucleons[ch]][ch][0]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
				nucleonsEkinSum_1Piplus[numberOutNucleons[ch]][ch][KToHist[origin]]->Fill(sumEkinN[ch],perweight/EkinBinWidth);
			    }
		        
		        }
		    }// for int ch=0 loop
               }

	// particle-loop-2  filling histograms of angular and kinetic energy distributions of nucleons 
	if (numberOutNucleons[0]>0 &&  numberOutNucleons[0]<=numMaxNucleon){  
	    for (Int_t  j=0; j<numberOutParticles; j++){
	        
		if( particleID[j]==1 && ((energy[j]-m_N) > thresEkinN) ){ // if nucleon above threshold
			
	                // --------------------  for all events
	                //angles
			nucleonsTheta[0][0][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsTheta[0][0][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, each origin separately
			nucleonsTheta[0][charge[j]+2][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsTheta[0][charge[j]+2][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, each origin separately
			    if (numberOutNucleons[charge[j]+2]>0){
				nucleonsTheta[numberOutNucleons[charge[j]+2]][charge[j]+2][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);		      // fixed number of nucleons of fixed charge, any origin
				nucleonsTheta[numberOutNucleons[charge[j]+2]][charge[j]+2][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);    // fixed number of nucleons of fixed charge, each origin separately
			    }
			//kinetic energy
			nucleonsEkin[0][0][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsEkin[0][0][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, each origin separately
			nucleonsEkin[0][charge[j]+2][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsEkin[0][charge[j]+2][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, each origin separately
			    if (numberOutNucleons[charge[j]+2]>0){
				nucleonsEkin[numberOutNucleons[charge[j]+2]][charge[j]+2][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);		      // fixed number of nucleons of fixed charge, any origin
				nucleonsEkin[numberOutNucleons[charge[j]+2]][charge[j]+2][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);    // fixed number of nucleons of fixed charge, each origin separately
			    }
			// ----------------------for events with 0 pions 
			//angles
			if (numberOutPions[0]==0){ 
			nucleonsTheta_0pions[0][0][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsTheta_0pions[0][0][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, each origin separately
			nucleonsTheta_0pions[0][charge[j]+2][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsTheta_0pions[0][charge[j]+2][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);  // X outgoing nucleons, any charge, each origin separately
			    if (numberOutNucleons[charge[j]+2]>0){
				nucleonsTheta_0pions[numberOutNucleons[charge[j]+2]][charge[j]+2][0]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);		      // fixed number of nucleons of fixed charge, any origin
				nucleonsTheta_0pions[numberOutNucleons[charge[j]+2]][charge[j]+2][KToHist[origin]]->Fill(angleRadian[j]/TMath::Pi()*180.,perweight/thetaBinWidth);    // fixed number of nucleons of fixed charge, each origin separately
			    }
			//kinetic energy
			nucleonsEkin_0pions[0][0][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsEkin_0pions[0][0][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, each origin separately
			nucleonsEkin_0pions[0][charge[j]+2][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, any origin
			nucleonsEkin_0pions[0][charge[j]+2][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);  // X outgoing nucleons, any charge, each origin separately
			    if (numberOutNucleons[charge[j]+2]>0){
				nucleonsEkin_0pions[numberOutNucleons[charge[j]+2]][charge[j]+2][0]->Fill((energy[j]-m_N),perweight/EkinBinWidth);		      // fixed number of nucleons of fixed charge, any origin
				nucleonsEkin_0pions[numberOutNucleons[charge[j]+2]][charge[j]+2][KToHist[origin]]->Fill((energy[j]-m_N),perweight/EkinBinWidth);    // fixed number of nucleons of fixed charge, each origin separately
			    }
			}
			
		}
	    }  
	}// particle-loop-2  filling histograms of angular energy distributions   





   } // end loop-for-entries

//#######################
//####   LOOP OVER ENTRIES   IS OVER
//#######################



// saving histograms to a file
TFile * fileHists = new TFile("NumberOutgoing.root","RECREATE");







///////////////////  stacked   number of pions

THStack * eventsAll_vs_numberOutPions_origin = new THStack("numberOutPions","number of outgoing pions ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutPions[0][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutPions[0][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutPions[0][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutPions[0][i]->Write();
                        eventsAll_vs_numberOutPions_origin->Add(eventsAll_vs_numberOutPions[0][i]); 
                       }
// eventsAll_vs_numberOutPions_origin->GetXaxis()->SetNdivisions(numMaxPion,kFALSE);


THStack * eventsAll_vs_numberOutPionsplus_origin = new THStack("numberOutPionsplus","number of outgoing #pi^{+} ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutPions[3][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutPions[3][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutPions[3][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutPions[3][i]->Write();
                        eventsAll_vs_numberOutPionsplus_origin->Add(eventsAll_vs_numberOutPions[3][i]); 
                       }

THStack * eventsAll_vs_numberOutPions0_origin = new THStack("numberOutPion0","number of outgoing #pi^{0} ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutPions[2][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutPions[2][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutPions[2][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutPions[2][i]->Write();
                        eventsAll_vs_numberOutPions0_origin->Add(eventsAll_vs_numberOutPions[2][i]); 
                       }


THStack * eventsAll_vs_numberOutPionsminus_origin = new THStack("numberOutPionminus","number of outgoing #pi^{-} ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutPions[1][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutPions[1][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutPions[1][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutPions[1][i]->Write();
                        eventsAll_vs_numberOutPionsminus_origin->Add(eventsAll_vs_numberOutPions[1][i]); 
                       }






THStack * eventsAll_vs_numberOutPions_charge = new THStack("numberOutPions","number of outgoing pions: charge ");
for(int c=-1; c<=1; c++){ eventsAll_vs_numberOutPions[c+2][0]->SetFillColor(c+4);
                         eventsAll_vs_numberOutPions[c+2][0]->SetMarkerStyle(21);
                         eventsAll_vs_numberOutPions[c+2][0]->SetMarkerColor(c+4);
                         eventsAll_vs_numberOutPions[c+2][0]->Write();
                         eventsAll_vs_numberOutPions_charge->Add(eventsAll_vs_numberOutPions[c+2][0]); 
                       }





///////////////////  stacked   number of nucleons


THStack * eventsAll_vs_numberOutNucleons_origin = new THStack("numberOutNucleons","number of outgoing nucleons ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutNucleons[0][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutNucleons[0][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutNucleons[0][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutNucleons[0][i]->Write();
                        eventsAll_vs_numberOutNucleons_origin->Add(eventsAll_vs_numberOutNucleons[0][i]); 
                       }


THStack * eventsAll_vs_numberOutProtons_origin = new THStack("numberOutProtons","number of outgoing protons ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutNucleons[3][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutNucleons[3][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutNucleons[3][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutNucleons[3][i]->Write();
                        eventsAll_vs_numberOutProtons_origin->Add(eventsAll_vs_numberOutNucleons[3][i]); 
                       }


THStack * eventsAll_vs_numberOutNeutrons_origin = new THStack("numberOutNeutrons","number of outgoing Neutrons ");
for(int i=1; i<=numOrigins; i++){eventsAll_vs_numberOutNucleons[2][i]->SetFillColor(i+2);
                        eventsAll_vs_numberOutNucleons[2][i]->SetMarkerStyle(21);
                        eventsAll_vs_numberOutNucleons[2][i]->SetMarkerColor(i+2);
                        eventsAll_vs_numberOutNucleons[2][i]->Write();
                        eventsAll_vs_numberOutNeutrons_origin->Add(eventsAll_vs_numberOutNucleons[2][i]); 
                       }









///////////////////  stacked   number of nucleons for events with 0 pions



THStack * events0pions_vs_numberOutNucleons_origin = new THStack("numberOutNucleons","number of outgoing nucleons ");
for(int i=1; i<=numOrigins; i++){events0pions_vs_numberOutNucleons[0][i]->SetFillColor(i+2);
                        events0pions_vs_numberOutNucleons[0][i]->SetMarkerStyle(21);
                        events0pions_vs_numberOutNucleons[0][i]->SetMarkerColor(i+2);
                        events0pions_vs_numberOutNucleons[0][i]->Write();
                        events0pions_vs_numberOutNucleons_origin->Add(events0pions_vs_numberOutNucleons[0][i]); 
                       }


THStack * events0pions_vs_numberOutProtons_origin = new THStack("numberOutProtons","number of outgoing protons ");
for(int i=1; i<=numOrigins; i++){events0pions_vs_numberOutNucleons[3][i]->SetFillColor(i+2);
                        events0pions_vs_numberOutNucleons[3][i]->SetMarkerStyle(21);
                        events0pions_vs_numberOutNucleons[3][i]->SetMarkerColor(i+2);
                        events0pions_vs_numberOutNucleons[3][i]->Write();
                        events0pions_vs_numberOutProtons_origin->Add(events0pions_vs_numberOutNucleons[3][i]); 
                       }


// output percentage in data  files
ofstream outFile_NumerOutProtons_0pions("numberOutProtons_0pions.dat",ios::out);
if (!outFile_NumerOutProtons_0pions){ cerr<<"File  numberOutProtons_0pions.dat  could not be opened"<<endl;
				      exit(1);
				    }
outFile_NumerOutProtons_0pions<<"# 1:numberOutProtons	2:all 3:QE 4:Delta 5:hihgRES  6:1piBGR 7:DIS 8:2p2h-NN     for proton Ekin threshold >"<<thresEkinN<<endl;
for (int k=1; k<=(events0pions_vs_numberOutNucleons[3][0]->GetSize()-2); k++){
    outFile_NumerOutProtons_0pions<<events0pions_vs_numberOutNucleons[3][0]->GetBinCenter(k)<<"\t";
    for(int i=0; i<=numOrigins; i++){
        outFile_NumerOutProtons_0pions<<events0pions_vs_numberOutNucleons[3][i]->GetBinContent(k)<<"\t";
    }
    outFile_NumerOutProtons_0pions<<endl;
}
outFile_NumerOutProtons_0pions.close();



THStack * events0pions_vs_numberOutNeutrons_origin = new THStack("numberOutNeutrons","number of outgoing neutrons ");
for(int i=1; i<=numOrigins; i++){events0pions_vs_numberOutNucleons[2][i]->SetFillColor(i+2);
                        events0pions_vs_numberOutNucleons[2][i]->SetMarkerStyle(21);
                        events0pions_vs_numberOutNucleons[2][i]->SetMarkerColor(i+2);
                        events0pions_vs_numberOutNucleons[2][i]->Write();
                        events0pions_vs_numberOutNeutrons_origin->Add(events0pions_vs_numberOutNucleons[2][i]); 
                       }





///////////////////  stacked   number of nucleons for events with 1 pi+



THStack * events1Piplus_vs_numberOutNucleons_origin = new THStack("numberOutNucleons","number of outgoing nucleons ");
for(int i=1; i<=numOrigins; i++){events1Piplus_vs_numberOutNucleons[0][i]->SetFillColor(i+2);
                        events1Piplus_vs_numberOutNucleons[0][i]->SetMarkerStyle(21);
                        events1Piplus_vs_numberOutNucleons[0][i]->SetMarkerColor(i+2);
                        events1Piplus_vs_numberOutNucleons[0][i]->Write();
                        events1Piplus_vs_numberOutNucleons_origin->Add(events1Piplus_vs_numberOutNucleons[0][i]); 
                       }


THStack * events1Piplus_vs_numberOutProtons_origin = new THStack("numberOutProtons","number of outgoing protons ");
for(int i=1; i<=numOrigins; i++){events1Piplus_vs_numberOutNucleons[3][i]->SetFillColor(i+2);
                        events1Piplus_vs_numberOutNucleons[3][i]->SetMarkerStyle(21);
                        events1Piplus_vs_numberOutNucleons[3][i]->SetMarkerColor(i+2);
                        events1Piplus_vs_numberOutNucleons[3][i]->Write();
                        events1Piplus_vs_numberOutProtons_origin->Add(events1Piplus_vs_numberOutNucleons[3][i]); 
                       }


THStack * events1Piplus_vs_numberOutNeutrons_origin = new THStack("numberOutNeutrons","number of outgoing neutrons ");
for(int i=1; i<=numOrigins; i++){events1Piplus_vs_numberOutNucleons[2][i]->SetFillColor(i+2);
                        events1Piplus_vs_numberOutNucleons[2][i]->SetMarkerStyle(21);
                        events1Piplus_vs_numberOutNucleons[2][i]->SetMarkerColor(i+2);
                        events1Piplus_vs_numberOutNucleons[2][i]->Write();
                        events1Piplus_vs_numberOutNeutrons_origin->Add(events1Piplus_vs_numberOutNucleons[2][i]); 
                       }







////////////////////////////////////////////////   legends 

TLegend *legendOrigin = new TLegend(0.7, 0.4, 0.92, 0.7);
legendOrigin->SetBorderSize(0);
for(int i=1; i<=numOrigins; i++){ legendOrigin->AddEntry(eventsAll_vs_numberOutPions[1][i], HistToName[i], "F");   }


TLegend * legendCharge = new TLegend(0.6,0.6,0.85,0.8);
legendCharge->SetBorderSize(0);
for(int c=-1; c<=1; c++){ legendCharge->AddEntry(eventsAll_vs_numberOutPions[c+2][0], CToCharge[c+2], "F");   }


legendOrigin->Clear();
legendOrigin->SetBorderSize(0);
for(int i=1; i<=numOrigins; i++){ legendOrigin->AddEntry(eventsAll_vs_numberOutPions[0][i], HistToName[i], "F");   }





///////////////////////////////   reaction title and final state
TLatex * reac = new TLatex();
reac->SetTextAlign(33);
reac->SetNDC();
TLatex * fs = new TLatex();
fs->SetTextAlign(33);
fs->SetNDC();
TString fs_0pions = "0#pi";
TString fs_1piplus = "1#pi^{+} 0#pi^{0} 0#pi^{-} ";





///////////////////////////////////////   plotting stacked histograms
/////////////// pions 


TCanvas * numberOutPions_origin = new TCanvas("numberOutPions_origin", "number outgoing pions: various origins");
numberOutPions_origin->Divide(3,2);

numberOutPions_origin->cd(1);
eventsAll_vs_numberOutPions_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);

numberOutPions_origin->cd(2);
eventsAll_vs_numberOutPionsplus_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);

numberOutPions_origin->cd(3);
eventsAll_vs_numberOutPions0_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);

numberOutPions_origin->cd(4);
eventsAll_vs_numberOutPionsminus_origin->Draw();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);


numberOutPions_origin->cd(5);
eventsAll_vs_numberOutPions_charge->Draw();
legendCharge->Draw();
reac->DrawLatex(0.91,0.9,reaction);

numberOutPions_origin->Update();
numberOutPions_origin->Print("numberOutPions_origin.eps");




/////////////// nucleons

TCanvas * numberOutNucleons_origin = new TCanvas("numberOutNucleons_origin", "number outgoing nucleons");
numberOutNucleons_origin->Divide(3,2);

numberOutNucleons_origin->cd(1);
eventsAll_vs_numberOutNucleons_origin->Draw();
eventsAll_vs_numberOutNucleons_origin->GetXaxis()->SetTitle("number of outgoing nucleons");
eventsAll_vs_numberOutNucleons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);

numberOutNucleons_origin->cd(2);
eventsAll_vs_numberOutProtons_origin->Draw();
eventsAll_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
eventsAll_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
eventsAll_vs_numberOutProtons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);


numberOutNucleons_origin->cd(3);
eventsAll_vs_numberOutNeutrons_origin->Draw();
eventsAll_vs_numberOutNeutrons_origin->GetXaxis()->SetLimits(-0.5,7.5);
eventsAll_vs_numberOutNeutrons_origin->GetXaxis()->SetTitle("number of outgoing neutrons");
eventsAll_vs_numberOutNeutrons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
numberOutNucleons_origin->Update();


numberOutNucleons_origin->cd(5);

for (int ch=0; ch<=3; ch++){ invArea[ch]=1./eventsAll_vs_numberOutNucleons[ch][0]->Integral();}


THStack * percentageAll_vs_numberOutProtons_origin = new THStack("percentageOutProtons","percentage of outgoing protons ");
for(int i=1; i<=numOrigins; i++){ eventsAll_vs_numberOutNucleons[3][i]->Scale(invArea[3]);
                         percentageAll_vs_numberOutProtons_origin->Add(eventsAll_vs_numberOutNucleons[3][i]); }

percentageAll_vs_numberOutProtons_origin->Draw();
percentageAll_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
percentageAll_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
percentageAll_vs_numberOutProtons_origin->GetYaxis()->SetTitle("percentage of total");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
numberOutNucleons_origin->Update();
numberOutNucleons_origin->Print("numberOutNucleons_origin.eps");






/////////////// nucleons  with 0 pions 

TCanvas * numberOutNucleons_0pions_origin = new TCanvas("0pions_numberOutNucleons_origin", "0 pions: number outgoing nucleons");

numberOutNucleons_0pions_origin->Divide(3,2);


numberOutNucleons_0pions_origin->cd(1);
events0pions_vs_numberOutNucleons_origin->Draw();
events0pions_vs_numberOutNucleons_origin->GetXaxis()->SetTitle("number of outgoing nucleons");
events0pions_vs_numberOutNucleons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

numberOutNucleons_0pions_origin->cd(2);
events0pions_vs_numberOutProtons_origin->Draw();
events0pions_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
events0pions_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
events0pions_vs_numberOutProtons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

numberOutNucleons_0pions_origin->cd(3);
events0pions_vs_numberOutNeutrons_origin->Draw();
events0pions_vs_numberOutNeutrons_origin->GetXaxis()->SetLimits(-0.5,7.5);
events0pions_vs_numberOutNeutrons_origin->GetXaxis()->SetTitle("number of outgoing neutrons");
events0pions_vs_numberOutNeutrons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


numberOutNucleons_0pions_origin->cd(5);
for (int ch=0; ch<=3; ch++){ invArea[ch]=1./events0pions_vs_numberOutNucleons[ch][0]->Integral();}
THStack * percentage0pions_vs_numberOutProtons_origin = new THStack("0pion_percentageOutProtons","0 pions: percentage of outgoing protons ");
for(int i=1; i<=numOrigins; i++){ events0pions_vs_numberOutNucleons[3][i]->Scale(invArea[3]);
                         percentage0pions_vs_numberOutProtons_origin->Add(events0pions_vs_numberOutNucleons[3][i]); }
percentage0pions_vs_numberOutProtons_origin->Draw();
percentage0pions_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
percentage0pions_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
percentage0pions_vs_numberOutProtons_origin->GetYaxis()->SetTitle("percentage of total");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
numberOutNucleons_0pions_origin->Update();
numberOutNucleons_0pions_origin->Print("0pions_numberOutNucleons_origin.eps");





/////////////// the same  nucleons  with 0 pions  with ArgoNeuT data

TCanvas * ArgoNeuT_numberOutNucleons_0pions_origin = new TCanvas("0pions_ArgoNeuT_numberOutNucleons_origin", "0 pions: number outgoing nucleons");

ArgoNeuT_numberOutNucleons_0pions_origin->Divide(3,2);


ArgoNeuT_numberOutNucleons_0pions_origin->cd(1);
events0pions_vs_numberOutNucleons_origin->Draw();
events0pions_vs_numberOutNucleons_origin->GetXaxis()->SetTitle("number of outgoing nucleons");
events0pions_vs_numberOutNucleons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

ArgoNeuT_numberOutNucleons_0pions_origin->cd(2);
events0pions_vs_numberOutProtons_origin->Draw();
events0pions_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
events0pions_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
events0pions_vs_numberOutProtons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");

// divide data points by A=40
Int_t A=40;
Float_t scale=1./35.; // to covert events to the xsec
////////// for numu flux in barnumode
/*
TGraphErrors * ArgoNeut_protMult = new TGraphErrors("events/ArgoNeuT-numu-barnumode-0pions-proton-mult.exp","%lg %lg %lg %lg");
Afor (int i=0;i<ArgoNeut_protMult->GetN();i++){ArgoNeut_protMult->GetX()[i] -=0.5;  ArgoNeut_protMult->GetY()[i] *= scale/A; ArgoNeut_protMult->GetEY()[i] *= scale/A; }
*/
//// for barnumu flux in barnumode
/*
TGraphErrors * ArgoNeut_protMult = new TGraphErrors("events/ArgoNeuT-barnumu-barnumode-0pions-proton-mult.exp","%lg %lg %lg %lg");
for (int i=0;i<ArgoNeut_protMult->GetN();i++){ArgoNeut_protMult->GetX()[i] -=0.5; ArgoNeut_protMult->GetY()[i] *= scale/A; ArgoNeut_protMult->GetEY()[i] *= scale/A; }

ArgoNeut_protMult->SetTitle("0 pions: ArgoNeut proton multiplicity");
ArgoNeut_protMult->SetMarkerStyle(kFullCircle);
ArgoNeut_protMult->SetMarkerColor(kBlack);
ArgoNeut_protMult->Draw("PSame");
*/
//////////
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);

ArgoNeuT_numberOutNucleons_0pions_origin->cd(3);
events0pions_vs_numberOutNeutrons_origin->Draw();
events0pions_vs_numberOutNeutrons_origin->GetXaxis()->SetLimits(-0.5,7.5);
events0pions_vs_numberOutNeutrons_origin->GetXaxis()->SetTitle("number of outgoing neutrons");
events0pions_vs_numberOutNeutrons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);


ArgoNeuT_numberOutNucleons_0pions_origin->cd(5);
percentage0pions_vs_numberOutProtons_origin->Draw();
percentage0pions_vs_numberOutProtons_origin->GetXaxis()->SetLimits(-0.5,7.5);
percentage0pions_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
percentage0pions_vs_numberOutProtons_origin->GetYaxis()->SetTitle("percentage of total");

//// for numu flux in numode

TGraphErrors * ArgoNeut_protMult = new TGraphErrors("events/ArgoNeuT-numu-numode-0pions-proton-mult.exp","%lg %lg %lg %lg");
for (int i=0;i<ArgoNeut_protMult->GetN();i++){ArgoNeut_protMult->GetX()[i] -=0.5;  ArgoNeut_protMult->GetY()[i] *= 1.; ArgoNeut_protMult->GetEY()[i] *= 1.; }
ArgoNeut_protMult->SetTitle("0 pions: ArgoNeut proton multiplicity");
ArgoNeut_protMult->SetMarkerStyle(kFullCircle);
ArgoNeut_protMult->SetMarkerColor(kBlack);
ArgoNeut_protMult->Draw("PSame");

////

legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
ArgoNeuT_numberOutNucleons_0pions_origin->Update();
ArgoNeuT_numberOutNucleons_0pions_origin->Print("0pions_ArgoNeuT_numberOutNucleons_origin.eps");










/////////////// nucleons  with 1 pi+


TCanvas * numberOutNucleons_1Piplus_origin = new TCanvas("1Piplus_numberOutNucleons_origin", "1 #pi^{+}: number outgoing nucleons");

numberOutNucleons_1Piplus_origin->Divide(3,2);

//for (int i=1; i<=numOrigins; i++){ numberOutNucleons_1Piplus_origin->cd(i);
//numberOutNucleons_1Piplus_origin->SetLeftMargin(0.15);
//numberOutNucleons_1Piplus_origin->SetRightMargin(0.05);
//}


numberOutNucleons_1Piplus_origin->cd(1);
events1Piplus_vs_numberOutNucleons_origin->Draw();
events1Piplus_vs_numberOutNucleons_origin->GetXaxis()->SetTitle("number of outgoing nucleons");
events1Piplus_vs_numberOutNucleons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
events1Piplus_vs_numberOutNucleons_origin->GetXaxis()->CenterTitle();
events1Piplus_vs_numberOutNucleons_origin->GetYaxis()->CenterTitle();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_1piplus);
numberOutNucleons_1Piplus_origin->Update();


numberOutNucleons_1Piplus_origin->cd(2);
events1Piplus_vs_numberOutProtons_origin->Draw();
events1Piplus_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
events1Piplus_vs_numberOutProtons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
events1Piplus_vs_numberOutProtons_origin->GetXaxis()->CenterTitle();
events1Piplus_vs_numberOutProtons_origin->GetYaxis()->CenterTitle();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_1piplus);
numberOutNucleons_1Piplus_origin->Update();


numberOutNucleons_1Piplus_origin->cd(3);
events1Piplus_vs_numberOutNeutrons_origin->Draw();
events1Piplus_vs_numberOutNeutrons_origin->GetXaxis()->SetTitle("number of outgoing neutrons");
events1Piplus_vs_numberOutNeutrons_origin->GetYaxis()->SetTitle("#sigma, 10^{-38} cm^{2}");
events1Piplus_vs_numberOutNeutrons_origin->GetXaxis()->CenterTitle();
events1Piplus_vs_numberOutNeutrons_origin->GetYaxis()->CenterTitle();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_1piplus);
numberOutNucleons_1Piplus_origin->Update();



numberOutNucleons_1Piplus_origin->cd(5);
for (int ch=0; ch<=3; ch++){ invArea[ch]=1./events1Piplus_vs_numberOutNucleons[ch][0]->Integral();}
THStack * percentage1Piplus_vs_numberOutProtons_origin = new THStack("1Piplus_percentageOutProtons","1 #pi^{+}: percentage of outgoing protons");
for(int i=1; i<=numOrigins; i++){ events1Piplus_vs_numberOutNucleons[3][i]->Scale(invArea[3]);
                         percentage1Piplus_vs_numberOutProtons_origin->Add(events1Piplus_vs_numberOutNucleons[3][i]); }
percentage1Piplus_vs_numberOutProtons_origin->Draw();
percentage1Piplus_vs_numberOutProtons_origin->GetXaxis()->SetTitle("number of outgoing protons");
percentage1Piplus_vs_numberOutProtons_origin->GetYaxis()->SetTitle("percentage of total");
percentage1Piplus_vs_numberOutProtons_origin->GetXaxis()->CenterTitle();
percentage1Piplus_vs_numberOutProtons_origin->GetYaxis()->CenterTitle();
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_1piplus);
numberOutNucleons_1Piplus_origin->Update();
numberOutNucleons_1Piplus_origin->Print("1Piplus_numberOutNucleons_origin.eps");















////////////////////////////////  kinetic energy distribution of nucleons (sum over nucleons)


TCanvas * sumKineticEnergy_outNucleons = new TCanvas("sumKineticEnergy_outNucleons", "Ekin  outgoing nucleons (sum)");

sumKineticEnergy_outNucleons->Divide(3,2);

for (int i=1; i<=6; i++){
nucleonsEkinSum[i][0][0]->SetLineColor(2);
nucleonsEkinSum[i][0][0]->SetLineWidth(2);
//nucleonsEkinSum[i]->SetMarkerStyle(21);
//nucleonsEkinSum[i]->SetMarkerColor(2);

sumKineticEnergy_outNucleons->cd(i);

nucleonsEkinSum[i][0][0]->Draw();
nucleonsEkinSum[i][0][0]->GetXaxis()->SetTitle("T_{N}, GeV");
nucleonsEkinSum[i][0][0]->GetXaxis()->CenterTitle();
nucleonsEkinSum[i][0][0]->GetYaxis()->SetTitle("d#sigma/dT_{N}, 10^{-38} cm^{2}/GeV");
nucleonsEkinSum[i][0][0]->GetYaxis()->CenterTitle();
//sprintf(title,"d#sigma/dT_{N}, 10^{-38} cm^{2}/GeV  for  %d nucleons",i);
//nucleonsEkinSum[i]->GetYaxis()->SetTitle(name);
reac->DrawLatex(0.91,0.9,reaction);
sumKineticEnergy_outNucleons->Update();
}

sumKineticEnergy_outNucleons->Print("sumKineticEnergy_outNucleons.eps");








TCanvas * sumKineticEnergy_outProtons = new TCanvas("sumKineticEnergy_outProtons", "Ekin outgoing protons (sum)");

sumKineticEnergy_outProtons->Divide(3,2);

for (int i=1; i<=6; i++){
nucleonsEkinSum[i][3][0]->SetLineColor(2);
nucleonsEkinSum[i][3][0]->SetLineWidth(2);

sumKineticEnergy_outProtons->cd(i);

nucleonsEkinSum[i][3][0]->Draw();
nucleonsEkinSum[i][3][0]->GetXaxis()->SetTitle("T_{p}, GeV");
nucleonsEkinSum[i][3][0]->GetXaxis()->CenterTitle();
nucleonsEkinSum[i][3][0]->GetYaxis()->SetTitle("d#sigma/dT_{p}, 10^{-38} cm^{2}/GeV");
nucleonsEkinSum[i][3][0]->GetYaxis()->CenterTitle();
reac->DrawLatex(0.91,0.9,reaction);
sumKineticEnergy_outProtons->Update();
}

sumKineticEnergy_outProtons->Print("sumKineticEnergy_outProtons.eps");






TCanvas * sumKineticEnergy_outProtons_0pions = new TCanvas("0pions_sumKineticEnergy_outProtons", "0pions: Ekin outgoing protons (sum)");

sumKineticEnergy_outProtons_0pions->Divide(3,2);

for (int i=1; i<7; i++){
nucleonsEkinSum_0pions[i][3][0]->SetLineColor(2);
nucleonsEkinSum_0pions[i][3][0]->SetLineWidth(2);

sumKineticEnergy_outProtons_0pions->cd(i);

nucleonsEkinSum_0pions[i][3][0]->Draw();
nucleonsEkinSum_0pions[i][3][0]->GetXaxis()->SetTitle("T_{p}, GeV");
nucleonsEkinSum_0pions[i][3][0]->GetXaxis()->CenterTitle();
nucleonsEkinSum_0pions[i][3][0]->GetYaxis()->SetTitle("d#sigma/dT_{p}, 10^{-38} cm^{2}/GeV");
nucleonsEkinSum_0pions[i][3][0]->GetYaxis()->CenterTitle();
nucleonsEkinSum_0pions[i][3][0]->GetYaxis()->SetTitleOffset(1.3);
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
sumKineticEnergy_outProtons_0pions->Update();
}

sumKineticEnergy_outProtons_0pions->Print("0pions_sumKineticEnergy_outProtons.eps");









TCanvas * sumKineticEnergy_outProtons_1Piplus = new TCanvas("1Piplus_sumKineticEnergy_outProtons", "1 #pi^{+}: Ekin outgoing protons (sum)");

sumKineticEnergy_outProtons_1Piplus->Divide(2,2);

for (int i=1; i<5; i++){
nucleonsEkinSum_1Piplus[i][3][0]->SetLineColor(2);
nucleonsEkinSum_1Piplus[i][3][0]->SetLineWidth(2);

sumKineticEnergy_outProtons_1Piplus->cd(i);

nucleonsEkinSum_1Piplus[i][3][0]->Draw();
nucleonsEkinSum_1Piplus[i][3][0]->GetXaxis()->SetTitle("T_{p}, GeV");
nucleonsEkinSum_1Piplus[i][3][0]->GetXaxis()->CenterTitle();
nucleonsEkinSum_1Piplus[i][3][0]->GetYaxis()->SetTitle("d#sigma/dT_{p}, 10^{-38} cm^{2}/GeV");
nucleonsEkinSum_1Piplus[i][3][0]->GetYaxis()->CenterTitle();
nucleonsEkinSum_1Piplus[i][3][0]->GetYaxis()->SetTitleOffset(1.3);
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_1piplus);
sumKineticEnergy_outProtons_1Piplus->Update();
}

sumKineticEnergy_outProtons_1Piplus->Print("1Piplus_sumKineticEnergy_outProtons.eps");











////////////////////////////////////////////
////////////////////////////////  kinetic energy distribution of protons (each proton)
////////////////////////////////////////////

TCanvas * kineticEnergy_outProtons = new TCanvas("kineticEnergy_outProtons", "Ekin outgoing protons");
//
kineticEnergy_outProtons->Divide(3,2);
//
for (int i=1; i<=6; i++){
nucleonsEkin[i][3][0]->SetLineColor(2);
nucleonsEkin[i][3][0]->SetLineWidth(3);
nucleonsEkin[i][3][0]->GetXaxis()->SetLimits(0.,0.9-i*0.1);
//
kineticEnergy_outProtons->cd(i);
(kineticEnergy_outProtons->GetPad(i))->SetLogy();
//
nucleonsEkin[i][3][0]->Draw("");
nucleonsEkin[i][3][0]->GetXaxis()->SetTitle("T_{p}, GeV");
nucleonsEkin[i][3][0]->GetXaxis()->CenterTitle();
nucleonsEkin[i][3][0]->GetYaxis()->SetTitle("d#sigma/dT_{p}, 10^{-38} cm^{2}/GeV");
nucleonsEkin[i][3][0]->GetYaxis()->CenterTitle();
    for (int ori=1; ori<=numOrigins; ori++){
    nucleonsEkin[i][3][ori]->SetLineColor(ori+2);
    nucleonsEkin[i][3][ori]->SetLineWidth(3);
    nucleonsEkin[i][3][ori]->GetXaxis()->SetLimits(0.,0.9-i*0.1);
    nucleonsEkin[i][3][ori]->Draw("Same");
    }
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
kineticEnergy_outProtons->Update();
}
//
kineticEnergy_outProtons->Print("KineticEnergy_outProtons.eps");




TCanvas * kineticEnergy_outProtons_0pions = new TCanvas("0pions_kineticEnergy_outProtons", "0pions: Ekin outgoing protons");
//
kineticEnergy_outProtons_0pions->Divide(3,2);
//
for (int i=1; i<=6; i++){
nucleonsEkin_0pions[i][3][0]->SetLineColor(2);
nucleonsEkin_0pions[i][3][0]->SetLineWidth(2);
nucleonsEkin_0pions[i][3][0]->GetXaxis()->SetLimits(0.,0.9-i*0.1);
//
kineticEnergy_outProtons_0pions->cd(i);
(kineticEnergy_outProtons_0pions->GetPad(i))->SetLogy();
//
nucleonsEkin_0pions[i][3][0]->Draw("");
nucleonsEkin_0pions[i][3][0]->GetXaxis()->SetTitle("T_{p}, GeV");
nucleonsEkin_0pions[i][3][0]->GetXaxis()->CenterTitle();
nucleonsEkin_0pions[i][3][0]->GetYaxis()->SetTitle("d#sigma/dT_{p}, 10^{-38} cm^{2}/GeV");
nucleonsEkin_0pions[i][3][0]->GetYaxis()->CenterTitle();
    for (int ori=1; ori<=numOrigins; ori++){
    nucleonsEkin_0pions[i][3][ori]->SetLineColor(ori+2);
    nucleonsEkin_0pions[i][3][ori]->SetLineWidth(3);
    nucleonsEkin_0pions[i][3][ori]->GetXaxis()->SetLimits(0.,0.9-i*0.1);
    nucleonsEkin_0pions[i][3][ori]->Draw("Same");
    }
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
legendOrigin->Draw();
kineticEnergy_outProtons_0pions->Update();
}
//
kineticEnergy_outProtons_0pions->Print("0pions_kineticEnergy_outProtons.eps");


















////////////////////////
/////////   angular distribution of protons 
///////////////////////



TCanvas * angle_outProtons = new TCanvas("angle_outProtons", "theta_{p} outgoing protons");
//
angle_outProtons->Divide(3,2);
//
for (int i=1; i<=6; i++){
nucleonsTheta[i][3][0]->SetLineColor(2);
nucleonsTheta[i][3][0]->SetLineWidth(3);
//
angle_outProtons->cd(i);
//
nucleonsTheta[i][3][0]->Scale(1000.);
nucleonsTheta[i][3][0]->Draw("");
nucleonsTheta[i][3][0]->GetXaxis()->SetTitle("#theta_{p}, degrees");
nucleonsTheta[i][3][0]->GetXaxis()->CenterTitle();
nucleonsTheta[i][3][0]->GetYaxis()->SetTitle("d#sigma/d#theta_{p}, 10^{-41} cm^{2}/GeV");
nucleonsTheta[i][3][0]->GetYaxis()->CenterTitle();
    for (int ori=1; ori<=numOrigins; ori++){
    nucleonsTheta[i][3][ori]->SetLineColor(ori+2);
    nucleonsTheta[i][3][ori]->SetLineWidth(3);
    nucleonsTheta[i][3][ori]->Scale(1000.);
    nucleonsTheta[i][3][ori]->Draw("Same");
    }
legendOrigin->Draw();
reac->DrawLatex(0.91,0.9,reaction);
angle_outProtons->Update();
}
//
angle_outProtons->Print("angle_outProtons.eps");






TCanvas * angle_outProtons_0pions = new TCanvas("0pions_angle_outProtons", "0pions: theta_{p} outgoing protons");
//
angle_outProtons_0pions->Divide(3,2);
//
for (int i=1; i<=6; i++){
nucleonsTheta_0pions[i][3][0]->SetLineColor(2);
nucleonsTheta_0pions[i][3][0]->SetLineWidth(3);
//
angle_outProtons_0pions->cd(i);
//
nucleonsTheta_0pions[i][3][0]->Scale(1000.);
nucleonsTheta_0pions[i][3][0]->Draw();
nucleonsTheta_0pions[i][3][0]->GetXaxis()->SetTitle("#theta_{p}, degrees");
nucleonsTheta_0pions[i][3][0]->GetXaxis()->CenterTitle();
nucleonsTheta_0pions[i][3][0]->GetYaxis()->SetTitle("d#sigma/d#theta_{p}, 10^{-41} cm^{2}/GeV");
nucleonsTheta_0pions[i][3][0]->GetYaxis()->CenterTitle();
    for (int ori=1; ori<=numOrigins; ori++){
    nucleonsTheta_0pions[i][3][ori]->SetLineColor(ori+2);
    nucleonsTheta_0pions[i][3][ori]->SetLineWidth(3);
    nucleonsTheta_0pions[i][3][ori]->Scale(1000.);
    nucleonsTheta_0pions[i][3][ori]->Draw("Same");
    }
reac->DrawLatex(0.91,0.9,reaction);
fs->DrawLatex(0.91,0.8,fs_0pions);
legendOrigin->Draw();
angle_outProtons_0pions->Update();
}
//
angle_outProtons_0pions->Print("0pions_angle_outProtons.eps");






fileHists->Close();

}   // end Analyze::LoopNumberOutgoing




