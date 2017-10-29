//11.03.2012  
// analyzes the TTree from the GiBUUEvents.root file  (for NEUTRINO events)
//to execute this file run  root and  .X  GiBUUEventsToClass.cpp


void GiBUUEventsToClass(){

if (gSystem-> AccessPathName("GiBUUEvents.root")){
    gROOT->ProcessLine(".x convert-events-to-simpleTTree.cpp");
} //end if(gSystem ...)


// input info
TFile * f = new TFile("GiBUUEvents.root");
T->Print();
TTree * T = (TTree*)f->Get("T");

//output info
//TFile * f2 = new TFile("GiBUU-class-Events.root","RECREATE");
//TTree * T2 = new TTree("T2", "GiBUU events in the form of class with all kinematical variables");


Int_t		eventNumber;
Int_t		origin;
Float_t		perweight;
Int_t		numberOutParticles;
Int_t *			particleID;
Int_t *			charge;
Float_t *			positionX;
Float_t	*		positionY;
Float_t	*		positionZ;
Float_t	*		energy;
Float_t	*		momentumX;
Float_t	*		momentumY;
Float_t	*		momentumZ;
Long64_t *		history;



T->SetBranchAddress("origin", 			&origin);
T->SetBranchAddress("perweight", 		&perweight);
T->SetBranchAddress("numberOutParticles", 	&numberOutParticles);
T->SetBranchAddress("particleID", 		particleID);
T->SetBranchAddress("charge", 			charge);
T->SetBranchAddress("positionX", 		positionX);
T->SetBranchAddress("positionY", 		positionY);
T->SetBranchAddress("positionZ", 		positionZ);
T->SetBranchAddress("energy", 			energy);
T->SetBranchAddress("momentumX", 		momentumX);
T->SetBranchAddress("momentumY", 	 	momentumY);
T->SetBranchAddress("momentumZ", 		momentumZ);

int const kMaxNumberParticles=69;


//create two Histograms

TH1I * hNumberOutParticles = new TH1I("hNumberOutParticles", "number of outgoing particles", kMaxNumberParticles, 1, kMaxNumberParticles-1);
TH1I * hOrigin 		   = new TH1I("hOrigin", "origin of events", 40, 1, 40-1);
TH2I * hOutParticlesVsOrigin 	= new TH2I("hOutParticlesVsOrigin", "number of outgoing particles versus the origin of events", kMaxNumberParticles, 1, kMaxNumberParticles-1, 40, 1, 40-1);


//read all entries and fill the histograms
Int_t nEntries = (Int_t)T->GetEntries();
cout<<"number of entries="<<nEntries<<endl;
for (Int_t i=0; i<nEntries; i++){
    //T->Show(i);
    T->GetEntry(i);
    hNumberOutParticles->Fill(numberOutParticles);
    hOrigin->Fill(origin);
    hOutParticlesVsOrigin->Fill(origin,numberOutParticles);
} // end for-loop


// check if graphic mode is available
if (gROOT->IsBatch()) return;


// draw several histograms
TCanvas * c1 = new TCanvas("c1","initial info about GiBUU events");
c1->Divide(2,2);
// draw a histogram with a number of particles
c1->cd(1);
hNumberOutParticles->SetLineColor(2);
hNumberOutParticles->SetLineWidth(2);
hNumberOutParticles->GetXaxis()->SetTitle("number of outgoing particles");
hNumberOutParticles->GetYaxis()->SetTitle("number of events");
hNumberOutParticles->Draw();
// draw a histogram with an origin  of the events
c1->cd(2);
hOrigin->SetLineColor(3);
hOrigin->SetLineWidth(3);
hOrigin->GetXaxis()->SetTitle("origin of the events");
hOrigin->GetYaxis()->SetTitle("number of events");
hOrigin->Draw();
//
c1->cd(3);
hOutParticlesVsOrigin->SetOption("lego2");
hOutParticlesVsOrigin->GetXaxis()->SetTitle("origin of the events");
hOutParticlesVsOrigin->GetYaxis()->SetTitle("number of outgoing particles");
hOutParticlesVsOrigin->GetZaxis()->SetTitle("number of events");
hOutParticlesVsOrigin->Draw();


//open a browser and TreeViewer
//new TBrowser;
//T->StartViewer();


//t2.Branch("origin", 			&origin, 	"origin/I");
//t2.Branch("perweight", 			&preweight, 	"perweight/F");
//t2.Branch("numberOutParticles", 	&numberOutParticles, 	"numberOutParticles/I");


}
