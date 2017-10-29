//11.03.2012  
// convert the GiBUU output file "FinalEvents.dat" to the simple TTree which is written in the GiBUUEvents.root file
// to execute this file run  root and  .X  convert-events-to-simpleTTree.cpp

//20.03.2012
//changes to that is also reads neutrino energy


//26.11.2012
// DONE list
//  !   if the first character is # the line is recognized as comment
//  !   how to combine events for several nucleus, for example to CH2:
//		run H two times, C 12 times and combine files with linux "cat"
//		since C12 runs slower, lower numEnsembles can be used
//		(!!!) the num_runs_sameEnergy MUST be the same for all runs for all nuclei


//2.11.2012
// TODO list
// ?    how to add information  about the target nucleus, neutrino flavor, proceess(EM,CC,NC)  to the root file ?
// ?    how to add info about the initial neutrino flux
// ???: in GiBUU neutrino runs all electron, muon and tau have ID=903 and charge=0
// ???: info about the struck nucleon is not in the GiBUU event output, so we can not check energy conservation
// ???: in each run GiBUU provides events either after FSI (numTimeSteps>0) or before FSI (if numTimeSteps=0);
//         how to compare them in one plot ?
//author: Olga Lalakulich




{ 

gROOT->Reset();
gSystem.Load("libPhysics.so");


const Int_t kMaxNumberParticles=115;


// this class holds the variables for the events
struct TGiBUUevent{
	Int_t		eventNumber;
	Int_t		origin;
	Float_t		perweight;
	Float_t		enu;
	Int_t		numberOutParticles;
		Int_t 		particleID[kMaxNumberParticles];
		Int_t 		charge[kMaxNumberParticles];
		Float_t		positionX[kMaxNumberParticles];
		Float_t		positionY[kMaxNumberParticles];
		Float_t		positionZ[kMaxNumberParticles];
		Float_t		energy[kMaxNumberParticles];
		Float_t		momentumX[kMaxNumberParticles];
		Float_t		momentumY[kMaxNumberParticles];
		Float_t		momentumZ[kMaxNumberParticles];
		Long64_t	history[kMaxNumberParticles];
};

TGiBUUevent event;

// check if the file FinalEvents.dat is in the current directory
if (gSystem->AccessPathName("FinalEvents.dat")){
    cout<<"file FinalEvents.dat is not found in the current directory.\n"
	<<"Put it there and delete the header line. \n"
	<<"*****  The TTree and Root file are NOT created *****"<<endl;
    return;
} //end if(gSystem ...)


// open the ASCII file
FILE *fp = fopen("FinalEvents.dat", "r");
char line[205];

//create a new ROOT file
TFile *f = new TFile("GiBUUEvents.root", "RECREATE");

//create a TTree
TTree *tree = new TTree("T", "events from the ascii file");

//create branch for the event
tree->Branch("eventNumber", &event.eventNumber, "eventNumber/I");
tree->Branch("origin",      &event.origin,      "origin/I");
tree->Branch("perweight",   &event.perweight,   "perweight/F");
tree->Branch("enu",   &event.enu,         "enu/F");
tree->Branch("numberOutParticles", &event.numberOutParticles, "numberOutParticles/I");

tree->Branch("particleID", event.particleID, "particleID[numberOutParticles]/I");
tree->Branch("charge",     event.charge,     "charge[numberOutParticles]/I");
tree->Branch("positionX",     event.positionX,   "positionX[numberOutParticles]/F");
tree->Branch("positionY",     event.positionY,   "positionY[numberOutParticles]/F");
tree->Branch("positionZ",     event.positionZ,   "positionZ[numberOutParticles]/F");
tree->Branch("energy",        event.energy,      "energy[numberOutParticles]/F");
tree->Branch("momentumX",     event.momentumX,   "momentumX[numberOutParticles]/F");
tree->Branch("momentumY",     event.momentumY,   "momentumY[numberOutParticles]/F");
tree->Branch("momentumZ",     event.momentumZ,   "momentumZ[numberOutParticles]/F");
tree->Branch("history",       event.history,     "history[numberOutParticles]/L");









//
//fill the tree from the ASCII file
//
Int_t num=1, i=0, intDummy=1; //num=previous event number;  i= number of outparticles found
Float_t         floatDummy;
char character; // dummy to identify a new line


// find out what was the number_runs_sameEnergy  by reading the last line of the file
// this is because the perewight should be divided by number_runs_sameEnergy
// I googled for a better method to read the last line (like using linux "tail" inside C++) but found nothing reasonable
// Kai G suggested a good solution with fseek (see commented lines below), but for some reason that does not work with ROOT
Int_t num_runs_sameEnergy;
rewind(fp);
while (fgets (&line,205,fp)){
sscanf(&line[1],"%d", &num_runs_sameEnergy);
}
cout<<"num_runs_sameEnergy="<<num_runs_sameEnergy<<endl;
rewind(fp);



// gives Error: Symbol SEEK_END is not defined in current scope  convert-events-to-simpleTTree.cpp:128
// an if i try #include <cstdio>  does not work
// found in the internet:
// <cstdio>  The stdio library is made mostly obsolete by the newer iostream library, but many programs still use it. 
//There are facilities for random access files and greater control over output format, error handling, and temporary files. 
//  !!!  Mixing both I/O libraries is not recommended. There are no facilities for string I/O. 
/*
Int_t num_runs_sameEnergy;
rewind(fp);
fseek(fp, -172, SEEK_END);
fgets (&line,8,fp) 
sscanf(&line[1],"%d", &num_runs_sameEnergy);
cout<<"num_runs_sameEnergy="<<num_runs_sameEnergy<<endl;
*/


intDummy=0;
while (fgets (&line,205,fp)){
//check if the  line is a comment; if yes, skip the rest of the while-loop
sscanf(&line[0],"%c", &character);  // cout<<"character is"<<character<<endl;  cin>>intDummy;
if (character == '#'){continue;}
//read event number
Int_t currentEventNumber;
sscanf(&line[8],"%d", &currentEventNumber);  if (currentEventNumber%2000== 0){cout<<"n=  "<<intDummy<<"   currentEventNumber="<<currentEventNumber<<endl; intDummy++;}
//if event number is different from the previous one then create array for the outgoing particles
    if (currentEventNumber == num){ i++;} else {event.eventNumber=num;
                                                event.numberOutParticles=i;
                                                //cout<<"eventNumber="<<event.eventNumber<<"   i="<<i<<"  numberOutParticles="<<event.numberOutParticles; cin>>intDummy;
                                                //for(Int_t j=0; j<event.numberOutParticles; j++){cout<<event.particleID[j]<<"   "<<event.history[j]<<"   "<<event.positionX[j]<<endl;}
						tree->Fill();
						for(Int_t j=0; j<=kMaxNumberParticles-1; j++){event.particleID[j]=0;
									    event.charge[j]=0;
									    event.history[j]=0;
									    event.positionX[j]=0;
									    event.positionY[j]=0;
									    event.positionZ[j]=0;
									    event.energy[j]=0;
									    event.momentumX[j]=0;
									    event.momentumY[j]=0;
									    event.momentumZ[j]=0;}
                                                i=1; 
                                                num=currentEventNumber;
                                               }
    //read particle ID and charge
    sscanf(&line[18],"%d", &event.particleID[i-1]); //cout<<"particle "<<i<<"  particleID="<<event.particleID[i-1]<<endl;
    sscanf(&line[26],"%d", &event.charge[i-1]);
//read perweight of the event
sscanf(&line[29],"%f", &floatDummy);   event.perweight=floatDummy/num_runs_sameEnergy;
//    cout<<"perweight="<<event.perweight<<endl;
    //read particle position
    sscanf(&line[43],"%f %f %f", &event.positionX[i-1], &event.positionY[i-1], &event.positionZ[i-1]);
//    cout<<"eventNumber="<<currentEventNumber<<"    particle "<< i<<"    position= "<<event.positionX[i-1]<<"   "<<event.positionY[i-1]<<"   "<<event.positionZ[i-1]<<endl;
    //read particle 4-momentum
    sscanf(&line[85],"%f %f %f %f", &event.energy[i-1], &event.momentumX[i-1], &event.momentumY[i-1], &event.momentumZ[i-1]);
    //read the history of the particle
    sscanf(&line[141],"%d", &event.history[i-1]);
//read the origin of the event
sscanf(&line[156],"%d", &event.origin);
//read the neutrino energy
sscanf(&line[160],"%f", &event.enu);
} //end while 

//write the last event
event.eventNumber=num;
event.numberOutParticles=i;
tree->Fill();


//check what the tree looks like
tree->Print();
fclose(fp);
f->Write();

T->StartViewer();

} //end of the program




