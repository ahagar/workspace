/***********************************************************
 * meri ugao iznedju rekonstruisanog i generisanog fotona, i ako je ugao manji ode 5, daje odnos energija
 * 	Read reconstructed photons from a  .slcio file, and draws histograms for different variables
 *
 *  Author:Goran Kačarević
 *
 *
 *
 ***********************************************************/



#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TVectorT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"
#include "Riostream.h"
#include "TGraph.h"


// LCIO includes
#include "lcio.h"
#include <IOIMPL/LCFactory.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "EVENT/ReconstructedParticle.h"
#include <UTIL/LCTOOLS.h>
#include <Exceptions.h>
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCRelation.h"
#include "IO/LCReader.h"
#endif

#include "stdlib.h"
#include "CandidateData.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cassert>


#include "varList.h"
#include "utils.h"

using namespace std;

const Double_t minPt1 = 15;			   //minimalna energija fotona
const Double_t minPt2 = 15;			   //minimalna energija fotona
const Double_t minPt3 = 18;			   //minimalna energija fotona
const Double_t minPt4 = 20;			   //minimalna energija fotona

const double mH = 126.0;               //Higgs mass
const double photonSize = 2;
const Double_t coneAngle = 2.5;		   //ugao konusa


Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif





Float_t electronEnergy=0;
Float_t positronEnergy=0;
Float_t thetaElectron=0 , thetaPositron=0;
Float_t noOfleptons =0, noOfElectrons = 0,noOfPositrons=0  ;


	TTree eventList("eventsSignal", "ILD event list");
	eventList.Branch("electronEnergy", &electronEnergy, "electronEnergy");
	eventList.Branch("positronEnergy", &positronEnergy, "positronEnergy");
	eventList.Branch("thetaElectron", &thetaElectron, "thetaElectron");
	eventList.Branch("thetaPositron", &thetaPositron, "thetaPositron");
	eventList.Branch("noOfleptons", &noOfleptons, "noOfleptons");
	eventList.Branch("noOfElectrons", &noOfElectrons, "noOfElectrons");
	eventList.Branch("noOfPositrons", &noOfPositrons, "noOfPositrons");



	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;




	//petlja koja iščitava .slcio podatke ukoliko ima više fajlova za jedan process
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{
		cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);
				try
				{
					lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));
				}
				catch(lcio::IOException &ex)
				{
					cout << ". Could not open.\n"; // Exception " << ex.what() << endl;
					continue;
				}
				cout << ". Reading.\n";


		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja
		cout << "Broj dogadjaja po fajlu je  : " << brojDogadjaja << endl;


		std::vector <TLorentzVector> allPhotons; //vektor koji prikuplja sve fotone
		Int_t brojac =0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
			brojac++;

			vector<TLorentzVector> recoElectrons;
			vector<TLorentzVector> mcElectrons;
			vector<TLorentzVector> mcPositrons;




				std::vector<std::string> colNames = *evt->getCollectionNames();
			//std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < (int)colNames.size(); ++i)
			{
			//std::cout << colNames[i] << endl;
			}



				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");
				bool signal = false;


				for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
				{
					IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

					if(mcParticle->getPDG() ==25 && abs(mcParticle->getDaughters()[0]->getPDG()) == 6 ) signal =true;

				}
				if (signal==false) continue;
			//	cout<< "signal je: " << signal<< endl;

				int noOfLeptons =0;
				int noOfElec =0;
				int noOfPos =0;


				for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
				{
					IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);


					if (mcParticle->getGeneratorStatus() != 1) continue;
					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

					const double *p = mcParticle->getMomentum(); // impuls čestice
					double e = mcParticle->getEnergy();	//energija čestice
					temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
					if(temp.E()<100) continue;

					Int_t particlePDG = mcParticle->getPDG();
					if(particlePDG == 11)
					{
							thetaElectron	= temp.Theta()  ;
							electronEnergy = temp.E();
							noOfLeptons++;
							noOfElec++;
							mcElectrons.push_back(temp);

					}
					if(particlePDG == -11)
					{
						thetaPositron = temp.Theta();
						positronEnergy = temp.E();
						noOfLeptons++;
						noOfPos++;
						mcPositrons.push_back(temp);

					}



				}
				noOfPositrons = noOfPos;
				noOfElectrons = noOfElec;
				noOfleptons = noOfLeptons;

				eventList.Fill();

				int noRecoElec =0, noRecoPos =0;

			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("PandoraPFANewPFOs");/*("PandoraPFANewPFOs");*/

			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici


				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
				if (temp.E()<100) continue;

				Int_t particlePDG = recParticle->getType();

				if(particlePDG == 11){

					noRecoElec++;
				}


				if(particlePDG == -11){
					noRecoPos++;

				}


			}





		//	cout << "broj finalnih elektrona u događaju: " << mcElectrons.size() << endl;


		} // End of event loop

		lcReader->close();

	} // End of file loop






	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	eventList.Write();





	rootFile.Close();



	return 0;
}



Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;
	UInt_t nFirstJob = 1;
	if(argc>iarg) nFirstJob = atoi(argv[iarg]); iarg++;
	UInt_t nLastJob = 10;
	if(argc>iarg) nLastJob   = atoi(argv[iarg]); iarg++;

	TString fName = "ee_rem_col_";
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "angleElectrons.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
