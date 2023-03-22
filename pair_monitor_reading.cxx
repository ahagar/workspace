/***********************************************************
 *
 * Read energy deposited from pairs in PairMonitor
 *
 *  Author:Goran Kačarević
 *  4.6.2015.
 *
 *
 ***********************************************************/

#include <EVENT/CalorimeterHit.h>
#include <RtypesCore.h>
#include <cmath>
#include <string>
#include <vector>



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
#include "EVENT/SimCalorimeterHit.h"
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


//function that calculates summ of energies from given collection
Double_t PairMonitorHitEnergies(EVENT::LCCollection* PairMonitorCollection)
{
	Double_t result = 0;

	for (Int_t i = 0; i < PairMonitorCollection->getNumberOfElements(); i++)
	{
		EVENT::CalorimeterHit* calHit = (EVENT::CalorimeterHit*)PairMonitorCollection->getElementAt(i);
		result += calHit->getEnergy();
		//double *pos = calHit->getPosition();
		// da li je pos NULL?
		//double x = pos[0];
	}

	return result;
}

double LenergyS1 = 0;
double LenergyS2 = 0;
double LenergyS3 = 0;
double LenergyS4 = 0;
double RenergyS1 = 0;
double RenergyS2 = 0;
double RenergyS3 = 0;
double RenergyS4 = 0;
int Lsektor1 = 0;
int Lsektor2 = 0;
int Lsektor3 = 0;
int Lsektor4 = 0;
int Rsektor1 = 0;
int Rsektor2 = 0;
int Rsektor3 = 0;
int Rsektor4 = 0;
int brojMC = 0;
double totalEnergy = 0; // ukupna Energija



Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif


	//gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/crtanje/histogrami/CLICdpStyle.C");
	//histogrami za različite kinematičke varijable
	TH1F histoTheta ("histoTheta", "Theta; #theta_{#gamma}", 90, 0, 180);
	TH1F histoPhi ("histoPhi", "Phi; #phi_{#gamma}", 50, 0, 3.14);
	TH1F energyOfPhotons ("energyOfPhotons", "PhotonEnergy; E_{#gamma} (GeV)", 3000, 0, 3000);
	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{#gamma#gamma} (GeV)", 3000, 0, 3000);
	TH1F histoCandidatePt ("CandidatePt", "CandidatePt; Pt_{#gamma#gamma} (GeV)", 3000, 0, 3000);
	TH1F histopT ("histopT", "pT; pT_{GeV}", 1000, 0, 1);




	TTree mcList("mcparticles", "mcparticles");
	varListPairMonitor vl; /* SL specific */
//	eventList.Branch("invM", &(vl.CandidateInvariantM), "invM");
	//eventList.Branch("CanE", &(vl.CandidateEnergy), "CanE");
	mcList.Branch("mcPt", &(vl.mcPt), "mcPt");
	mcList.Branch("mcE", &(vl.mcE), "mcE");
	mcList.Branch("mcTheta", &(vl.mcTheta), "mcTheta");

	TTree pairList("pairMonitor", "pairMonitor");
	//varListPairMonitor vl; /* SL specific */
	pairList.Branch("pairE", &(vl.pairE), "pairE");
	pairList.Branch("LenergyS1", &(vl.LenergyS1), "E (GeV)");
	pairList.Branch("LenergyS2", &(vl.LenergyS2), "E (GeV)");
	pairList.Branch("LenergyS3", &(vl.LenergyS3), "E (GeV)");
	pairList.Branch("LenergyS4", &(vl.LenergyS4), "E (GeV)");
	pairList.Branch("RenergyS1", &(vl.RenergyS1), "E (GeV)");
	pairList.Branch("RenergyS2", &(vl.RenergyS2), "E (GeV)");
	pairList.Branch("RenergyS3", &(vl.RenergyS3), "E (GeV)");
	pairList.Branch("RenergyS4", &(vl.RenergyS4), "E (GeV)");
	pairList.Branch("x_axis", &(vl.x_axis), "mm");
	pairList.Branch("y_axis", &(vl.y_axis), "mm");
	pairList.Branch("ro", &(vl.ro), "cm");
	pairList.Branch("azimuth", &(vl.azimuth), "rad");
	pairList.Branch("pdgPM", &(vl.pdgPM), "pdgPM");


	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;


	Double_t particlesTotal = 0;
	Double_t leftparticlesTotal = 0;
	Double_t rightparticlesTotal = 0;
	double eventEnergy = 0;



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

		int particleCounter = 0;
		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 )
		{

			cout <<"u eventu";
			std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
			IMPL::LCCollectionVec* pairParticles = (IMPL::LCCollectionVec*)evt->getCollection("PairMonitorCollection");
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");/*("PandoraPFANewPFOs");("PandoraPFOCollection");*/
			//IMPL::LCCollectionVec* beamCalParticles = (IMPL::LCCollectionVec*)evt->getCollection("BeamCalCollection");

			//	ene += PairMonitorHitEnergies(PairMonitorCollection);

			//petlja po MC česticama
			for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
			{

				brojMC++;
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

				TLorentzVector mcTemp;

				const double *p_mc = mcParticle->getMomentum(); // impuls čestice
				double e_mc = mcParticle->getEnergy();	//energija čestice
				mcTemp.SetPxPyPzE(p_mc[0], p_mc[1], p_mc[2], e_mc);  	//zapisujemo vrednosti energije i impulsa u četvorovektor


			/*	double mcPt = mcTemp.Pt();
				double mcE = mcTemp.E();
				double mcTheta = mcTemp.Theta();*/
				histopT.Fill(mcTemp.Pt());
				vl.mcPt = mcTemp.Pt();
				vl.mcE = mcTemp.E();
				vl.mcTheta = mcTemp.Theta();
				mcList.Fill();

			}//end MCParticles

			//Petlja po PairČesticama
			for (Int_t i = 0; i < pairParticles->getNumberOfElements() ; i++)
			{

				particleCounter++;
				EVENT::SimCalorimeterHit* pairParticle = (EVENT::SimCalorimeterHit*) pairParticles->getElementAt(i);







				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

			//	const double *p = pairParticle->; // impuls čestice
				double e = pairParticle->getEnergy();	//energija čestice
				//temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

		//		Int_t particlePDG = fabs(recParticle->getType());

				 particlesTotal += pairParticle->getNMCContributions();

				//cout <<"pairParticle->getNMCParticles() " << 	pairParticle->getPDGCont(i)<<endl;

				totalEnergy +=e;

				eventEnergy += e;

				vl.pairE += e;


				const float* coords = pairParticle->getPosition();

				double x = coords[0];
				double y = coords [1];
				double z = coords [2];

				double r = TMath::Sqrt(x*x+y*y);
				double phi = TMath::ATan2(y,x);

		//		cout <<"cellID: " << pairParticle->getCellID1() <<"; r : "<<r<<", phi: " <<phi<<endl;

			//	cout <<"x: "<< x<<", y: "<< y<< ", z: "<< z<< endl;


				// --------------------------------------
				int numContribs = pairParticle->getNMCContributions();

				// sve vektore ili "zbirne" promenljive koje ti trebaju deklarisi pre loopa
				std::vector<EVENT::MCParticle*> contribParticles;
				for (int k = 0; k < numContribs; k++)
				{
					EVENT::MCParticle* contrib = pairParticle->getParticleCont(k);

					// ako ti treba sama cestica, ubaci u neki vektor
					contribParticles.push_back(contrib);

					// ako ti treba neki property
					int contribPDG = contrib->getPDG();
					double contribE = contrib->getEnergy();
					const double* contribP = contrib->getMomentum();    // contribP je array




				// TODO:
				if(z > 0 && contribPDG == -11)
				{
				//	cout <<" pdg cestice deponovane u z > 0 je : "<< contribPDG<<endl;//
					rightparticlesTotal += 1;
					vl.ro = r;
					vl.azimuth = phi;
					vl.pdgPM = contribPDG;

// sad radi sta hoces :)

				// sektori
				if (x > 0 && y > 0 )
				{
					Rsektor1++;
					RenergyS1 +=e;
				}

				if (x > 0 && y < 0)
				{
					Rsektor2++;
					RenergyS2 +=e;

				}


				if (x < 0 && y < 0)
				{
					Rsektor3++;
					RenergyS3 +=e;

				}

				if (x < 0 && y > 0)
				{
					Rsektor4++;
					RenergyS4 +=e;

				}
				}//end z > 0


				//**************************************************
				if(z < 0 && contribPDG == 11)
				{
				leftparticlesTotal +=pairParticle->getNMCContributions();
				vl.x_axis = x;
				vl.y_axis = y;

				if (x > 0 && y > 0 )
				{
					Lsektor1++;
					LenergyS1 +=e;
				}

				if (x > 0 && y < 0)
				{
					Lsektor2++;
					LenergyS2 +=e;

				}


				if (x < 0 && y < 0)
				{
					Lsektor3++;
					LenergyS3 +=e;

				}

				if (x < 0 && y > 0)
				{
					Lsektor4++;
					LenergyS4 +=e;
				}
				pairList.Fill();

				}//end z < 0
				}//kraj contrib particles
		//	if (e>1) cout <<e<<endl;

			//cellID0()



			}  // end of particle loop



			cout <<"broj cestica: "<< particlesTotal<<endl;


			cout<<"za z > 0,   U-D = "<< (Rsektor1 + Rsektor2)- ( Rsektor3 + Rsektor4)<<", L-R = "<<(Rsektor1 - Rsektor4)-(Rsektor2+Rsektor3) <<endl;

		} // End of event loop

		lcReader->close();
		//cout<<"ukupno cestica je: "<<particleCounter<<endl;


	/*	vl.energyS1 = energyS1;
		vl.energyS2 = energyS2;
		vl.energyS3 = energyS3;
		vl.energyS4 = energyS4;*/


	} // End of file loop


	//cout <<"sektor1: "<<sektor1<<", sektor2: "<<sektor2<<" ,sektor1: "<<sektor3<< "sektor4: "<<sektor4<<endl;
	//cout <<"z > 0,   sektor1: "<<RenergyS1<<", sektor2: "<<RenergyS2<<" ,sektor3: "<<RenergyS3<< ", sektor4: "<<RenergyS4<<endl;

	double RupEnergy = RenergyS1+RenergyS4;
	double RdownEnergy = RenergyS2+RenergyS3;
	double RleftEnergy = RenergyS4 + RenergyS3;
	double RrightEnergy = RenergyS1+RenergyS2;
	double RsumEnergy= RenergyS1+RenergyS2+RenergyS3+RenergyS4;
	cout <<"broj monte karlo cestica: "<<brojMC<<endl;


	double LupEnergy = LenergyS1+LenergyS4;
	double LdownEnergy = LenergyS2+LenergyS3;
	double LleftEnergy = LenergyS4 + LenergyS3;
	double LrightEnergy = LenergyS1+LenergyS2;
	double LsumEnergy = LenergyS1 + LenergyS2 + LenergyS3 + LenergyS4 ;
	cout<<"prekoSektora energija : "<<RsumEnergy + LsumEnergy << ", i ukupna "<<totalEnergy<<endl;


//	cout<<"za z > 0,   U-D = "<< (Rsektor1 + Rsektor2)- ( Rsektor3 + Rsektor4)<<", L-R = "<<(Rsektor1 - Rsektor4)-(Rsektor2+Rsektor3) <<endl;
	cout<<"za z < 0,   U-D = "<< (LupEnergy-LdownEnergy)/LsumEnergy<<", L-R = "<<(LleftEnergy - LrightEnergy)/LsumEnergy<<endl;

	cout<<"energija po eventu je: "<<eventEnergy*1000<<" KeV"<<endl;

	cout<<"energija po čestici: "<<eventEnergy*1000/particlesTotal<<" KeV"<<endl;
	cout <<"pogoci u z < 0: "<< leftparticlesTotal<<", pogoci u z>0: "<< rightparticlesTotal<<endl;




	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	mcList.Write();
	pairList.Write();

	histopT.Write();

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

	TString rfName = "app.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
