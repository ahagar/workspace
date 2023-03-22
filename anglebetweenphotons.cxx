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
const Double_t maxConeEnergy = 20 ;	   //maksimalna energija konusa
const Double_t leptonPt = 20000;		   //energija leptona posle koje zadovoljavaju uslov leptonFound
const Double_t EnergyCenterMass = 3000; // energija u sistemu centra masa, koju koristim za missing energy
const Double_t minInvMass = 0;
const Double_t maxInvMass = 1400000;
const Double_t cutEremaining = 100;
const Double_t cutHiggsPt = 20;
const Double_t cutEhiggsMin = 100;
const Double_t cutEhiggsMax = 1000;
const Double_t cutCosinusHelicity = 0.9;
const Double_t minhiggsInvM = 110;
const Double_t maxhiggsInvM = 140;
const Double_t minTheta = 2;


Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif







	TTree eventList("eventsSignal", "ILD event list");
	varListGoran vl; /* SL specific */
	eventList.Branch("thetaMCe", &(vl.thetaMCe), "thetaMCe");
	eventList.Branch("thetaRecoe", &(vl.thetaRecoe), "thetaRecoe");
	eventList.Branch("angleE", &(vl.angleE), "angleE");
	eventList.Branch("photonsRatioEnergy", &(vl.photonsRatioEnergy), "photonsRatioEnergy");
	eventList.Branch("electronsRatioEnergy", &(vl.electronsRatioEnergy), "electronsRatioEnergy");

	eventList.Branch("photonsPerEvent", &(vl.NumberPhotonsbyEvent), "photonsPerEvent");
	eventList.Branch("energyPhotons", &(vl.energyOfAllPhotons), "energyPhotons");
	eventList.Branch("thetaPhotons", &(vl.thetaOfAllPhotons), "thetaPhotons");
	eventList.Branch("ptPhotons", &(vl.ptOfAllPhotons), "ptPhotons");
	eventList.Branch("highestPt", &(vl.HighestPhotonPt), "highestPt");
	eventList.Branch("2ndHighestPt", &(vl.secondHighestPhotonPt), "2ndHighestPt");







	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	int CounterHiggs =0;



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
			vector<TLorentzVector> photons;
			vector<TLorentzVector> HiggsPhotons;

			vector<TLorentzVector> recoElectrons;
			vector<TLorentzVector> mcElectrons;

			vector<Double_t> photonsECal;
			vector<Double_t> photonsHCal;

			vector<TLorentzVector> particles;   // other than photons
			vector <Double_t> PtPhotons;
			vector <Double_t> ConeEnergyOfPhotons;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			//std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < (int)colNames.size(); ++i)
			{
			//std::cout << colNames[i] << endl;
			}



				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");
				for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
				{
					IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

					if (mcParticle->getGeneratorStatus() != 1) continue;
					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

					const double *p = mcParticle->getMomentum(); // impuls čestice
					double e = mcParticle->getEnergy();	//energija čestice
					temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

					Int_t particlePDG = fabs(mcParticle->getPDG());
				//	cout<<"particle PDG is: "<<particlePDG<<endl;
					Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice


					if(particlePDG == 22)	//rad sa fotonima (PDG=22)
					{


						const EVENT::MCParticleVec & parents = mcParticle->getParents();
						if (parents.size() > 0)
						{
							IMPL::MCParticleImpl* parent = (IMPL::MCParticleImpl*)parents[0];

							if (parent->getPDG() == 25)
							{

								if(Pt>minPt1)
								{
								HiggsPhotons.push_back(temp);
								CounterHiggs++;
								/*vl.ptHiggsPhotons=temp.Pt();
								vl.energyHiggsPhotons=temp.E();
								vl.thetaHiggsPhotons = temp.Theta();*/
								}

							}
						}
					}

					if(abs(particlePDG) == 11)	//rad sa elektronima (PDG=11)
					{
						mcElectrons.push_back(temp);
					}


				}

			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("PandoraPFANewPFOs");/*("PandoraPFANewPFOs");*/


			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici



				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(recParticle->getType());
				//cout<<"particle PDG is: "<<particlePDG<<endl;



				if (abs(particlePDG) == 11)
				{
					recoElectrons.push_back(temp);
					//cout<<"elektron spoted"<<endl;

				}



				if(particlePDG == 22)	//rad sa fotonima (PDG=22)
				{


					photons.push_back(temp);

				} //pdg = 22



			}  // end of particle loop


	//ovde proveravam ugao izmedju fotona
			if (photons.size() == photonSize && HiggsPhotons.size()==photonSize)
			{
				Double_t ratioOfEnergy;
				Double_t anglephotons;
				for (int i = 0; i < (int)HiggsPhotons.size(); ++i)
				{
					for (int j = 0; j < (int)photons.size(); ++j)
					{
						anglephotons = HiggsPhotons[i].Angle(photons[j].Vect())* 180/M_PI;
						vl.anglePhotons = anglephotons ;

							ratioOfEnergy = photons[j].E()/HiggsPhotons[i].E();
							vl.photonsRatioEnergy = ratioOfEnergy ;

							eventList.Fill();

					}

				}
				//cout<<"brojac1:  "<<brojac1<<endl;
			}



			//cout <<"reco "<< recoElectrons.size()<<",    MC: "<< mcElectrons.size()<<endl;
			//Double_t ratioOfEnergyelectons;
			Double_t angleelectrons;
			for (int i = 0; i < (int)mcElectrons.size(); ++i)
			{

				for (int j = 0; j < (int)recoElectrons.size(); ++j)
				{
					angleelectrons = mcElectrons[i].Angle(recoElectrons[j].Vect())* 180/M_PI;
					vl.angleE = angleelectrons ;
					vl.electronsRatioEnergy = recoElectrons[j].E()/mcElectrons[i].E();
					vl.thetaMCe = mcElectrons[i].Theta();
					vl.thetaRecoe = photons[j].Theta();
					//eventList.Fill();

				}

			}




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

	TString rfName = "angleElectronsPhotons.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
