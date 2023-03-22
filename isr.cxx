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




	Float_t bsenergy = 0;
	Float_t bspt = 0;
	Float_t bstheta = 0;
	Float_t isrenergy = 0;
	Float_t isrpt = 0;
	Float_t isrtheta = 0;
	Float_t higgs_photon_pt = 0;
	Float_t reco_photon_pt = 0;




	TTree photons("photonTree", "ILD event list");
	photons.Branch("bsenergy", &bsenergy, "bsenergy");
	photons.Branch("bspt", &bspt, "bspt");
	photons.Branch("bstheta", &bstheta, "bstheta");
	photons.Branch("isrenergy", &isrenergy, "isrenergy");
	photons.Branch("isrpt", &isrpt, "isrpt");
	photons.Branch("isrtheta", &isrtheta, "isrtheta");

	TTree secondhighestpt("secondhighestpt", "ILD event list");
	secondhighestpt.Branch("higgs_photon_pt", &higgs_photon_pt, "higgs_photon_pt");
	secondhighestpt.Branch("reco_photon_pt", &reco_photon_pt, "reco_photon_pt");


 // elektronsko6i7Tree.Branch("pt_KG_elektrona_6", &pt_KG_elektrona_6, "pt_KG_elektrona_6");










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

			vector<TLorentzVector> HiggsPhotons;

			vector<TLorentzVector> recoElectrons;
			vector<TLorentzVector> mcElectrons;
			vector<TLorentzVector> mcPositrons;

			vector<TLorentzVector> bs_photon;
			vector<TLorentzVector> isr_photon;
			vector<Float_t> recPhotonsPt;
			bs_photon.clear();



				std::vector<std::string> colNames = *evt->getCollectionNames();
			//std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < (int)colNames.size(); ++i)
			{
			//std::cout << colNames[i] << endl;
			}



				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");
				IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("SelectedPandoraPFANewPFOs");


				for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
				{
					IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

					const double *p = recParticle->getMomentum(); // impuls čestice
					double e = recParticle->getEnergy();	//energija čestice
					temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
					Int_t particlePDG = recParticle->getType();

					if (particlePDG !=22) continue;

					IMPL::LCCollectionVec* mcPar = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticlesSkimmed");
						EVENT::LCRelation* link_to_rec_photon = 0;
	 					EVENT::MCParticle* pointer_to_mc_photon = 0;
	 					EVENT::ReconstructedParticle* pointer_to_rec_photon = 0;

	 				//	pointer_to_rec_photon =niz_fotoni[i];

	 					bool provera_isr = false;

	 					EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");
	 					for (int j = 0; j < links -> getNumberOfElements(); j++){

		 					EVENT::LCRelation* linki = (EVENT::LCRelation*) links -> getElementAt(j);
	 						EVENT::MCParticle* mcpj = (EVENT::MCParticle*) linki -> getTo();
							EVENT::ReconstructedParticle* reco = (EVENT::ReconstructedParticle*) linki -> getFrom();
						//	cout <<"našao čestice sa istom energijom."<<endl;
							if (reco ==0) continue;
							TLorentzVector temp (TVector3 (reco -> getMomentum()), reco -> getEnergy());

							//provera da li je foton ISR foton
							if (reco == recParticle){
							//	cout<<"linkovana cestica je: "<<rpi1->getType()<<endl;
								if (mcpj == 0 || mcpj->getPDG() !=22 || linki->getWeight()<0.985) continue;
								if(mcpj->getParents().size()== 0) continue;
								const EVENT::MCParticleVec & parent = mcpj -> getParents();


								//if (fabs(mcpj->getParents()[0]->getPDG())==22 &&  fabs(mcpj->getParents()[0]->getParents()[0]->getPDG())==11  && mcpj->getParents()[0]->getParents()[0]->getEnergy() > 1499.9 )
								if ( parent[0]->getPDG()==22 && parent[0]->getParents()[0]->getParents().size()==0 /*&& parent[0]->getParents()[0]->getParents()[0]->getEnergy()>1499*/)
								{
									provera_isr = true;
							//	cout<< "gen cestica ima koliko roditelja: "<<mcpj->getParents()[0]->getPDG()<<endl;
							//	cout << "energija roditelja je : "<< mcpj->getParents()[0]->getParents()[0]->getEnergy()<<endl;
								//	cout<< "baba je: "<< parent[0]->getParents()[0]->getPDG()<<", energija babe je: " <<parent[0]->getParents()[0]->getEnergy()<<endl;
							//		cout<<"_______________________________"<<endl;

								}//end of provera ISR foton
							}//end of if (reco == recParticle){

	 					}// end of LINKS
						if (provera_isr == true ) continue;

						recPhotonsPt.push_back(temp.Pt());



				}// end of recParticle loop



				for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
				{
					IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

					//if (mcParticle->getGeneratorStatus() != 1) continue;
					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

					const double *p = mcParticle->getMomentum(); // impuls čestice
					double e = mcParticle->getEnergy();	//energija čestice
					temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

					Int_t particlePDG = mcParticle->getPDG();
				//	cout<<"particle PDG is: "<<particlePDG<<endl;
					Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice

					const EVENT::MCParticleVec & parent = mcParticle ->getParents(); // foton potomak


					if (particlePDG == 11 && temp.E()>1499 && mcParticle->getParents().size() == 0){
						//nasao elektron ili pozitron koji se sudara
						const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters(); // foton potomak
					//	cout <<"pocetni elektron ima koliko cerki: "<< daughter.size()<<", redni broj događaja: "<< brojac<<endl;
						for (int j = 0; j < (int) daughter.size(); j++){
							if (daughter[j]->getPDG() == 22){ // ovo su beamstrahlung fotoni
								TLorentzVector bs_temp (TVector3 (daughter[j] -> getMomentum()), daughter[j] -> getEnergy());
					//			cout << "pdg cerke: "<<daughter[j]->getPDG()<<endl;
								bs_photon.push_back(bs_temp);
							}
						}
					}// end of (particlePDG== 11 && temp.E()>1499)


					if(particlePDG == 11 && temp.E()>1499 && mcParticle->getParents().size() == 0){
						const EVENT::MCParticleVec & daughter_isr = mcParticle -> getDaughters(); // foton potomak
						const EVENT::MCParticleVec & granddaughter_isr = daughter_isr[0] -> getDaughters(); // foton potomak

						for (int k = 0; k < (int) granddaughter_isr.size(); k++){
						if (granddaughter_isr[k]->getPDG() == 22){ // ovo su beamstrahlung fotoni
							TLorentzVector isr_temp (TVector3 (granddaughter_isr[k] -> getMomentum()), granddaughter_isr[k] -> getEnergy());
				//			cout << "pdg cerke: "<<daughter[j]->getPDG()<<endl;
							isr_photon.push_back(isr_temp);
						}
					}
					}


					//uzimamo foton sa nižim pT
					if (particlePDG == 25 && mcParticle->getDaughters().size()==2){
					//	cout << "našao Higgsa"<<endl;
						const EVENT::MCParticleVec & higgs_daughter = mcParticle -> getDaughters(); // foton potomak
						TLorentzVector higgs_photon1 (TVector3 (higgs_daughter[0] -> getMomentum()), higgs_daughter[0] -> getEnergy());
						TLorentzVector higgs_photon2 (TVector3 (higgs_daughter[1] -> getMomentum()), higgs_daughter[1] -> getEnergy());

						if (higgs_photon1.Pt() < higgs_photon2.Pt() ){
							higgs_photon_pt = higgs_photon1.Pt();
						}else{
							higgs_photon_pt = higgs_photon2.Pt();
						}


					}



				}// end of MC particle loop

				if (bs_photon.size()==2 && isr_photon.size()==2){
			//		cout << "bs fotona: "<< bs_photon.size()<<endl;
				//	cout << "isr fotona: "<< isr_photon.size()<<endl;
					for (int l = 0 ; l < 2 ; l++){
					 bsenergy=bs_photon[l].E();
					 isrenergy=isr_photon[l].E();
					 bspt=bs_photon[l].Pt();
					 isrpt=isr_photon[l].Pt();
					 isrtheta=isr_photon[l].Theta()*180/M_PI;
					 bstheta=bs_photon[l].Theta()*180/M_PI;
					 photons.Fill();


					}

				}

				sort(recPhotonsPt.begin(), recPhotonsPt.end(), greater<int>());

				reco_photon_pt = recPhotonsPt[1];

				secondhighestpt.Fill();



		} // End of event loop





		lcReader->close();

	} // End of file loop






	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	photons.Write();
	secondhighestpt.Write();




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

	TString rfName = "isr.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
