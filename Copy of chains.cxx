/***********************************************************
 *
 * 	Read reconstructed photons from a  .slcio file, and draws histograms for different variables
 *
 *  Author:Goran Kačarević
 *  15.10.2015.
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



//function that calculates summ of energies from given collection
Double_t SumCalorimeterHitEnergies(EVENT::LCCollection* calHitCollection)
{
	Double_t result = 0;
	Double_t theta= 0;

	for (Int_t i = 0; i < calHitCollection->getNumberOfElements(); i++)
	{
		EVENT::CalorimeterHit* calHit = (EVENT::CalorimeterHit*)calHitCollection->getElementAt(i);
		const float *pos = calHit->getPosition();
		//calHit->getCellID0()
		// da li je pos NULL?
		double x = pos[0];
		double y = pos[1];
		double z = pos[2];
		theta =  TMath::ATan2(sqrt(x*x+y*y),z)*180/M_PI;

		result += calHit->getEnergy();

	}

	return theta;
}
EVENT::MCParticle*  chains (EVENT::MCParticle*  finalElectron)
{


return 0;
}


Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif



		TH1F thetaEcalBarrel ("thetaEcalBarrel", "#theta_{ECAL Barrel}; #theta_{hit}", 180, 0, 180);
		TH1F thetaEcalEndCap("thetaEcalEndCap", "#theta_{ECAL End Cap}; #theta_{hit}", 180, 0, 180);



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


		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");

		   	for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		   	{

				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				if (mcParticle -> getPDG() == 25 && mcParticle -> getParents()[0] -> getPDG() == 11)
				{
					const EVENT::MCParticleVec & parent = mcParticle -> getParents();
					for (Int_t h = 0; h < parent.size(); h++)
					{
						//TLorentzVector e_starac1, p_starac1; //četvorovektor u koji sakupljamo informacije o svakoj čestici

						Int_t PDG_starac = parent[h] -> getPDG();

						const double *pgen = parent[h] -> getMomentum(); // impuls čestice
						double egen = parent[h] -> getEnergy();	//energija čestice

						// cout << "Energija: " << parent[h] -> getEnergy() << endl;

					}//end za roditelje 163


					// KG traže se potomci Higsa, koristi se za traženje signala u liniji 1617
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					//const EVENT::MCParticleVec & parent = mcParticle -> getParents(); // Traži se roditelj Higsa
					const EVENT::MCParticleVec & cerke = parent[0] -> getDaughters(); // Traže se potomci roditelja Higsa

					for (int l = 0; l < (int) cerke.size(); l++) // Petlja po potomcima Higsa koji su naš signal
					{

						int PDG_cerka = cerke[l] -> getPDG();
						if (cerke[l]->getGeneratorStatus() == 1) cout <<"elektron 6 ili 7 je finalan!!! "<<endl;

						TLorentzVector cerkaTemp (TVector3 (cerke[l] -> getMomentum()), cerke[l] -> getEnergy());

						const double *pgen = cerke[l] -> getMomentum(); // impuls čestice
						double egen = cerke[l] -> getEnergy();	//energija čestice

						if (abs (PDG_cerka) == 11){
							//const EVENT::MCParticleVec & finEl =	cerke[l]->getDaughters();
							if(cerke[l]->getDaughters().size()==1){
								if(cerke[l]->getPDG()==11 && cerke[l]->getGeneratorStatus()==1 ){
									cout<<"nasao sam finalni elektron!!!"<<endl;
								}

							if (cerke[l]->getGeneratorStatus()!=1){
								cout <<"naći način da se sačuva čestica i vrati na početak fje"<<endl;
								}
							}
							if (cerke[l]->getDaughters().size()>1){
								for(int i=0; i < (int) cerke[l]->getDaughters().size(); i++){
									if(cerke[l]->getDaughters()[i]->getPDG()==11 && cerke[l]->getDaughters()[i]->getGeneratorStatus()==1){
										cout << "nasao sam finalnu cesticu"<<endl;
									}else{
										cout<<"nastavi dalje"<<endl;
									}
								}
							}

						}//end probe funkcije


						// Ako je potomak elektron 6
						if (PDG_cerka == 11) // promenljive za elektron 6
						{



						}

						// Ako je potomak pozitron 7
						if (PDG_cerka == -11) // promenljive za pozitron 7
						{

 							//N_KG_pozitrona_7++;


						}

					}


				}///end if za čestice koje se sudaraju


			}

		} // End of event loop





		lcReader->close();

	} // End of file loop



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

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

	TString rfName = "barelendcap.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
