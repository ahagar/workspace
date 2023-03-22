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

const Double_t minPt1 = 10;			   //minimalna energija fotona
const Double_t minPt2 = 10;			   //minimalna energija fotona
const Double_t minPt3 = 11;			   //minimalna energija fotona
const Double_t minPt4 = 12;			   //minimalna energija fotona

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


Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif



		TH1F thetaEcalBarrel ("thetaEcalBarrel", "#theta_{ECAL Barrel}; #theta_{hit}", 180, 0, 180);
		TH1F thetaEcalEndCap("thetaEcalEndCap", "#theta_{ECAL End Cap}; #theta_{hit}", 180, 0, 180);





	TTree eventList("eventsSignal", "ILD event list");
	varListGoran vl; /* SL specific */

	eventList.Branch("hCalSum", &(vl.hCalSum), "hCalSum");
	eventList.Branch("eCalSum", &(vl.eCalSum), "eCalSum");
	eventList.Branch("hCalEndCap", &(vl.hCalEndCap), "hCalEndCap");
	eventList.Branch("eCalEndCap", &(vl.eCalEndCap), "eCalEndCap");
	eventList.Branch("eCalBarel", &(vl.eCalBarel), "eCalBarel");
	eventList.Branch("hCalBarel", &(vl.hCalBarel), "hCalBarel");
	eventList.Branch("eCalOther", &(vl.eCalOther), "eCalOther");
	eventList.Branch("hCalOther", &(vl.hCalOther), "hCalOther");
	eventList.Branch("calSum", &(vl.calSum), "calSum");
	eventList.Branch("hCalSum_NotRepr", &(vl.hCalSum_NotRepr), "hCalSum_NotRepr");
	eventList.Branch("eCalSum_NotRepr", &(vl.eCalSum_NotRepr), "eCalSum_NotRepr");
	eventList.Branch("hCalEndCap_NotRepr", &(vl.hCalEndCap_NotRepr), "hCalEndCap_NotRepr");
	eventList.Branch("eCalEndCap_NotRepr", &(vl.eCalEndCap_NotRepr), "eCalEndCap_NotRepr");
	eventList.Branch("eCalBarel_NotRepr", &(vl.eCalBarel_NotRepr), "eCalBarel_NotRepr");
	eventList.Branch("hCalBarel_NotRepr", &(vl.hCalBarel_NotRepr), "hCalBarel_NotRepr");
	eventList.Branch("eCalOther_NotRepr", &(vl.eCalOther_NotRepr), "eCalOther_NotRepr");
	eventList.Branch("hCalOther_NotRepr", &(vl.hCalOther_NotRepr), "hCalOther_NotRepr");
	eventList.Branch("calSum_NotRepr", &(vl.calSum_NotRepr), "calSum_NotRepr");







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


			IMPL::LCCollectionVec* eCalBarrelHits = (IMPL::LCCollectionVec*)evt->getCollection("ECALBarrel_Reprocess");
			IMPL::LCCollectionVec* eCalEndcapHits = (IMPL::LCCollectionVec*)evt->getCollection("ECALEndcap_Reprocess");
			IMPL::LCCollectionVec* eCalOtherHits = (IMPL::LCCollectionVec*)evt->getCollection("ECALOther_Reprocess");


			IMPL::LCCollectionVec* hCalBarrelHits = (IMPL::LCCollectionVec*)evt->getCollection("HCALBarrel_Reprocess");
			IMPL::LCCollectionVec* hCalEndcapHits = (IMPL::LCCollectionVec*)evt->getCollection("HCALEndcap_Reprocess");
			IMPL::LCCollectionVec* hCalOtherHits = (IMPL::LCCollectionVec*)evt->getCollection("HCALOther_Reprocess");

			IMPL::LCCollectionVec* eCalBarrelHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("ECALBarrel");
			IMPL::LCCollectionVec* eCalEndcapHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("ECALEndcap");
			IMPL::LCCollectionVec* eCalOtherHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("ECALOther");


			IMPL::LCCollectionVec* hCalBarrelHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("HCALBarrel");
			IMPL::LCCollectionVec* hCalEndcapHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("HCALEndcap");
			IMPL::LCCollectionVec* hCalOtherHits_NotRepr = (IMPL::LCCollectionVec*)evt->getCollection("HCALOther");



			Double_t eCalSum = 0;
			Double_t eCalBarel = 0;
			Double_t eCalEndCap = 0;
			Double_t eCalOther = 0;

			Double_t hCalSum = 0;
			Double_t hCalBarel = 0;
			Double_t hCalEndCap = 0;
			Double_t hCalOther = 0;

			// reprocesed
			eCalBarel += SumCalorimeterHitEnergies(eCalBarrelHits);
			eCalEndCap += SumCalorimeterHitEnergies(eCalEndcapHits);
		//	eCalOther += SumCalorimeterHitEnergies(eCalOtherHits);
			eCalSum = eCalBarel + eCalEndCap ;

			hCalBarel += SumCalorimeterHitEnergies(hCalBarrelHits);
			hCalEndCap += SumCalorimeterHitEnergies(hCalEndcapHits);
		//	hCalOther += SumCalorimeterHitEnergies(hCalOtherHits);
			hCalSum = hCalBarel + hCalEndCap ;



			Double_t calSum = eCalSum + hCalSum;

			//not Reprocessed
			Double_t eCalSum_NotRepr = 0;
			Double_t eCalBarel_NotRepr = 0;
			Double_t eCalEndCap_NotRepr = 0;
			Double_t eCalOther_NotRepr = 0;

			Double_t hCalSum_NotRepr = 0;
			Double_t hCalBarel_NotRepr = 0;
			Double_t hCalEndCap_NotRepr = 0;
			Double_t hCalOther_NotRepr = 0;



			eCalBarel_NotRepr += SumCalorimeterHitEnergies(eCalBarrelHits_NotRepr);
			eCalEndCap_NotRepr += SumCalorimeterHitEnergies(eCalEndcapHits_NotRepr);
			eCalOther_NotRepr += SumCalorimeterHitEnergies(eCalOtherHits_NotRepr);
			eCalSum_NotRepr= eCalBarel_NotRepr + eCalEndCap_NotRepr +eCalOther_NotRepr;

			hCalBarel_NotRepr += SumCalorimeterHitEnergies(hCalBarrelHits_NotRepr);
			hCalEndCap_NotRepr += SumCalorimeterHitEnergies(hCalEndcapHits_NotRepr);
			hCalOther_NotRepr += SumCalorimeterHitEnergies(hCalOtherHits_NotRepr);
			hCalSum_NotRepr = hCalBarel_NotRepr + hCalEndCap_NotRepr + hCalOther_NotRepr;

			Double_t calSum_NotRepr = eCalSum_NotRepr + hCalSum_NotRepr;


			thetaEcalBarrel.Fill(eCalBarel_NotRepr);
			thetaEcalEndCap.Fill(eCalEndCap_NotRepr);



			vl.eCalSum = eCalSum;
			vl.hCalSum = hCalSum;
			vl.calSum = calSum;
			vl.eCalBarel = eCalBarel;
			vl.eCalEndCap = eCalEndCap;
			vl.eCalOther = eCalOther;
			vl.hCalBarel = hCalBarel;
			vl.hCalEndCap = hCalEndCap;
			vl.hCalOther = hCalOther;

			vl.eCalSum_NotRepr = eCalSum_NotRepr;
			vl.hCalSum_NotRepr = hCalSum_NotRepr;
			vl.calSum_NotRepr = calSum_NotRepr;
			vl.eCalBarel_NotRepr = eCalBarel_NotRepr;
			vl.eCalEndCap_NotRepr = eCalEndCap_NotRepr;
			vl.eCalOther_NotRepr = eCalOther_NotRepr;
			vl.hCalBarel_NotRepr = hCalBarel_NotRepr;
			vl.hCalEndCap_NotRepr = hCalEndCap_NotRepr;
			vl.hCalOther_NotRepr = hCalOther_NotRepr;




			//cout << "ukupna enerrgija u kalorimetrima je : "<< hCalSum << ", " << hCalBarel << ", "<< hCalEndCap <<",  "<< hCalOther<<endl;


			eventList.Fill();


		} // End of event loop





		lcReader->close();

	} // End of file loop



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");


	eventList.Write();





	thetaEcalBarrel.Write();
	thetaEcalEndCap.Write();



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
