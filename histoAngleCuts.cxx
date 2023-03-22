//Testira broj dogadjaja koji prodju cutovea ako imamo 2 fotona koji ispunjavaju uslov za PT

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
#include <UTIL/LCTOOLS.h>
#include <Exceptions.h>
#endif

#include "stdlib.h"
#include "CandidateData.h"
#include <sstream>
#include <iostream>
#include <iterator>

using namespace std;

const Double_t minPt = 0;			   //minimalna energija fotona
const double mH = 126.0;               //Higgs mass
const double minTheta = 5;			   //minimalni polarni ugao
const double maxTheta = 175;			   //miximalni polarni ugao
const double minInvMass = 0;
const double maxInvMass = 150000;



Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif

	//gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/test/test/CLICdpStyle.C");
	Float_t thetaE = 0;
	Float_t thetaP = 0;	
	Float_t elecE = 0;
	Float_t posE = 0;
	Float_t isrE = 0;
	Float_t isrPt = 0;
	Float_t isrTheta = 0;



	
		TTree eventList("treeE", "ILD event list");
		eventList.Branch("thetaE", &thetaE, "thetaE");
		eventList.Branch("thetaP", &thetaP, "thetaP");
		eventList.Branch("elecE", &elecE, "elecE");
		eventList.Branch("posE", &posE, "posE");
		eventList.Branch("isrE", &isrE, "isrE");
		eventList.Branch("isrPt", &isrPt, "isrPt");
		eventList.Branch("isrTheta", &isrTheta, "isrTheta");

		
		



	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	  Int_t totalCutEvents = 0;
	  Int_t totalEvents = 0;
	  int br = 0 ;
	  bool electron, positron;

	//petlja koja iščitava .slcio podatke ukoliko ima više fajlova za jedan process
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{

		cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob); //3i
				try
				{
					lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));// 3i
				}
				catch(lcio::IOException &ex)
				{
					cout << ". Could not open.\n"; // Exception " << ex.what() << endl;
					continue;
				}
	/*	cout << "Opening " << Form("%s.%.3i.slcio", fName.Data(), iJob);
		try
		{
			lcReader->open(Form("%s.%.3i.slcio", fName.Data(), iJob));
		}
		catch(lcio::IOException &ex)
		{
			cout << ". Could not open.\n";
			continue;
		}*/
		cout << "... Reading ...\n";

		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja
		cout << "Broj dogadjaja po fajlu je  : " << brojDogadjaja << endl;




		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while((evt = lcReader->readNextEvent()) != 0 )
		{
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");

    //       bool leptonFound = false;
			electron = false;
			positron = false;
            totalEvents++;
            if (totalEvents >5)continue;


			for (Int_t i = 0; i < mcParticles->getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);
				Int_t particlePDG = mcParticle->getPDG();
				//cout<<"simple UID: "<<mcParticle->simpleUID()<<endl;
				TLorentzVector temp;

				const double *p = mcParticle->getMomentum();
				double e = mcParticle->getEnergy();
				temp.SetPxPyPzE(p[0], p[1], p[2], e);

				if (particlePDG ==22){
					isrE=temp.E();
					isrPt=temp.Pt();
					isrTheta=temp.Theta();
				}

		//		if (mcParticle->getGeneratorStatus() != 1) continue;
			//	if (e < 100) continue;

				double ugao = temp.Theta() * 180/M_PI;
				

				
				if (particlePDG == 11)
				{

				thetaE = ugao;
				elecE = e;

				}
								
				if (particlePDG == -11)
				{

				thetaP = ugao;
				posE = e;

				}			//		cout << "usao u mc:"<<mcParticle -> getParents()[0] -> getPDG())<<endl;
			/*
				if (mcParticle -> getPDG() == 25 && fabs(mcParticle -> getParents()[0] -> getPDG()) == 11 )
				{
					// KG traže se potomci Higsa, koristi se za traženje signala u liniji 1617
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					const EVENT::MCParticleVec & parent = mcParticle -> getParents(); // Traži se roditelj Higsa
					const EVENT::MCParticleVec & cerke = parent[0] -> getDaughters(); // Traže se potomci roditelja Higsa
					cout << "usao u mc:"<<endl;

					for (int l = 0; l < (int) cerke.size(); l++) // Petlja po potomcima Higsa koji su naš signal
					{
						int PDG = cerke[l] -> getPDG();

						TLorentzVector cerkaTemp (TVector3 (cerke[l] -> getMomentum()), cerke[l] -> getEnergy());
						cout << "e: "<< cerkaTemp.E();

						// Ako je potomak elektron 6
						if (PDG == 11) // promenljive za elektron 6
						{
    
				thetaE = cerkaTemp.Theta();
				elecE = cerkaTemp.E();

			
						}

						// Ako je potomak pozitron 7
						if (PDG == -11) // promenljive za pozitron 7
						{

				thetaP = cerkaTemp.Theta();
				posE = cerkaTemp.E();;

							
						}

					}

				
				} // Kraj IF uslova za Higsov bozon koji se raspada na ma koja dva kvarka
*/

				eventList.Fill();


			}  // end of particle loop
			



		} // End of event loop



		lcReader->close();



	} // End of file loop
	cout << "minimalni Pt:     " << minPt << endl;

	cout << "Ukupan broj dogadjaja:     " << totalEvents << endl;
 //   cout << "Post-cut events in file:  " << totalCutEvents << endl;
   // cout << "Percentage:               " << 100.0 * (double)totalCutEvents/totalEvents << endl;


	/*ofstream results ;
	results.open("histoAngleCuts.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt << ", " << minTheta<<" < Theta <  "<< maxTheta  <<", "<< minInvMass<<" < M < " << maxInvMass <<  ", efikasnost = "<< 100.0 * (double)totalCutEvents/totalEvents  <<endl;
	results.close();
*/
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

	TString fName = "whizard";//"h_nunu_dst_6265_";
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "histoAngleCuts.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
