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
#include <UTIL/LCTOOLS.h>
#include <Exceptions.h>
#endif

#include "stdlib.h"
#include "CandidateData.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include "varList.h"


using namespace std;

const Double_t minPt1 = 0;			   //minimalna energija fotona
const Double_t minPt2 = 12;			   //minimalna energija fotona
const Double_t minPt3 = 15;			   //minimalna energija fotona
const Double_t minPt4 = 20;			   //minimalna energija fotona

const double mH = 126.0;               //Higgs mass
const Double_t coneAngle = 2.5;		   //ugao konusa
const Double_t maxConeEnergy = 20 ;	   //maksimalna energija konusa
const Double_t leptonPt = 200000;		   //energija leptona posle koje zadovoljavaju uslov leptonFound
const Double_t EnergyCenterMass = 3000; // energija u sistemu centra masa, koju koristim za missing energy
const Double_t minInvMass = 0;
const Double_t maxInvMass = 1400000;


//funkcija koja proverava koji par fotona je najbolji kandidat za rekonstrukciju Higgsa
bool IsBetterHiggsCandidate(CandidateData a, CandidateData b)
{
	TLorentzVector higgsA = a.Higgs();
	TLorentzVector higgsB = b.Higgs();

	Double_t distanceA = fabs(higgsA.M() - mH);
	Double_t distanceB = fabs(higgsB.M() - mH);
	if (distanceA < distanceB)
	{
		return true;
	}

	return false;
}

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif


	//gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/crtanje/histogrami/CLICdpStyle.C");
	//histogrami za različite kinematičke varijable


	TTree photonTree ("photonTree", "Generator particle tree");
	varListGoran vl;
	photonTree.Branch("CandidateM", &(vl.CandidateM), "CandidateM/F");

	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{#gamma#gamma} (GeV)", 1600, -100, 1500);



	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	Int_t totalCutEvents = 0;
	Int_t totalEvents = 0;
	Int_t neUslov = 0;


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
		Int_t numCutEvents = 0; 	//broj dogadjaja posle određenih cut-ova

		std::vector <TLorentzVector> allPhotons; //vektor koji prikuplja sve fotone
		Int_t brojac =0;


		//Int_t hardPhotons = 0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
			brojac++;

			vector <TLorentzVector> sviFotoni; //vektor koji prikuplja sve fotone
			vector <float> kandidati;


				std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");/*("PandoraPFANewPFOs"); MCParticles*/


			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
			{
			//	IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

				//if (recParticle->getGeneratorStatus() != 1) continue;

				if (mcParticle->getPDG() == 25 && mcParticle->getDaughters()[0]->getPDG() ==22){
					//cout<<"nasao sam Higsa!!!"<<endl;
					TLorentzVector photon1;
					TLorentzVector photon2;
					TLorentzVector higs;

					const double *pH = mcParticle->getMomentum(); // impuls čestice
					double eH = mcParticle->getEnergy();	//energija čestice


					const double *p1 = mcParticle->getDaughters()[0]->getMomentum(); // impuls čestice
					const double *p2 = mcParticle->getDaughters()[1]->getMomentum(); // impuls čestice
					double e1 = mcParticle->getDaughters()[0]->getEnergy();	//energija čestice
					double e2 = mcParticle->getDaughters()[1]->getEnergy();	//energija čestice
					photon1.SetPxPyPzE(p1[0],p1[1], p1[2],e1);
					photon2.SetPxPyPzE(p2[0],p2[1], p2[2],e2);
					higs = photon1+photon2;
				//	higs.SetPxPyPzE(pH[0], pH[1],pH[2],eH);

					TVector3 BoostToHiggs = -(higs.BoostVector() );
					higs.Boost(BoostToHiggs);
					photon1.Boost(BoostToHiggs);
					photon2.Boost(BoostToHiggs);

					TVector3 ph1_3 (photon1.X(), photon1.Y(), photon1.Z());
					TVector3 ph2_3 (photon2.X(), photon2.Y(), photon2.Z());
					Double_t angleP = ph1_3.Angle(ph2_3);
					cout << "ugao izmedju 2 fotona u Higsovom sistemu reference je: "<< angleP<<endl;
					cout <<"vrednosti Higsa u Higsovom sistemu reference su"<<endl;
					higs.Print();
					cout <<"vrednosti fotona 1 u Higsovom sistemu reference su"<<endl;
					photon1.Print();
					cout <<"vrednosti fotona2 u Higsovom sistemu reference su"<<endl;
					photon2.Print();
					cout <<"__________________________________________________________________________"<<endl;




				}




				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici


				const double *p = mcParticle->getMomentum(); // impuls čestice
				double e = mcParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

			//	Int_t particlePDG = fabs(recParticle->getType());
				Int_t particlePDG = fabs(mcParticle->getPDG());

			//	cout<<"particle PDG is: "<<particlePDG<<endl;


			//	cout << "testiranje1 " <<  endl;

				//if(particlePDG ==11 &&(temp.Theta()*180/M_PI > 8 && temp.Theta()*180/M_PI <172) && temp.E()>500.) cout<< "ugao elektrona je: "<< temp.Theta()*180/M_PI<<endl;

				if(particlePDG == 22 && temp.Pt() > 0)	//rad sa fotonima (PDG=22)
				{

					sviFotoni.push_back(temp);
				}


			}  // end of particle loop

		//	cout << "broj fotona za Pt > 15 po događaju je: " <<sviFotoni.size()<<endl;


			for (int m = 0; m < (int)sviFotoni.size()-1; m++){
				for(int l = 1; l < (int)sviFotoni.size();l++){
					float invM1 = (sviFotoni[m]+sviFotoni[l]).M();
					vl.CandidateM = invM1;
					photonTree.Fill();
					//histoCandidateM.Fill(invM1);
				if (invM1>100 && invM1<150)	kandidati.push_back(invM1);
					//if (invM1 > 100) cout << "test"<< endl;

				}
			}

			if (kandidati.size() ==0){
				//	cout <<"ne postoji kandidat u događaju"<<endl;
						neUslov++;
					//	cout<< "broj fotona u događaju je:" << sviFotoni.size()<<endl;
						for (int i = 0; i < (int) sviFotoni.size()-1; i++){
							for (int j = 1; j < (int) sviFotoni.size(); j++){
						//		cout <<"inv masa para fotona je : "<< (sviFotoni[i]+sviFotoni[j]).M()<<endl;
								float invM2 = (sviFotoni[i]+sviFotoni[j]).M();
								histoCandidateM.Fill(invM2);
							//	candidateM = invM2;
							//	if (invM2 < -1) cout <<"candidateM= "<<invM2<< endl;
					//			photonTree.Fill();

							}//end of j

						}//end of i allPhotons
			}


		} // End of event loop


		totalCutEvents+=numCutEvents;//ukupan broj dogadjaja koji prodje cutove iz svih fajlova
		totalEvents+=brojDogadjaja; //ukupan broj dogadjaja iz svih fajlova


		cout<<"broj događaja koji ne prodje uslov je: "<<neUslov<< endl;

		lcReader->close();

	} // End of file loop


	cout<<"ukupan broj događaja je : "<< totalEvents<< endl;




	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	//eventList.Write();

	photonTree.Write();
	histoCandidateM.Write();



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

	TString fName = "h_nunu_dst_6265_";
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "MCParticlesSkimmed.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
