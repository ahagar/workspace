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
const double minTheta = 0;			   //minimalni polarni ugao
const double maxTheta = 180;			   //miximalni polarni ugao
const double minInvMass = 0;
const double maxInvMass = 150000;



Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif

	gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/test/test/CLICdpStyle.C");


	//Deklaracija histograma
	TH1F histoPtMax ("HighestPhotonPt", "Max Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPt2Max ("2ndHighestPhotonPt", "2nd Max Photon Pt; Pt (GeV)", 1500, 0, 1500);


	Float_t ptPhoton, ePhoton ; //Deklaracija promenljivih koje se koriste u drvetu

	TTree eventList("eventsSignal", "ILD event list"); //Deklaracija drveta
	eventList.Branch("ptPhoton", &ptPhoton, "pT_{GeV}");//deklaracija grana u drvetu
	eventList.Branch("ePhoton", &ePhoton, "ePhoton");


	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	  Int_t totalCutEvents = 0;
	  Int_t totalEvents = 0;

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

		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja po fajlu
		cout << "Broj dogadjaja po fajlu je  : " << brojDogadjaja << endl;




		// Prolazimo kroz svaki događaj
		EVENT::LCEvent* evt = 0;
		while((evt = lcReader->readNextEvent()) != 0)
		{
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");//koju kolekciju uzima

            totalEvents++; // brojač koji se uvećava kroz svaki događaj

			vector<TLorentzVector> photons;//deklaracija niza koji prima kao članove promenljive TLorentzVector(4-vektor)
			vector <Double_t> PtPhotons; // deklaracija niza koji prima kao članove pt fotona, tipa double
			vector <Double_t> tng;
			vector <Double_t> angle;
			vector <Double_t> thetaPhotons;


			//petlja koja učitava česticu u datoj kolekciji
			for (Int_t i = 0; i < mcParticles->getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);//uzima i-tu česticu


				TLorentzVector temp; //deklaracija četvorovektora kojim možemo dalje da manipulišemo

				const double *p = mcParticle->getMomentum(); // trovektor impulsa
				double e = mcParticle->getEnergy(); // energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e); //upis impulsa i energije u četvorovektor
				if (mcParticle->getGeneratorStatus() != 1) continue; // proverava da li je čestica finalna

				Int_t ugao = temp.Theta() * 180/M_PI; // računa polarni ugao čestice

				Int_t particlePDG = fabs(mcParticle->getPDG()); // PDG čestice

				if(particlePDG == 22)	// fotoni (PDG == 22)
				{
					PtPhotons.push_back(temp.Pt()); // upis u niz
					sort(PtPhotons.begin(), PtPhotons.end(), greater<int>()); // sortiranje ččestica u nizu od najveće vrednosti ka najvanjoj

					Double_t tangens = temp.Pt()/temp.E();
					tng.push_back(tangens);
					sort(tng.begin(), tng.end(), greater<int>());

					ePhoton=e;
					ptPhoton=temp.Pt();

					eventList.Fill();

					if (ugao > 90 )
						{
							ugao = 180- ugao;
						}
					angle.push_back(ugao);
					sort(angle.begin(), angle.end(), greater<int>());



					if(temp.Pt() > minPt && temp.Theta()*180/M_PI > minTheta && temp.Theta()*180/M_PI < maxTheta)
					{
						photons.push_back(temp);
						PtPhotons.push_back(temp.Pt());
						sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());


						Double_t promenljiva = TMath::ATan(temp.Pt()/fabs(temp.Pz()))*180/M_PI;//promenljiva je theta
						thetaPhotons.push_back(promenljiva);
						sort(thetaPhotons.begin(), thetaPhotons.end(), greater<int>());
					}

				}
			/*	else if (particlePDG != 22)
				{
				    if(temp.Pt() > 200000)
					    leptonFound = true;
				}*/
			}  // end of particle loop

			//if (leptonFound) continue;
		//	sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());

            if (photons.size() > 1)
            	{
            	//totalCutEvents++;
				vector<CandidateData> candidates;


				histoPtMax.Fill(PtPhotons[0]);




				//cout << "theta je: " << thetaPhotons[1] << endl;


				for (int i = 0; i < (int)photons.size() - 1; ++i)
					{
					bool eventFound = false;

						for (int j = i + 1; j < photons.size(); ++j)
						{
							CandidateData current(photons[i], photons[j]);//uzima dva fotona

							TLorentzVector pair = current.Higgs();//sabira četvorovektore


							if (minInvMass < pair.M() && pair.M() < maxInvMass)//uslov da bi par bio higs kandidat
							{
								candidates.push_back(current);
				            	totalCutEvents++;

								histoPt2Max.Fill(PtPhotons[1]);
								eventFound = true;
								break;

							}
						}
						if (eventFound) break;
					}

            	}





		} // End of event loop



		lcReader->close();



	} // End of file loop
	cout << "minimalni Pt:     " << minPt << endl;

	cout << "Ukupan broj dogadjaja:     " << totalEvents << endl;
    cout << "Post-cut events in file:  " << totalCutEvents << endl;
    cout << "Percentage:               " << 100.0 * (double)totalCutEvents/totalEvents << endl;

	ofstream results ;
	results.open("histoNoCuts.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt << ", " << minTheta<<" < Theta <  "<< maxTheta  <<", "<< minInvMass<<" < M < " << maxInvMass <<  ", efikasnost = "<< 100.0 * (double)totalCutEvents/totalEvents  <<endl;
	results.close();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");


	histoPtMax.Write();
	histoPt2Max.Write();
	eventList.Write();


	rootFile.Close();


	TCanvas c1;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c1.SetCanvasSize(1000,650);
	c1.Divide(1,1,0.01,0.01);
	c1.cd(1);
//	c1.GetPad(1)->SetLogy();
	histoPtMax.Draw();
	c1.Print("noCutshistoMaxPt.eps");



	TCanvas c2;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c2.SetCanvasSize(1000,650);
	c2.Divide(1,1,0.01,0.01);
	c2.cd(1);
	//c2.GetPad(1)->SetLogy();
	histoPt2Max.Draw();
	c2.Print("noCutshisto2ndMax.eps");





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

	TString rfName = "histoNoCuts.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
