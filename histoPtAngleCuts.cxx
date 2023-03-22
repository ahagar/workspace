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

const Double_t minPt = 10;			   //minimalna energija fotona
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

	TH1F histoPtMax ("HighestPhotonPt", "Max Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPt2Max ("2ndHighestPhotonPt", "2nd Max Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPtMaxZoom ("Max Photon Pt Zoomed", "Max Photon Pt Zoomed; Pt (GeV)", 50, 0, 50);
	TH1F histoPt2MaxZoom ("2nd Max Photon Pt Zoomed", "2nd Max Photon Pt Zoomed; Pt (GeV)", 150, 0, 150);
	TH1F histoPhotonPt ("Photon Pt", "Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPhotonTheta ("Photon Theta", "Photon Theta", 90, 0, 90);
	TH1F histoPhotonTg ("tangens", "tangens", 10000, 0, 0.1);
	TH1F histo2ndHighestThetaPhotons ("2ndHighestThetaPhotons", "2nd Highest Theta of Photons; theta* (deg)", 90, 0 , 90);


	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	  Int_t totalCutEvents = 0;
	  Int_t totalEvents = 0;

	//petlja koja iščitava .slcio podatke ukoliko ima više fajlova za jedan process
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{

		cout << "Opening " << Form("%s.%.3i.slcio", fName.Data(), iJob);
				try
				{
					lcReader->open(Form("%s.%.3i.slcio", fName.Data(), iJob));
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
		while((evt = lcReader->readNextEvent()) != 0)
		{
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");

    //       bool leptonFound = false;
            int photonsPassed = 0;
            totalEvents++;

			vector<TLorentzVector> photons;
			vector <Double_t> PtPhotons;
			vector <Double_t> tng;
			vector <Double_t> angle;
			vector <Double_t> thetaPhotons;


			for (Int_t i = 0; i < mcParticles->getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);


				TLorentzVector temp;

				const double *p = mcParticle->getMomentum();
				double e = mcParticle->getEnergy();
				temp.SetPxPyPzE(p[0], p[1], p[2], e);
				if (mcParticle->getGeneratorStatus() != 1) continue;

				Int_t ugao = temp.Theta() * 180/M_PI;

				Int_t particlePDG = fabs(mcParticle->getPDG());

				if(particlePDG == 22)	// fotoni (PDG == 22)
				{
					//PtPhotons.push_back(temp.Pt());
					//sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());
					histoPhotonPt.Fill(temp.Pt());

					Double_t tangens = temp.Pt()/temp.E();
					tng.push_back(tangens);
					sort(tng.begin(), tng.end(), greater<int>());

					if (ugao > 90 )
						{
							ugao = 180- ugao;
						}
					angle.push_back(ugao);
					sort(angle.begin(), angle.end(), greater<int>());



					if(temp.Pt() > minPt && temp.Theta()*180/M_PI > minTheta && temp.Theta()*180/M_PI < maxTheta)
					{
						photonsPassed++;
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

				histoPtMaxZoom.Fill(PtPhotons[0]);

				histoPhotonTg.Fill(tng[1]);
				histoPhotonTheta.Fill(angle[1]);




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

								histo2ndHighestThetaPhotons.Fill(thetaPhotons[1]);
								histoPt2Max.Fill(PtPhotons[1]);
								histoPt2MaxZoom.Fill(PtPhotons[1]);
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
	results.open("pTAngleCuts.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt << ", " << minTheta<<" < Theta <  "<< maxTheta  <<", "<< minInvMass<<" < M < " << maxInvMass <<  ", efikasnost = "<< 100.0 * (double)totalCutEvents/totalEvents  <<endl;
	results.close();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	histoPtMaxZoom.Write();
	histoPt2MaxZoom.Write();
	histoPtMax.Write();
	histoPt2Max.Write();
	histoPhotonPt.Write();
	histoPhotonTg.Write();
	histoPhotonTheta.Write();
	histo2ndHighestThetaPhotons.Write();


	rootFile.Close();


	TCanvas c1;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c1.SetCanvasSize(1000,650);
	c1.Divide(1,1,0.01,0.01);
	c1.cd(1);
//	c1.GetPad(1)->SetLogy();
	histoPtMax.Draw();
	c1.Print("pTAngleCutshistoMaxPt.eps");



	TCanvas c2;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c2.SetCanvasSize(1000,650);
	c2.Divide(1,1,0.01,0.01);
	c2.cd(1);
	//c2.GetPad(1)->SetLogy();
	histoPt2Max.Draw();
	c2.Print("pTAngleCutshisto2ndMax.eps");

	TCanvas c3;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c3.SetCanvasSize(1000,650);
	c3.Divide(1,1,0.01,0.01);
	c3.cd(1);
	c3.GetPad(1)->SetLogy();
	histoPtMaxZoom.Draw();
	c3.Print("pTAngleCutshistoMaxPtZoom.eps");



	TCanvas c4;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c4.SetCanvasSize(1000,650);
	c4.Divide(1,1,0.01,0.01);
	c4.cd(1);
	//c4.GetPad(1)->SetLogy();
	histoPt2MaxZoom.Draw();
	c4.Print("pTAngleCutshisto2ndMaxZoom.eps");

	TCanvas c5;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c5.SetCanvasSize(1000,650);
	c5.Divide(1,1,0.01,0.01);
	c5.cd(1);
	c5.GetPad(1)->SetLogy();
	histoPhotonPt.Draw();
	c5.Print("pTAngleCutsPhotonPt.eps");

	TCanvas c6;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c6.SetCanvasSize(1000,650);
	c6.Divide(1,1,0.01,0.01);
	c6.cd(1);
	c6.GetPad(1)->SetLogy();
	histoPhotonTheta.Draw();
	c6.Print("pTAngleCutsPhotonTheta.eps");

	TCanvas c7;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c7.SetCanvasSize(1000,650);
	c7.Divide(1,1,0.01,0.01);
	c7.cd(1);
	c7.GetPad(1)->SetLogy();
	histoPhotonTg.Draw();
	c7.Print("pTAngleCutsTangens.eps");

	TCanvas c8;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c8.SetCanvasSize(1000,650);
	c8.Divide(1,1,0.01,0.01);
	c8.cd(1);
	//c8.GetPad(1)->SetLogy();
	histo2ndHighestThetaPhotons.Draw();
	c8.Print("pTAngleCutsThetaStar.eps");

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

	TString rfName = "histoPtAngleCuts.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
