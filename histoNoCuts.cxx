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

	//gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/test/test/CLICdpStyle.C");

	TH1F histoPtMax ("HighestPhotonPt", "Max Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPt2Max ("2ndHighestPhotonPt", "2nd Max Photon Pt; Pt (GeV)", 1500, 0, 1500);
	TH1F histoPtMaxZoom ("Max Photon Pt Zoomed", "Max Photon Pt Zoomed; Pt (GeV)", 50, 0, 50);
	TH1F histoPt2MaxZoom ("2nd Max Photon Pt Zoomed", "2nd Max Photon Pt Zoomed; Pt (GeV)", 150, 0, 150);
	TH1F ptOfPhotons ("ptOfPhotons", "ptOfPhotons; Pt (GeV)", 750, 0, 1500);
	TH1F histoPhotonTheta ("Photon Theta", "Photon Theta", 90, 0, 90);
	TH1F histoPhotonTg ("tangens", "tangens", 10000, 0, 0.1);
	TH1F histo2ndHighestThetaPhotons ("2ndHighestThetaPhotons", "2nd Highest Theta of Photons; theta* (deg)", 90, 0 , 90);
	TH1F histoPhotonEnergy ("histoPhotonEnergy", "PhotonEnergy; E_{#gamma} (GeV)", 150, 0, 150);
	TH1F histoLeptonEnergy ("histoLeptonEnergy", "histoLeptonEnergy; E_{l} (GeV)", 110, 0, 550);
	TH1F histoLeptonTheta ("histoLeptonTheta", "histoLeptonTheta; #theta_{l} (GeV)", 180, 0, 180);



	//gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/test/test/CLICdpStyle.C");
	Float_t thetaE = 0;
	Float_t thetaP = 0;
	Float_t elecE = 0;
	Float_t posE = 0;
	int br=0;
	int brCestica = 0;
	int brElec= 0;
	int brPos = 0;


	TTree eventList("treeE", "ILD event list");
	eventList.Branch("thetaE", &thetaE, "thetaE");
	eventList.Branch("thetaP", &thetaP, "thetaP");
	eventList.Branch("elecE", &elecE, "elecE");
	eventList.Branch("posE", &posE, "posE");


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

		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja
		cout << "Broj dogadjaja po fajlu je  : " << brojDogadjaja << endl;




		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while((evt = lcReader->readNextEvent()) != 0)
		{
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");

    //       bool leptonFound = false;
            totalEvents++;

			vector<TLorentzVector> photons;
			vector<TLorentzVector> electrons;

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

				Float_t ugao = temp.Theta() * 180/M_PI;

				Int_t particlePDG = mcParticle->getPDG();

				if (particlePDG == 11) brElec++;
				if (particlePDG == -11) brPos++;
				if (fabs(particlePDG)==11) electrons.push_back(temp);


			/*	if (particlePDG ==11 && mcParticle->getGeneratorStatus() ==3){

				//	cout<< "našao početni elektron"<<endl;
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

					for (int l = 0; l < (int) daughter.size(); l++) // Petlja po potomcima Higsa koji su naš signal
								{
									int PDG = daughter[l] -> getPDG();

									TLorentzVector cerkaTemp (TVector3 (daughter[l] -> getMomentum()), daughter[l] -> getEnergy());
								//	cout << "e: "<< cerkaTemp.E();

									// Ako je potomak elektron 6
									if (PDG == 11) // promenljive za elektron 6
									{

							thetaE = cerkaTemp.Theta()*180./M_PI;
							elecE = cerkaTemp.E();
							histoLeptonTheta.Fill(cerkaTemp.Theta()*180/M_PI);
							histoLeptonEnergy.Fill(cerkaTemp.E());


									}

									// Ako je potomak pozitron 7
									if (PDG == -11) // promenljive za pozitron 7
									{

							thetaP = cerkaTemp.Theta()*180/M_PI;
							posE = cerkaTemp.E();
							histoLeptonTheta.Fill(cerkaTemp.Theta()*180/M_PI);
							histoLeptonEnergy.Fill(cerkaTemp.E());


									}



								}

				}*/



			/*	if (fabs(particlePDG == 11)) {
				//	histoLeptonEnergy.Fill(temp.E());
			//		histoLeptonTheta.Fill(temp.Theta()*180/M_PI);
				}

				if (particlePDG == -11){
					posE = temp.E();
					thetaP = temp.Theta() * 180/M_PI;
					histoLeptonTheta.Fill(ugao);
					histoLeptonEnergy.Fill(temp.E());
				}
				if (particlePDG ==11){
					elecE = temp.E();
					thetaE = temp.Theta() * 180 / M_PI;
					histoLeptonTheta.Fill(ugao);
					histoLeptonEnergy.Fill(temp.E());

				}*/



			}  // end of particle loop

			if (electrons.size()!=2){
				br++;
				cout << totalEvents<<endl;
			}

			eventList.Fill();
		} // End of event loop



		lcReader->close();



	} // End of file loop
	cout << "brojac:     " << br << endl;
	cout << "broj finalnih elektrona:     " << brElec << endl;
	cout << "broj finalnih posi:     " << brPos << endl;
	cout << "broj cestica:     " << brCestica << endl;
	cout << "Ukupan broj dogadjaja:     " << totalEvents << endl;
  //  cout << "Post-cut events in file:  " << totalCutEvents << endl;
    cout << "Percentage:               " << 100.0 * (double)totalCutEvents/totalEvents << endl;

/*	ofstream results ;
	results.open("histoNoCuts.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt << ", " << minTheta<<" < Theta <  "<< maxTheta  <<", "<< minInvMass<<" < M < " << maxInvMass <<  ", efikasnost = "<< 100.0 * (double)totalCutEvents/totalEvents  <<endl;
	results.close();*/

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	histoPtMaxZoom.Write();
	histoPt2MaxZoom.Write();
	histoPtMax.Write();
	histoPt2Max.Write();
	ptOfPhotons.Write();
	histoPhotonTg.Write();
	histoPhotonTheta.Write();
	histo2ndHighestThetaPhotons.Write();
	histoLeptonEnergy.Write();
	eventList.Write();
	histoLeptonTheta.Write();

	rootFile.Close();


	/*TCanvas c1;
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

	TCanvas c3;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c3.SetCanvasSize(1000,650);
	c3.Divide(1,1,0.01,0.01);
	c3.cd(1);
	c3.GetPad(1)->SetLogy();
	histoPtMaxZoom.Draw();
	c3.Print("noCutshistoMaxPtZoom.eps");



	TCanvas c4;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c4.SetCanvasSize(1000,650);
	c4.Divide(1,1,0.01,0.01);
	c4.cd(1);
	//c4.GetPad(1)->SetLogy();
	histoPt2MaxZoom.Draw();
	c4.Print("noCutshisto2ndMaxZoom.eps");

	TCanvas c5;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c5.SetCanvasSize(1000,650);
	c5.Divide(1,1,0.01,0.01);
	c5.cd(1);
	c5.GetPad(1)->SetLogy();
	ptOfPhotons.Draw();
	c5.Print("noCutsPhotonPt.eps");

	TCanvas c6;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c6.SetCanvasSize(1000,650);
	c6.Divide(1,1,0.01,0.01);
	c6.cd(1);
	c6.GetPad(1)->SetLogy();
	histoPhotonTheta.Draw();
	c6.Print("noCutsPhotonTheta.eps");

	TCanvas c7;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c7.SetCanvasSize(1000,650);
	c7.Divide(1,1,0.01,0.01);
	c7.cd(1);
	c7.GetPad(1)->SetLogy();
	histoPhotonTg.Draw();
	c7.Print("noCutsTangens.eps");

	TCanvas c8;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(1);
	c8.SetCanvasSize(1000,650);
	c8.Divide(1,1,0.01,0.01);
	c8.cd(1);
	//c8.GetPad(1)->SetLogy();
	histo2ndHighestThetaPhotons.Draw();
	c8.Print("noCutsThetaStar.eps");

	TCanvas c9;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c9.SetCanvasSize(1000,650);
	c9.Divide(1,1,0.01,0.01);
	c9.cd(1);
	//c4.GetPad(1)->SetLogy();
	histoPhotonEnergy.Draw();
	c9.Print("histoPhotonEnergy.pdf");
	return 0;*/



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
