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
#include <EVENT/LCRelation.h>
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

const Double_t minPt1 = 10;			   //minimalna energija fotona
const Double_t minPt2 = 10;			   //minimalna energija fotona
const Double_t minPt3 = 15;			   //minimalna energija fotona
const Double_t minPt4 = 20;			   //minimalna energija fotona

const double mH = 126.0;               //Higgs mass
const Double_t coneAngle = 2.5;		   //ugao konusa
const Double_t EnergyCenterMass = 3000; // energija u sistemu centra masa, koju koristim za missing energy
const Double_t minInvMass = 0;
const Double_t maxInvMass = 1500000;
const Double_t cutEremaining = 400;
const Double_t cutHiggsPt = 20;
const Double_t cutEhiggsMin = 100;
const Double_t cutEhiggsMax = 1000;
const Double_t cutCosinusHelicity = 0.9;


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


	gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/crtanje/histogrami/CLICdpStyle.C");

	//histogrami za različite kinematičke varijable
	TH1F histoTheta ("histoTheta", "Theta; #theta_{#gamma}", 50, 0, 180);
	TH1F histoPhi ("histoPhi", "Phi; #phi_{#gamma}", 50, 0, 3.14);
	TH1F histoPhotonEnergy ("histoPhotonEnergy", "PhotonEnergy; E_{#gamma} (GeV)", 150, 0, 1500);
	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{#gamma#gamma} (GeV)", 30, 110, 140);
	TH1F histoCandidatePt ("CandidatePt", "CandidatePt; Pt_{#gamma#gamma} (GeV)", 150, 0, 1550);
	TH1F histoPhotonPt ("PhotonPt", "PhotonPt; Pt_{#gamma} (GeV)", 50, 0, 50);
	TH1F histoRemainingEnergy ("Remaining_Visible_Energy", "Remaining_Visible_Energy; E_{remaining visible} (GeV)", 150, 0, 3000);
	TH1F histoCandidateEnergy ("CandidateE", "CandidateE; E_{#gamma#gamma} (GeV)",150, 0, 1550);
	TH1F histoCanditateTheta ("CandidateTheta", "CandidateTheta; #theta_{H}", 90, 0, 180);
	TH1F histoCandidatePhi ("CandidatePhi", "CandidatePhi' #phi_{#gamma}", 50, 0, 3.14);
	TH1F histoBoost ("CandidateBoost", "CandidateBoost; #beta_{#gamma#gamma}", 180, 0, 1);
	TH1F histoZbirPt ("ZbirPt", "ZbirPt; Pt_{1} + Pt_{2}", 300, 0, 3000);
	TH1F histoHigherEnergyPhoton ("HigherEnergyPhoton", "HigherEnergyPhoton; E_{#gamma1}", 150, 0, 1500);
	TH1F histoLowerEnergyPhoton ("LowerEnergyPhoton", "LowerEnergyPhoton; E_{#gamma2}", 150, 0, 1500);
	TH1F histoCosinusHelicityAngle ("cosinus_helicity_angle", "cosinus_helicity_Angle; cos_{#theta}", 50, 0, 1);
	TH1F histoHelicityAngle ("helicity_angle", "helicity_angle; #theta", 90, 0, 3.15);
	TH2F histogram ("test","title", 20, 0, 180, 20, 0, 180);
	TH1F histoangleBetweenPhotons ("angleBetweenPhotons", "angleBetweenPhotons; #alpha", 50, 0, 3.14);
	TH1F histoNumberPhotonsbyEvent ("Number of Photons by Event", "Number of Photons by Event", 10, 0, 10);
	TH1F anglePhotonParticle ("Angle Between photon and particle", "angle photon particle", 50, 0, 20);
	TH1F histoconeEnergy ("Cone Energy", "Cone Energy; E{cone} (GeV)", 100, 0, 100);
	TH1F anglePhotonParticlewE ("AngleBetweenphotonandparticle wE", " photonParticleAngleWE", 50, 0, 20);//otezinjeno sa Energijom
	TH1F histoconeEnergyFilter ("ConeEnergy", "ConeEnergy; E_{cone} (GeV)", 100, 0, 10);
	TH1F histoHardPhotonsByEvent ("HardPhotonsbyEvent", "HardPhotonsbyEvent", 10, 0, 10);
	TH1F histoPtOtherParticles ("PtOtherParticles", "PtOtherParticles; Pt (GeV)", 100, 0, 20);
	TH1F histoMissingEnergy ("Missing Energy", "MissingEnergy; E_{miss} (GeV)", 100, 0, 3000);
	TH1F histoTestPt ("TestPt", "TestPt; Pt (GeV)", 50, 0, 50);
	TH1F histoTheta1Photon ("histoTheta", "ThetaofcandidatePhoton; #theta_{#gamma}", 50, 0, 3.14);
	TH1F histoTheta2Photon ("histoTheta2", "Theta; #theta_{#gamma}", 50, 0, 3.14);
    TH1F histoPhoton1CandidateTheta ("Thetaof1stPhotonofCandidate", "Thetaof1stPhotonofCandidate; #theta_{#gamma}", 180, 0, 180 );
    TH1F histoPhoton2CandidateTheta ("Thetaof2ndPhotonofCandidate", "Thetaof2ndPhotonofCandidate; #theta_{#gamma}", 180, 0, 180 );
	TH1F histoPtof1stPhotonofCandidate ("Ptof1stPhotonofCandidate", "Ptof1stPhotonofCandidate; Pt (GeV)", 250, 0, 500);
	TH1F histoPtof2ndPhotonofCandidate ("Ptof2ndPhotonofCandidate", "Pt of 2nd Photon of Candidate; Pt (GeV)", 250, 0, 500);
	TH1F histoHighestPhotonPt ("HighestPhotonPt", "HighestPhotonPt; Pt (GeV)", 250, 0, 500);//najveći pT od fotona koju prodju Pt cut
	TH1F histo2ndHighestPhotonPt ("2ndHighestPhotonPt", "2ndHighestPhotonPt; Pt (GeV)", 250, 0, 500);//drugi najveći Pt od fotona koji prodju pt cut
	TH1F histo2ndHighestPhotonPtZoomed ("2ndHighestPhotonPtZoomed", "2ndHighestPhotonPtZoomed; Pt (GeV)", 50, 0, 50);//drugi najveći Pt od fotona koji prodju pt cut
	TH1F histoVisibleEnergy ("Visible_Energy", "Visible_Energy; E_{vis}(GeV)", 150, 0, 3000);
	TH1F histoZbirImpulsa ("zbir_impulsa", "zbir_impulsa; Pt_{hel}(GeV)", 100, 0, 1);
	TH1F histoCandidateMFineBinning ("histoCandidateMFineBinning", "histoCandidateMFineBinning; M_{H} (GeV)", 50, 125, 135);
	TH1F histoEnergyRatio ("Energy_Ratio", "Energy_Ratio; Ratio of reconstructed and generated Energy", 40, 0, 4);
	TH2F Energies ("Energies","Energies", 150, 0, 1500, 150, 0, 1500);


	TTree eventList("events", "ILD event list");
	varListGoran vl; /* SL specific */
	eventList.Branch("photonEnergyRatio", &(vl.photonEnergyRatio), "photonEnergyRatio");


	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	Int_t totalCutEvents = 0;
	Int_t totalEvents = 0;
	Int_t counterPhotonsPt1 = 0;
	Int_t counterEventsPt2 = 0;
	Int_t counterEventsPt3 = 0;
	Int_t counterEventsPt4 = 0;

	Int_t counterEremaining = 0;
	Int_t counterPreselectionCuts= 0;
	Int_t counterCosinusHelicity =0;
	Int_t counterEhiggs = 0;
	Int_t counterHiggsPt = 0;

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


		Int_t counterAllCuts = 0;
		//Int_t hardPhotons = 0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
		//	Int_t hardPhotonsByEvent = 0;
			brojac++;
			vector<TLorentzVector> reconstructedPhotons;
			vector<TLorentzVector> particles;   // other than photons
			vector <Double_t> PtPhotonsRec;
			vector <Double_t> PtPhotonsGen;

			vector <Double_t> ePhotonsRec;
			vector <Double_t> ePhotonsGen;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("PandoraPFANewPFOs");/*("PandoraPFOCollection");*/
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");

			//double_t Evis = 0;	//Ukupna energija po dogadjaju
			Double_t ERecTotal = 0;
			Double_t EGenTotal = 0;


			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t j = 0; j < recParticles->getNumberOfElements() ; j++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(j);

				TLorentzVector temprec; //četvorovektor u koji sakupljamo informacije o svakoj čestici


				const double *prec = recParticle->getMomentum(); // impuls čestice
				double erec = recParticle->getEnergy();	//energija čestice
				temprec.SetPxPyPzE(prec[0], prec[1], prec[2], erec);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDGRec = fabs(recParticle->getType());
			//	cout<<"particle PDG is: "<<particlePDG<<endl;

				if(particlePDGRec == 22)	//rad sa fotonima (PDG=22)
				{
					//Double_t theta = temp.Theta(); //promenljiva koja nam daje Theta čestice
					//Double_t phi = temp.Phi();		//promenljiva koja nam daje Phi
					Double_t PtRec = temprec.Pt();		//promenljiva koja nam daje Pt čestice
					Double_t ERec = temprec.E();		//promenljiva koja nam daje E čestice
					ERecTotal +=ERec;
					ePhotonsRec.push_back(ERec);


					histoPhotonPt.Fill(PtRec);

					//uzimamo u obzir samo one fotone koji nam prodju uslove
					if(temprec.Pt() > minPt1)
					{
						PtPhotonsRec.push_back(temprec.Pt());
						sort(PtPhotonsRec.begin(), PtPhotonsRec.end(), greater<int>());


						counterPhotonsPt1++;

						allPhotons.push_back(temprec);//sakupili smo sve fotone
						reconstructedPhotons.push_back(temprec);
					}

				}




			}  // end of particle loop


			for (Int_t i = 0; i < mcParticles->getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);


				TLorentzVector tempgen;

				const double *pgen = mcParticle->getMomentum();
				double egen = mcParticle->getEnergy();
				tempgen.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen);
				if (mcParticle->getGeneratorStatus() != 1) continue;


				Int_t particlePDGGen = fabs(mcParticle->getPDG());

				if(particlePDGGen == 22)	// fotoni (PDG == 22)
				{
					PtPhotonsGen.push_back(tempgen.Pt());
					sort(PtPhotonsGen.begin(), PtPhotonsGen.end(), greater<int>());
					histoPhotonPt.Fill(tempgen.Pt());
					Double_t PtGen = tempgen.Pt();		//promenljiva koja nam daje Pt čestice
					Double_t EGen = tempgen.E();		//promenljiva koja nam daje E čestice
					EGenTotal+=EGen;
					ePhotonsGen.push_back(EGen);


					if(tempgen.Pt() > minPt1)
					{
						PtPhotonsGen.push_back(tempgen.Pt());
						sort(PtPhotonsGen.begin(), PtPhotonsGen.end(), greater<int>());

					}

				}
			} //End of generated particle loop


			histoEnergyRatio.Fill(ERecTotal/EGenTotal);
			//vl.photonEnergyRatio = ERecTotal/EGenTotal;

		//	cout<<"broj rekonstruisanih fotona je: "<<ePhotonsRec.size()<<endl;
			//cout<<"broj generisanih fotona je: "<<ePhotonsGen.size()<<endl;


            EVENT::LCCollection* links = evt->getCollection("RecoMCTruthLink");
            for (int i = 0; i < links->getNumberOfElements(); i++){
                EVENT::LCRelation* link = (EVENT::LCRelation*) links->getElementAt(i);
                IMPL::ReconstructedParticleImpl* reconstructed = (IMPL::ReconstructedParticleImpl*) link->getFrom();
                EVENT::MCParticle* generated = (EVENT::MCParticle*) link->getTo();

				double recoEnergy = reconstructed->getEnergy();
				if (generated !=NULL)
				{
					double genEnergy = generated->getEnergy();
					Energies.Fill(recoEnergy, genEnergy);
					vl.photonEnergyRatio =recoEnergy/genEnergy;

				}
			//	double genEnergy = generated->getEnergy();

			  //Energies.Fill(recoEnergy, 1000);
            }

            eventList.Fill();
		} // End of event loop


		totalCutEvents+=numCutEvents;//ukupan broj dogadjaja koji prodje cutove iz svih fajlova
		totalEvents+=brojDogadjaja; //ukupan broj dogadjaja iz svih fajlova

		cout << "Broj event-a posle cut-a je: " << numCutEvents << endl;
		cout << "Broj fotona posle primarnog cut-a  je: " << counterPhotonsPt1 << endl;
		cout << "Broj fotona posle svih cut-ova  je: " << counterAllCuts << endl;
//		cout << "broj ogadjaja posle cosinus cuta  "<< counterCosinusHelicity<<endl;


		lcReader->close();

	} // End of file loop

	Double_t percentageEremaining = (Double_t) counterEremaining/totalEvents * 100;
	Double_t percentageCosinusHelicity = (Double_t) counterCosinusHelicity/totalEvents * 100;
	Double_t percentageEhiggs = (Double_t) counterEhiggs/totalEvents * 100;
	Double_t percentageHiggsPt = (Double_t) counterHiggsPt/totalEvents * 100;
	Double_t percentagePreselectionCuts = (Double_t) counterPreselectionCuts/totalEvents * 100;


/*	Double_t percentagePt1 = (Double_t) totalCutEvents/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt1 << " %"<< endl;

	Double_t percentagePt2 = (Double_t) counterEventsPt2/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt2 << " %"<< endl;

	Double_t percentagePt3 = (Double_t) counterEventsPt3/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt3 << " %"<< endl;

	Double_t percentagePt4 = (Double_t) counterEventsPt4/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt4 << " %"<< endl;*/

	cout <<"Ukupan broj dogadjaja koji prodju katove je: " << totalCutEvents << endl;
	cout <<"Ukupan broj dogadjaja je: " << totalEvents << endl;


/*	ofstream file ;
	file.open("efikasnostPt.txt", ios_base::out | ios_base::app);

	file << minPt1 << "\t" << percentagePt1 << endl;
	file << minPt2 << "\t" << percentagePt2 << endl;
	file << minPt3 << "\t" << percentagePt3 << endl;
	file << minPt4 << "\t" << percentagePt4 << endl;
	file << "__________________________"<<endl;
	file.close();*/

/*	ofstream results ;
	results.open("rezultati.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt1 <<", "<< minInvMass<<" < M < " << maxInvMass <<  "  ,efikasnost = "<<percentagePt1 <<endl;
	results.close();*/
/*

	ofstream events;
	events.open( "events.txt " ,ios_base::out | ios_base::app);
	events << "ukupan broj događaja:  "<< totalEvents<< ",  broj događaja posle cut-va:  "<< totalCutEvents<<endl;
	events.close();


	ofstream preselectionEfficiency;
	preselectionEfficiency.open( "preselection.txt " ,ios_base::out | ios_base::app);
	preselectionEfficiency << "preselection efficiency: " << percentagePreselectionCuts <<endl;
	preselectionEfficiency << "efficiency ERemaining: " << percentageEremaining <<endl;
	preselectionEfficiency << "efficiency cosinus helicity: " << percentageCosinusHelicity <<endl;
	preselectionEfficiency << "efficiency Higgs Energy: " << percentageEhiggs <<endl;
	preselectionEfficiency << "efficiency Higgs Pt: " << percentageHiggsPt <<endl;
	preselectionEfficiency << "________________________________________________"<<endl;
	preselectionEfficiency.close();
*/


	/*TGraph efikasnostproba ("proba.txt", "%lg %lg", "\t");
	efikasnostproba.GetXaxis()->SetTitle("Pt_{#gamma} (GeV)");
	efikasnostproba.GetYaxis()->SetTitle("Efficiency");
	efikasnostproba.SetTitle("Signal");*/



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	eventList.Write();
    Energies.Write();



	rootFile.Close();

/*	TCanvas c1;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c1.SetCanvasSize(1000,650);
	c1.Divide(1,1,0.01,0.01);
	c1.cd(1);
	histoEnergyRatio.Draw();
	c1.Print("histoEnergyRatio.eps");

	TCanvas c2;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c2.SetCanvasSize(1000,650);
	c2.Divide(1,1,0.01,0.01);
	c2.cd(1);
	Energies.Draw();
	c2.Print("histoEnergije.eps");*/


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

	TString rfName = "reco_gen_photons.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
