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
const Double_t maxConeEnergy = 20000 ;	   //maksimalna energija konusa
const Double_t leptonPt = 20000;		   //energija leptona posle koje zadovoljavaju uslov leptonFound
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
	TH1F histoTheta ("histoTheta", "Theta; #theta_{#gamma}", 90, 0, 180);
	TH1F histoPhi ("histoPhi", "Phi; #phi_{#gamma}", 50, 0, 3.14);
	TH1F energyOfPhotons ("energyOfPhotons", "PhotonEnergy; E_{#gamma} (GeV)", 150, 0, 1500);
	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{#gamma#gamma} (GeV)", 30, 110, 140);
	TH1F histoCandidatePt ("CandidatePt", "CandidatePt; Pt_{#gamma#gamma} (GeV)", 150, 0, 1550);
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
	TH1F histoNumberPhotonsbyEvent ("Number of Photons by Event", "Number of Photons by Event", 100, 0, 100);
	TH1F anglePhotonParticle ("Angle Between photon and particle", "angle photon particle", 50, 0, 20);
	TH1F histoconeEnergy ("ConeEnergy", "ConeEnergy; E_{cone} (GeV)", 100, 0, 100);
	TH1F anglePhotonParticlewE ("AngleBetweenphotonandparticle wE", " photonParticleAngleWE", 50, 0, 20);//otezinjeno sa Energijom
	TH1F histoconeEnergyFilter ("ConeEnergyFilter", "ConeEnergyFilter; E_{cone} (GeV)", 100, 0, 20);
	TH1F histoHardPhotonsByEvent ("HardPhotonsbyEvent", "HardPhotonsbyEvent", 10, 0, 10);
	TH1F histoPtOtherParticles ("PtOtherParticles", "PtOtherParticles; Pt (GeV)", 100, 0, 20);
	TH1F histoMissingEnergy ("MissingEnergy", "MissingEnergy; E_{miss} (GeV)", 100, 0, 3000);
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
	TH1F ptOfPhotons ("ptOfPhotons", "ptOfPhotons; Pt_{#gamma} (GeV)", 250, 0, 500);
	TH1F thetaOfPhotons ("thetaOfPhotons", "thetaOfPhotons; #theta_{#gamma}", 90, 0, 180);
	TH1F thetaOfHiggsPhotons ("thetaOfHiggsPhotons", "thetaOfHiggsPhotons; #theta_{#gamma}", 90, 0, 180);
	//TH1F energyOfPhotons ("energyOfPhotons", "energyOfPhotons; E_{#gamma} (GeV)", 150, 0, 1500);
	TH1F pLOfPhotons ("pLOfPhotons", "pLOfPhotons; Pz_{#gamma} (GeV)", 250, 0, 500);
	//TH1F thetaOfPhotons ("thetaOfPhotons", "thetaOfPhotons; #theta_{#gamma}", 90, 0, 180);
	TH1F PtofCandidatePhotons ("PtofCandidatePhotons", "PtofCandidatePhotons; Pt_{#gamma} (GeV)", 100, 0, 1500);
	TH1F energyofCandidatePhotons ("energyofCandidatePhotons", "energyofCandidatePhotons; E_{#gamma} (GeV)", 100, 0, 1500);
	TH1F thetaOfCandidatePhotons ("thetaOfCandidatePhotons", "thetaOfCandidatePhotons; #theta_{#gamma}", 90, 0, 180);
	TH1F energyofCandidatePhotonsZoomed ("energyofCandidatePhotonsZoomed", "energyofCandidatePhotonsZoomed; E_{#gamma} (GeV)", 6, 0, 30);
	TH1F PtofCandidatePhotonsZoomed ("PtofCandidatePhotonsZoomed", "PtofCandidatePhotonsZoomed; Pt_{#gamma} (GeV)", 6, 0, 30);
	TH1F histoconeEnergyfor2Photons ("ConeEnergy2Photons", "ConeEnergy; E_{cone} (GeV)", 100, 0, 100);

	TH1F histoconeEnergy15 ("ConeEnergy15", "ConeEnergy15; E_{cone} (GeV)", 100, 0, 100);
	TH1F histoconeEnergy1 ("ConeEnergy1", "ConeEnergy1; E_{cone} (GeV)", 100, 0, 100);
	TH1F histoconeEnergy10 ("ConeEnergy10", "ConeEnergy10; E_{cone} (GeV)", 100, 0, 100);
	TH1F histoconeEnergy5 ("ConeEnergy5", "ConeEnergy5; E_{cone} (GeV)", 100, 0, 100);
	TH1F histoconeEnergy2 ("ConeEnergy2", "ConeEnergy5; E_{cone} (GeV)", 100, 0, 100);
	TH1F histoconeEnergy3 ("ConeEnergy3", "ConeEnergy5; E_{cone} (GeV)", 100, 0, 100);



	TTree eventList("eventsSignal", "ILD event list");
	varListGoran vl; /* SL specific */
	eventList.Branch("invM", &(vl.CandidateInvariantM), "invM");
	/*eventList.Branch("CanE", &(vl.CandidateEnergy), "CanE");
	eventList.Branch("CanTheta", &(vl.CandidateTheta), "CanTHeta");
	eventList.Branch("CanPt", &(vl.CandidatePt), "CanPt");
	eventList.Branch("CanPhi", &(vl.CandidatePhi), "CanPhi");
	eventList.Branch("photonsPerEvent", &(vl.NumberPhotonsbyEvent), "photonsPerEvent");
	eventList.Branch("cos_hel_angle", &(vl.cosHelAngle), "cos_hel_angle");
	eventList.Branch("energyPhotons", &(vl.energyOfAllPhotons), "energyPhotons");
	eventList.Branch("thetaPhotons", &(vl.thetaOfAllPhotons), "thetaPhotons");
	eventList.Branch("ptPhotons", &(vl.ptOfAllPhotons), "ptPhotons");
	eventList.Branch("highestPt", &(vl.HighestPhotonPt), "highestPt");
	eventList.Branch("2ndHighestPt", &(vl.secondHighestPhotonPt), "2ndHighestPt");
	eventList.Branch("MissingE", &(vl.MissingEnergy), "MissingE");
	eventList.Branch("ERem", &(vl.ERemaining), "ERem");
	eventList.Branch("Evis", &(vl.EVisible), "Evis");
	eventList.Branch("EcallE", &(vl.EcallEnergy), "EcallE");
	eventList.Branch("HcallEnergy", &(vl.HcallEnergy), "HcallEnergy");
	eventList.Branch("TotalcallEnergy", &(vl.TotalcallEnergy), "TotalcallEnergy");
	eventList.Branch("pLCandidate", &(vl.pLCandidate), "pLCandidate");
	eventList.Branch("pLCandidatePhotons1", &(vl.pLCandidatePhotons1), "pLCandidatePhotons1");
	eventList.Branch("pLCandidatePhotons2", &(vl.pLCandidatePhotons2), "pLCandidatePhotons2");
	eventList.Branch("thetaCandidatePhotons1", &(vl.thetaCandidatePhotons1), "thetaCandidatePhotons1");
	eventList.Branch("thetaCandidatePhotons2", &(vl.thetaCandidatePhotons2), "thetaCandidatePhotons2");
	eventList.Branch("pTCandidatePhotons1", &(vl.pTCandidatePhotons1), "pTCandidatePhotons1");
	eventList.Branch("pTCandidatePhotons2", &(vl.pTCandidatePhotons2), "pTCandidatePhotons2");*/
	eventList.Branch("coneEnergy", &(vl.coneEnergy), "coneEnergy");
	eventList.Branch("coneAngle", &(vl.coneAngle), "coneAngle");





	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	Int_t totalCutEvents = 0;
	Int_t totalEvents = 0;
	Int_t counterPhotonsPt1 = 0;
	Int_t counterEventsPt2 = 0;
	Int_t counterEventsPt3 = 0;
	Int_t counterEventsPt4 = 0;

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
			Int_t hardPhotonsByEvent = 0;
			brojac++;
			vector<TLorentzVector> photons;
			vector<TLorentzVector> particles;   // other than photons
			vector <Double_t> PtPhotons;
			vector <Double_t> ConeEnergyOfPhotons;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("SelectedPandoraPFANewPFOs_Reprocess");/*("PandoraPFANewPFOs");*/

			double_t Evis = 0;	//Ukupna energija po dogadjaju PandoraPFOsDefault_Reprocess, LooseSelectedPandoraPFANewPFOs_Reprocess,
			//SelectedPandoraPFANewPFOs_Reprocess, TightSelectedPandoraPFANewPFOs_Reprocess

			bool leptonFound = false;	//uslov da nemamo leptone
			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici


				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(recParticle->getType());
			//	cout<<"particle PDG is: "<<particlePDG<<endl;

				if(particlePDG == 22)	//rad sa fotonima (PDG=22)
				{
					Double_t theta = temp.Theta()*180/M_PI; //promenljiva koja nam daje Theta čestice
					//Double_t phi = temp.Phi();		//promenljiva koja nam daje Phi
				//	Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice
					ptOfPhotons.Fill(temp.Pt());
					energyOfPhotons.Fill(temp.E());
					thetaOfPhotons.Fill(theta);

					PtPhotons.push_back(temp.Pt());




					//uzimamo u obzir samo one fotone koji nam prodju uslove
					if(temp.Pt() >minPt1 /*minPt1*/)
					{
	histo2ndHighestPhotonPtZoomed.Fill(PtPhotons[1]);

						counterPhotonsPt1++;

						allPhotons.push_back(temp);//sakupili smo sve fotone
						TLorentzVector currentConeAxis = temp; //foton oko koga pravimo konus
						Double_t coneEnergy = 0; //Energija konusa oko fotona
						Double_t coneEnergy15 = 0; //Energija konusa oko fotona
						Double_t coneEnergy10 = 0; //Energija konusa oko fotona
						Double_t coneEnergy5 = 0; //Energija konusa oko fotona
						Double_t coneEnergy1 = 0; //Energija konusa oko fotona
						Double_t coneEnergy2 = 0; //Energija konusa oko fotona
						Double_t coneEnergy3 = 0; //Energija konusa oko fotona




						for (Int_t k = 0; k < recParticles->getNumberOfElements() ; k++) /*&& coneEnergy <= maxConeEnergy -zaustavlja fot loop kada predje maxE */
						{
							IMPL::ReconstructedParticleImpl* mcParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(k);
							const double *impuls = mcParticle->getMomentum();
							double energija = mcParticle->getEnergy();
							TLorentzVector otherParticle;
							otherParticle.SetPxPyPzE(impuls[0], impuls[1], impuls[2],energija);//četvorovektor drugih čestica
							if (currentConeAxis == otherParticle) continue; //da ne bi ubrojali foton oko kog pravimo konus
                            Double_t particleAngle = currentConeAxis.Angle(otherParticle.Vect())* 180 / TMath::Pi() ;//ugao izmedju fotona oko kog pravimo konus i čestice
							//anglePhotonParticle.Fill(particleAngle);//histogram uglova izmedju fotona i čestice
							//anglePhotonParticlewE.Fill(particleAngle, otherParticle.E());//histogram uglova izmedju fotona i čestice otežinjen energijom te čestice
                            vl.coneAngle = particleAngle;

							//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
							if (particleAngle <=2.5)
							{
								coneEnergy += otherParticle.E();
								vl.coneEnergy = coneEnergy;
							}



							if (particleAngle <=15)
							{
								coneEnergy15 += otherParticle.E();
							}


							if (particleAngle <=10)
							{
								coneEnergy10 += otherParticle.E();
							}

							if (particleAngle <=5)
							{
								coneEnergy5 += otherParticle.E();
							}

							if (particleAngle <=1)
							{
								coneEnergy1 += otherParticle.E();
							}


							if (particleAngle <=2)
							{
								coneEnergy2 += otherParticle.E();
							}



							if (particleAngle <=3)
							{
								coneEnergy3 += otherParticle.E();
							}








						}
						ConeEnergyOfPhotons.push_back(coneEnergy);
						histoconeEnergy.Fill(coneEnergy); //histogram koji iscrtava energiju konusa

						histoconeEnergy10.Fill(coneEnergy10); //histogram koji iscrtava energiju konusa
						histoconeEnergy15.Fill(coneEnergy15); //histogram koji iscrtava energiju konusa
						histoconeEnergy5.Fill(coneEnergy5); //histogram koji iscrtava energiju konusa
						histoconeEnergy1.Fill(coneEnergy1); //histogram koji iscrtava energiju konusa
						histoconeEnergy3.Fill(coneEnergy2); //histogram koji iscrtava energiju konusa
						histoconeEnergy2.Fill(coneEnergy3); //histogram koji iscrtava energiju konusa


						//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
						if (coneEnergy <= maxConeEnergy)
							{
								counterAllCuts++;//brojač za fotone koji prodju sve cutove
								photons.push_back(temp);
								histoconeEnergyFilter.Fill(coneEnergy);//histogram koji iscrtava energiju konusa

							}

					}

				}
				//dogadjaj ne prolazi ako se u dogadjaju detektuje lepton
				else if (particlePDG != 22)
					{

					if(temp.Pt() > leptonPt)
						leptonFound = true;
					histoPtOtherParticles.Fill(temp.Pt());
					}
				if( temp.Pt() > 5)
				{
					Evis += e ;//ukupna enrgija po dogadjaju
				}



			}  // end of particle loop

			sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());

			histoHighestPhotonPt.Fill(PtPhotons[0]);
			histo2ndHighestPhotonPt.Fill(PtPhotons[1]);

			if (ConeEnergyOfPhotons.size()==2)
						{
						histoconeEnergyfor2Photons.Fill(ConeEnergyOfPhotons[0]);
						histoconeEnergyfor2Photons.Fill(ConeEnergyOfPhotons[1]);
						}


			histoNumberPhotonsbyEvent.Fill(photons.size());
			Int_t Emiss = EnergyCenterMass - Evis; // missing energy
			histoMissingEnergy.Fill(Emiss);

			//cout << "lepton found: " <<leptonFound <<endl;
			if(leptonFound) continue; //ako je u finalnom stanju imamo lepton ili kvark (sem neutrina), događaj se preskače


			// At this point, all relevant photons have been collected.

		//	cout << "photon size is: " << photons.size()<<endl;
			//proveravamo da li u eventu ima detektovano više od jednog fotona koji bi bili kandidat za higsov bozon
			if (photons.size() == photonSize)// photonSize)
			{
				vector<CandidateData> candidates;

				//histoNumberPhotonsbyEvent.Fill(photons.size());//histogram koji nam pokazuje broj fotona po događaju



				for (int i = 0; i < (int)photons.size() - 1; ++i)
				{
					for (int j = i + 1; j < photons.size(); ++j)
					{
						CandidateData current(photons[i], photons[j]);//uzima dva fotona

						TLorentzVector pair = current.Higgs();//sabira četvorovektore

						if (minInvMass < pair.M() && pair.M() < maxInvMass)//uslov da bi par bio higs kandidat
						{
							candidates.push_back(current);

						}
					}
				}

				// Pt tests
				bool pt2Found = false;
				bool pt3Found = false;
				bool pt4Found = false;
				for (int i = 0; i < candidates.size(); ++i)
				{
					CandidateData candidate = candidates[i];
					TLorentzVector photon1 = candidate.Photon1;
					TLorentzVector photon2 = candidate.Photon2;

					if (!pt2Found && (photon1.Pt() > minPt2 && photon2.Pt() > minPt2))
					{
						counterEventsPt2++;
						pt2Found = true;
					}

					if (!pt3Found && (photon1.Pt() > minPt3 && photon2.Pt() > minPt3))
					{
						counterEventsPt3++;
						pt3Found = true;
					}

					if (!pt4Found && (photon1.Pt() > minPt4 && photon2.Pt() > minPt4))
					{
						counterEventsPt4++;
						pt4Found = true;
					}
				}

				//proverava koji kandidat ima najbližu vrednost higsovom bozonu
				vector<CandidateData>::iterator candidate = min_element(candidates.begin(), candidates.end(), IsBetterHiggsCandidate);
				if (candidate != candidates.end())

				{
					numCutEvents++;//brojimo dogadjaje koji prodju uslove za energiju, invarijantnu masu, uglove,Pt

					TLorentzVector higgs = candidate->Higgs();
					Double_t Ehiggs = higgs.Energy();//energija kandidata
					Double_t Eremaining = Evis - Ehiggs;//vidljiva energija posle higgsa

					vl.CandidateInvariantM = higgs.M();

					TLorentzVector photon = candidate->Photon1;
					TLorentzVector photon2 = candidate->Photon2;
				    Double_t theta1 = photon.Theta()*180/M_PI;
				    Double_t theta2 = photon2.Theta()*180/M_PI;

				    PtofCandidatePhotons.Fill(photon.Pt());
				    PtofCandidatePhotons.Fill(photon2.Pt());

				    energyofCandidatePhotons.Fill(photon.E());
				    energyofCandidatePhotons.Fill(photon2.E());


				    PtofCandidatePhotonsZoomed.Fill(photon.Pt());
				    PtofCandidatePhotonsZoomed.Fill(photon2.Pt());

				    energyofCandidatePhotonsZoomed.Fill(photon.E());
				    energyofCandidatePhotonsZoomed.Fill(photon2.E());


				    thetaOfCandidatePhotons.Fill(theta1);
				    thetaOfCandidatePhotons.Fill(theta2);



					histoPhoton1CandidateTheta.Fill(theta1);
					histoPhoton2CandidateTheta.Fill(theta2);

					Double_t Px11 = candidate->Photon1.X();
					Double_t Py11 = candidate->Photon1.Y();
					Double_t Px22 = candidate->Photon2.X();
					Double_t Py22 = candidate->Photon2.Y();
					Double_t higherPtPhoton = TMath::Max(TMath::Sqrt(Px11*Px11 +Py11*Py11), TMath::Sqrt(Px22*Px22 +Py22*Py22));//Energija energičnijeg fotona koji je kandidat
					Double_t lowerPtPhoton = TMath::Min(TMath::Sqrt(Px11*Px11 +Py11*Py11), TMath::Sqrt(Px22*Px22 +Py22*Py22));//Energija manje energičnog fotona koji je kandidat

					histoPtof1stPhotonofCandidate.Fill(higherPtPhoton);
					histoPtof2ndPhotonofCandidate.Fill(lowerPtPhoton);

					TVector3 boosttoparent = -(higgs.BoostVector());//prelazimo u sitem CM


					photon.Boost(boosttoparent);//prebacujemo fotone u sistem CM
					photon2.Boost(boosttoparent);
					TVector3 photon_3v=photon.Vect();
					TVector3 photon2_3v = photon2.Vect();
					TVector3 higgs_3v = higgs.Vect();

					Double_t anglePhoton1 = photon.Angle(higgs_3v);
					Double_t anglePhoton2 = photon2.Angle(higgs_3v);

					Double_t lowerAngle = TMath::Min(anglePhoton1, anglePhoton2) ;
					//Double_t higherAngle = TMath::Max(anglePhoton1, anglePhoton2) ;
					Double_t cosinusHelicity = TMath::Cos(lowerAngle);


					Double_t higherEnergyPhotonE = TMath::Max(candidate->Photon1.E(), candidate->Photon2.E());//Energija energičnijeg fotona koji je kandidat
					Double_t lowerEnergyPhotonE = TMath::Min(candidate->Photon1.E(), candidate->Photon2.E());//Energija manje energičnog fotona koji je kandidat

					Double_t Px1 = photon.X();
					Double_t Py1 = photon.Y();
					Double_t Pz1 = photon.Z();
					Double_t Px2 = photon2.X();
					Double_t Py2 = photon2.Y();
					Double_t Pz2 = photon2.Z();

					Double_t zbirPt = TMath::Sqrt(Px1*Px1 +Py1*Py1) + TMath::Sqrt(Px2*Px2 +Py2*Py2); // zbir transverzalnih impulsa
					Double_t Theta1 =candidate->Photon1.Theta()*180/M_PI;
					Double_t Theta2 = candidate->Photon2.Theta()*180/M_PI;
					Double_t Phi1 = candidate->Photon1.Phi();
					Double_t Phi2 = candidate->Photon2.Phi();
					Double_t CMTheta2 =photon2.Theta()*180/M_PI;
					Double_t CMTheta1 =photon.Theta()*180/M_PI;
					Double_t zbirImpulsa = TMath::Sqrt(Px1*Px1 +Py1*Py1 + Pz1*Pz1) - TMath::Sqrt(Px2*Px2 +Py2*Py2 + Pz2*Pz2);

				//	cout << "zbir impulsa je: "<<zbirImpulsa<<endl;// proveravamo da li su fotoni back to back


				//	Double_t angleBetweenPotons = candidate->Photon1.Angle(candidate->Photon2.Vect());

					//energyOfPhotons.Fill(higherEnergyPhotonE, lowerEnergyPhotonE);
					histoTheta.Fill(Theta1);
					histoTheta.Fill(Theta2);
					histoPhi.Fill(Phi1);
					histoPhi.Fill(Phi2);
					histoCandidateM.Fill(higgs.M());
					histoCandidateMFineBinning.Fill(higgs.M());
					histoCandidatePt.Fill(higgs.Pt());
					histoTestPt.Fill(higgs.Pt());
					histoCandidateEnergy.Fill(higgs.Energy());
					histoCanditateTheta.Fill(higgs.Theta()*180/M_PI);
					histoCandidatePhi.Fill(higgs.Phi());
					histoRemainingEnergy.Fill(Eremaining) ;
					histoBoost.Fill(higgs.BoostVector().Mag());
					histoHigherEnergyPhoton.Fill(higherEnergyPhotonE);
					histoLowerEnergyPhoton.Fill(lowerEnergyPhotonE);
					histoZbirPt.Fill(zbirPt);
					histogram.Fill(CMTheta1, CMTheta2);
					histoHardPhotonsByEvent.Fill(hardPhotonsByEvent);
					histoVisibleEnergy.Fill(Evis);
					histoCosinusHelicityAngle.Fill(cosinusHelicity);
					histoHelicityAngle.Fill(anglePhoton1);
					histoHelicityAngle.Fill(anglePhoton2);
					histoZbirImpulsa.Fill(zbirImpulsa);


				}
			}


		} // End of event loop


		totalCutEvents+=numCutEvents;//ukupan broj dogadjaja koji prodje cutove iz svih fajlova
		totalEvents+=brojDogadjaja; //ukupan broj dogadjaja iz svih fajlova

		cout << "Broj event-a posle cut-a je: " << numCutEvents << endl;
		cout << "Broj fotona posle primarnog cut-a  je: " << counterPhotonsPt1 << endl;
		cout << "Broj fotona posle svih cut-ova  je: " << counterAllCuts << endl;






		/*for (Int_t i = 0; i < allPhotons.size()-1 ; i++)
			{

		for (Int_t j=i+1; j <  allPhotons.size(); j++)
			{
				Double_t tempAngle = allPhotons[i].Angle(allPhotons[j].Vect());
				histoangleBetweenPhotons.Fill(tempAngle);

			}

			}*/



		lcReader->close();

	} // End of file loop

	Double_t percentagePt1 = (Double_t) totalCutEvents/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt1 << " %"<< endl;

	Double_t percentagePt2 = (Double_t) counterEventsPt2/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt2 << " %"<< endl;

	Double_t percentagePt3 = (Double_t) counterEventsPt3/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt3 << " %"<< endl;

	Double_t percentagePt4 = (Double_t) counterEventsPt4/totalEvents * 100;
	cout << "Efikasnost je: " << percentagePt4 << " %"<< endl;

	cout <<"Ukupan broj dogadjaja koji prodju katove je: " << totalCutEvents << endl;
	cout <<"Ukupan broj dogadjaja je: " << totalEvents << endl;


	ofstream file ;
	file.open("proba.txt", ios_base::out | ios_base::app);

	file << minPt1 << "\t" << percentagePt1 << endl;
	file << minPt2 << "\t" << percentagePt2 << endl;
	file << minPt3 << "\t" << percentagePt3 << endl;
	file << minPt4 << "\t" << percentagePt4 << endl;
	file.close();

	ofstream results ;
	results.open("rezultati.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt1 <<", "<< minInvMass<<" < M < " << maxInvMass << ", leptonPt > "<< leptonPt << ", cone Angle < "<< 2*coneAngle << ", efikasnost = "<<percentagePt1 <<endl;
	results.close();

	ofstream events;
	events.open( "events.txt " ,ios_base::out | ios_base::app);
	events << "ukupan broj događaja:  "<< totalEvents<< ",  broj događaja posle cut-va:  "<< totalCutEvents;
	events.close();

	TGraph efikasnostproba ("proba.txt", "%lg %lg", "\t");
	efikasnostproba.GetXaxis()->SetTitle("Pt_{#gamma} (GeV)");
	efikasnostproba.GetYaxis()->SetTitle("Efficiency");
	efikasnostproba.SetTitle("Signal");



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	//eventList.Write();

	histoconeEnergy.Write();
	histoconeEnergy15.Write();
	histoconeEnergy10.Write();
	histoconeEnergy5.Write();
	histoconeEnergy1.Write();
	histoconeEnergy2.Write();
	histoconeEnergy3.Write();


	histoCanditateTheta.Write();
	histoCandidatePt.Write();
	histoCandidateM.Write();
	histoCandidateEnergy.Write();
	histoCandidatePhi.Write();
	histoBoost.Write();
	histoRemainingEnergy.Write();
	histoHelicityAngle.Write();
	histoHigherEnergyPhoton.Write();
	histoLowerEnergyPhoton.Write();
	histoZbirPt.Write();
	energyOfPhotons.Write();
	//histoTheta.Write();
	//histoPhi.Write();
	//histoangleBetweenPhotons.Write();
	histoMissingEnergy.Write();
	//efikasnost.Write();
	histoPhoton1CandidateTheta.Write();
	histoPhoton2CandidateTheta.Write();
	histoPtof1stPhotonofCandidate.Write();
	histoPtof2ndPhotonofCandidate.Write();
	histoconeEnergyFilter.Write();
	histoHighestPhotonPt.Write();
	histo2ndHighestPhotonPt.Write();
	histoCosinusHelicityAngle.Write();
	histoNumberPhotonsbyEvent.Write();
	histoCandidateMFineBinning.Write();
	ptOfPhotons.Write();
	PtofCandidatePhotons.Write();
	energyofCandidatePhotons.Write();
	thetaOfCandidatePhotons.Write();
	energyofCandidatePhotonsZoomed.Write();
	PtofCandidatePhotonsZoomed.Write();
    histoconeEnergyfor2Photons.Write();
	eventList.Write();


	rootFile.Close();
/*
	TCanvas c1;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c1.SetCanvasSize(1000,650);
	c1.Divide(1,1,0.01,0.01);
	c1.cd(1);
	histoCandidateM.Draw();
	c1.Print("histoCandidateM.pdf");

	TCanvas c2;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c2.SetCanvasSize(1000,650);
	c2.Divide(1,1,0.01,0.01);
	c2.cd(1);
	c2.GetPad(1)->SetLogy();
	histoCandidatePt.Draw();
	c2.Print("histoCandidatePt.pdf");

	TCanvas c3;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c3.SetCanvasSize(1000,650);
	c3.Divide(1,1,0.01,0.01);
	c3.cd(1);
	histoTheta.Draw();
	c3.Print("PhotonTheta.pdf");

	TCanvas c4;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c4.SetCanvasSize(1000,650);
	c4.Divide(1,1,0.01,0.01);
	c4.cd(1);
	c4.SetLogy();
	histoRemainingEnergy.Draw();
	c4.Print("histoRemainingEnergy.pdf");

	TCanvas c5;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c5.SetCanvasSize(1000,650);
	c5.Divide(1,1,0.01,0.01);
	c5.cd(1);
	histoCandidateEnergy.Draw();
	c5.Print("histoCandidateEnergy.pdf");

	TCanvas c6;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c6.SetCanvasSize(1000,650);
	c6.Divide(1,1,0.01,0.01);
	c6.cd(1);
	c6.GetPad(1)->SetLogy();
	histoCanditateTheta.Draw();
	c6.Print("histoCandidateTheta.pdf");

	TCanvas c7;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c7.SetCanvasSize(1000,650);
	c7.Divide(1,1,0.01,0.01);
	c7.cd(1);
	c7.GetPad(1)->SetLogy();
	histoCandidatePhi.Draw();
	c7.Print("histoCandidatePhi.pdf");

	TCanvas c8;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c8.SetCanvasSize(1000,650);
	c8.Divide(1,1,0.01,0.01);
	c8.cd(1);
	c8.GetPad(1)->SetLogy();
	histoBoost.Draw();
	c8.Print("histoCandidateBoost.pdf");

	TCanvas c9;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c9.SetCanvasSize(1000,650);
	c9.Divide(1,1,0.01,0.01);
	c9.cd(1);
	//c9.GetPad(1)->SetLogy();
	histoHigherEnergyPhoton.Draw();
	c9.Print("histoHigherEnergyPhoton.pdf");

	TCanvas c10;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c10.SetCanvasSize(1000,650);
	c10.Divide(1,1,0.01,0.01);
	c10.cd(1);
	//c10.GetPad(1)->SetLogy();
	histoLowerEnergyPhoton.Draw();
	c10.Print("histoLowerEnergyPhoton.pdf");

	TCanvas c11;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c11.SetCanvasSize(1000,650);
	c11.Divide(1,1,0.01,0.01);
	c11.cd(1);
	c11.GetPad(1)->SetLogy();
	histoZbirPt.Draw();
	c11.Print("histoZbirPt.pdf");

	TCanvas c12;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c12.SetCanvasSize(1000,650);
	c12.Divide(1,1,0.01,0.01);
	c12.cd(1);
	//c12.GetPad(1)->SetLogy();
	histoHelicityAngle.Draw();
	c12.Print("histoHelicityAngle.pdf");

	TCanvas c13;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c13.SetCanvasSize(1000,650);
	c13.Divide(1,1,0.01,0.01);
	c13.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoconeEnergy.Draw();
	c13.Print("Cone_Energy.pdf");

	TCanvas c14;
	gStyle->SetOptStat(111111);
	gStyle->SetPalette( 1 );
	c14.SetCanvasSize(1000,650);
	c14.Divide(1,1,0.01,0.01);
	c14.cd(1);
	//c13.GetPad(1)->SetLogy();
	energyOfPhotons.Draw();
	c14.Print("energyOfPhotons.pdf");

	TCanvas c15;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c15.SetCanvasSize(1000,650);
	c15.Divide(1,1,0.01,0.01);
	c15.cd(1);
	//c15.GetPad(1)->SetLogy();
	histoPhi.Draw();
	c15.Print("PhotonPhi.pdf");

	TCanvas c16;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c16.SetCanvasSize(1000,650);
	c16.Divide(1,1,0.01,0.01);
	c16.cd(1);
	//c15.GetPad(1)->SetLogy();
	histogram.Draw();//colz
	c16.Print("ISR.pdf");

	TCanvas c17;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c17.SetCanvasSize(1000,650);
	c17.Divide(1,1,0.01,0.01);
	c17.cd(1);
	//c15.GetPad(1)->SetLogy();
	histoangleBetweenPhotons.Draw();
	c17.Print("angleBetweenPhotons.pdf");

	TCanvas c18;
	gStyle->SetPalette( 1 );
	c18.SetCanvasSize(1000,650);
	c18.Divide(1,1,0.01,0.01);
	c18.cd(1);
	//c18.GetPad(1)->SetLogy();
	gStyle->SetOptStat(111111);
	histoNumberPhotonsbyEvent.Draw();
	c18.Print("NumberPhotonsbyEvnt.pdf");

	TCanvas c19;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c19.SetCanvasSize(1000,650);
	c19.Divide(1,1,0.01,0.01);
	c19.cd(1);
	//c19.GetPad(1)->SetLogy();
	anglePhotonParticle.Draw();
	c19.Print("PhotonParticleAngle.pdf");

	TCanvas c20;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c20.SetCanvasSize(1000,650);
	c20.Divide(1,1,0.01,0.01);
	c20.cd(1);
	//c20.GetPad(1)->SetLogy();
	anglePhotonParticlewE.Draw();
	c20.Print("PhotonParticleAngleWE.pdf");

	TCanvas c21;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c21.SetCanvasSize(1000,650);
	c21.Divide(1,1,0.01,0.01);
	c21.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoconeEnergyFilter.Draw();
	c21.Print("Cone_Energy_Filter.pdf");

	TCanvas c22;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c22.SetCanvasSize(1000,650);
	c22.Divide(1,1,0.01,0.01);
	c22.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoHardPhotonsByEvent.Draw();
	c22.Print("HardPhotonsByEvent.pdf");

	TCanvas c23;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c23.SetCanvasSize(1000,650);
	c23.Divide(1,1,0.01,0.01);
	c23.cd(1);
	c23.GetPad(1)->SetLogy();
	histoPtOtherParticles.Draw();
	c23.Print("histoPtOtherParticles.pdf");

	TCanvas c24;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c24.SetCanvasSize(1000,650);
	c24.Divide(1,1,0.01,0.01);
	c24.cd(1);
	//c24.GetPad(1)->SetLogy();
	histoMissingEnergy.Draw();
	c24.Print("histoMissingEnergy.pdf");

	TCanvas c25;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c25.SetCanvasSize(1000,650);
	c25.Divide(1,1,0.01,0.01);
	c25.cd(1);
	//c24.GetPad(1)->SetLogy();
	histoTestPt.Draw();
	c25.Print("histoTestPt.pdf");


	TCanvas c26;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c26.SetCanvasSize(1000,650);
	c26.Divide(1,1,0.01,0.01);
	c26.cd(1);
	c26.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	efikasnostproba.Draw("AL*");
	c26.Print("efikasnostPhotonPT.pdf");

	TCanvas c27;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c27.SetCanvasSize(1000,650);
	c27.Divide(1,1,0.01,0.01);
	c27.cd(1);
	//c27.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histoPhoton1CandidateTheta.Draw();
	c27.Print("Thetaof1stPhotonofCandidate.pdf");

	TCanvas c28;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c28.SetCanvasSize(1000,650);
	c28.Divide(1,1,0.01,0.01);
	c28.cd(1);
	//c28.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histoPhoton2CandidateTheta.Draw();
	c28.Print("Thetaof2ndPhotonofCandidate.pdf");

	TCanvas c29;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c29.SetCanvasSize(1000,650);
	c29.Divide(1,1,0.01,0.01);
	c29.cd(1);
	//c28.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histoPtof1stPhotonofCandidate.Draw();
	c29.Print("Ptof1stPhotonofCandidate.pdf");


	TCanvas c30;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c30.SetCanvasSize(1000,650);
	c30.Divide(1,1,0.01,0.01);
	c30.cd(1);
	//c28.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histoPtof2ndPhotonofCandidate.Draw();
	c30.Print("Ptof2ndPhotonofCandidate.pdf");

	TCanvas c31;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c31.SetCanvasSize(1000,650);
	c31.Divide(1,1,0.01,0.01);
	c31.cd(1);
	//c28.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histoPtof2ndPhotonofCandidate.Draw();
	c31.Print("histoHighestPhotonPt.pdf");

	TCanvas c32;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c32.SetCanvasSize(1000,650);
	c32.Divide(1,1,0.01,0.01);
	c32.cd(1);
	//c28.SetTitle("Signal");
		//c24.GetPad(1)->SetLogy();
	histo2ndHighestPhotonPt.Draw("");
	c32.Print("histo2ndHighestPhotonPt.pdf");

	TCanvas c33;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c33.SetCanvasSize(1000,650);
	c33.Divide(1,1,0.01,0.01);
	c33.cd(1);
	//c33.SetTitle("Signal");
		//c33.GetPad(1)->SetLogy();
	histo2ndHighestPhotonPtZoomed.Draw();
	c33.Print("histo2ndHighestPhotonPtZoomed.pdf");

	TCanvas c34;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c34.SetCanvasSize(1000,650);
	c34.Divide(1,1,0.01,0.01);
	c34.cd(1);
	//c33.SetTitle("Signal");
		//c33.GetPad(1)->SetLogy();
	histoVisibleEnergy.Draw();
	c34.Print("histoVisibleEnergy.pdf");

	TCanvas c35;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(0);
	c35.SetCanvasSize(1000,650);
	c35.Divide(1,1,0.01,0.01);
	c35.cd(1);
	//c33.SetTitle("Signal");
		//c33.GetPad(1)->SetLogy();
	histoZbirImpulsa.Draw();
	c35.Print("histoZbirImpulsa.pdf");


	TCanvas c36;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c36.SetCanvasSize(1000,650);
	c36.Divide(1,1,0.01,0.01);
	c36.cd(1);
	//c33.SetTitle("Signal");
	c36.GetPad(1)->SetLogy();
	histoCosinusHelicityAngle.Draw();
	c36.Print("cosinusHelicityAngle.pdf");


	TCanvas c37;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c37.SetCanvasSize(1000,650);
	c37.Divide(1,1,0.01,0.01);
	c37.cd(1);
	//c33.SetTitle("Signal");
	//c37.GetPad(1)->SetLogy();
	histoCandidateMFineBinning.Draw();
	c37.Print("CandidateMFineBinning.pdf");

	TCanvas c38;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c38.SetCanvasSize(1000,650);
	c38.Divide(1,1,0.01,0.01);
	c38.cd(1);
	//c33.SetTitle("Signal");
	c38.GetPad(1)->SetLogy();
	ptOfPhotons.Draw();
	c38.Print("ptOfPhotons.pdf");

	TCanvas c39;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c39.SetCanvasSize(1000,650);
	c39.Divide(1,1,0.01,0.01);
	c39.cd(1);
	//c33.SetTitle("Signal");
	c39.GetPad(1)->SetLogy();
	PtofCandidatePhotons.Draw();
	c39.Print("PtofCandidatePhotons.pdf");

	TCanvas c40;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c40.SetCanvasSize(1000,650);
	c40.Divide(1,1,0.01,0.01);
	c40.cd(1);
	//c33.SetTitle("Signal");
	c40.GetPad(1)->SetLogy();
	energyofCandidatePhotons.Draw();
	c40.Print("energyofCandidatePhotons.pdf");

	TCanvas c41;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c41.SetCanvasSize(1000,650);
	c41.Divide(1,1,0.01,0.01);
	c41.cd(1);
	//c33.SetTitle("Signal");
	c41.GetPad(1)->SetLogy();
	thetaOfCandidatePhotons.Draw();
	c41.Print("ThetaofCandidatePhotons.pdf");

	TCanvas c42;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c42.SetCanvasSize(1000,650);
	c42.Divide(1,1,0.01,0.01);
	c42.cd(1);
	//c33.SetTitle("Signal");
	c42.GetPad(1)->SetLogy();
	energyofCandidatePhotonsZoomed.Draw();
	c42.Print("energyofCandidatePhotonsZoomed.pdf");

	TCanvas c43;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c43.SetCanvasSize(1000,650);
	c43.Divide(1,1,0.01,0.01);
	c43.cd(1);
	//c33.SetTitle("Signal");
	c43.GetPad(1)->SetLogy();
	PtofCandidatePhotonsZoomed.Draw();
	c43.Print("PtofCandidatePhotonsZoomed.pdf");

	TCanvas c44;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c44.SetCanvasSize(1000,650);
	c44.Divide(1,1,0.01,0.01);
	c44.cd(1);
	//c33.SetTitle("Signal");
	c44.GetPad(1)->SetLogy();
	histoconeEnergyfor2Photons.Draw();
	c44.Print("histoconeEnergyfor2Photons.pdf");
*/

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

	TString rfName = "InvM3TeV_Cone_angle.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
