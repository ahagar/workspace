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


#include "varList.h"
#include "utils.h"

using namespace std;

const Double_t minPt1 = 15;			   //minimalna energija fotona
const Double_t minPt2 = 20;			   //minimalna energija fotona
const Double_t minPt3 = 25;			   //minimalna energija fotona
const Double_t minPt4 = 30;			   //minimalna energija fotona

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
	




	TTree eventList("eventsSignal", "ILD event list");
	varListGoran vl; /* SL specific */
	eventList.Branch("CandidateM", &(vl.CandidateM), "CandidateM");
	eventList.Branch("CandidateEnergy", &(vl.CandidateEnergy), "CandidateEnergy");
	eventList.Branch("CandidateTheta", &(vl.CandidateTheta), "CandidateTheta");
	eventList.Branch("CandidatePt", &(vl.CandidatePt), "CandidatePt");
	eventList.Branch("CandidatePhi", &(vl.CandidatePhi), "CandidatePhi");
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
	eventList.Branch("pTCandidatePhotons2", &(vl.pTCandidatePhotons2), "pTCandidatePhotons2");
	eventList.Branch("energyCandidatePhoton1", &(vl.energyCandidatePhotons1), "energyCandidatePhoton1");
	eventList.Branch("energyCandidatePhoton2", &(vl.energyCandidatePhotons2), "energyCandidatePhoton2");
	eventList.Branch("neutrinoE", &(vl.neutrinoE), "neutrinoE");
	eventList.Branch("antineutrinoE", &(vl.antineutrinoE), "antineutrinoE");

	Float_t isr_photon_energy=0;
	Float_t isr_photon_pt=0;
	Float_t isr_photon_theta=0;

	Float_t real_photon_energy=0;
	Float_t real_photon_pt=0;
	Float_t real_photon_theta=0;

	Float_t photon_x_energy = 0;
	Float_t photon_x_pt = 0;
	Float_t photon_x_theta = 0;

	TTree isr("isr", "ILD event list");
	isr.Branch("isr_photon_energy", &isr_photon_energy, "isr_photon_energy");
	isr.Branch("isr_photon_pt", &isr_photon_pt, "isr_photon_pt");
	isr.Branch("isr_photon_theta", &isr_photon_theta, "isr_photon_theta");


	TTree photonReal("photonReal", "ILD event list");
	photonReal.Branch("real_photon_energy", &real_photon_energy, "real_photon_energy");
	photonReal.Branch("real_photon_pt", &real_photon_pt, "real_photon_pt");
	photonReal.Branch("real_photon_theta", &real_photon_theta, "real_photon_theta");

	TTree photon_x("photon_x", "ILD event list");
	photon_x.Branch("photon_x_energy", &photon_x_energy, "photon_x_energy");
	photon_x.Branch("photon_x_pt", &photon_x_pt, "photon_x_pt");
	photon_x.Branch("photon_x_theta", &photon_x_theta, "photon_x_theta");

	TH1F histoPhotons7 ("histoPhotons7", "histoPhotons7", 10, 0, 10);
	TH1F histoPhotons15 ("histoPhotons15", "histoPhotons15", 10, 0, 10);







	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	Double_t coneEnergy = 0; //Energija konusa oko fotona
	Double_t theta = 0;
	Double_t angle_Between_Photons=0;


	Int_t totalCutEvents = 0;
	Int_t totalCutEvents_ecal = 0;
	Int_t totalCutEvents_hcal = 0;
	Int_t primaryCutEvents=0;
	Int_t totalEvents = 0;
	Int_t counterPhotonsPt1 = 0;
	Int_t counterEventsPt2 = 0;
	Int_t counterEventsPt3 = 0;
	Int_t counterEventsPt4 = 0;


	Int_t CounterStrike = 0;



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
		Int_t numCutEvents = 0; 	//broj dogadjaja posle preselekcije cut-ova
		Int_t numCutEvents_ecal = 0; 	//broj dogadjaja posle ecal preselekcije cut-ova
		Int_t numCutEvents_hcal = 0; 	//broj dogadjaja posle hcal preselekcije cut-ova



		Int_t eventsPrimaryCut = 0; 	//broj dogadjaja posle određenih cut-ova


		std::vector <TLorentzVector> allPhotons; //vektor koji prikuplja sve fotone
		Int_t brojac =0;


		Int_t counterAllCuts = 0;
		//Int_t hardPhotons = 0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{


			Double_t ecal = 0;
			Double_t ecalEvt = 0;
			Double_t ecalPht = 0;

			Double_t hcal = 0;
			Double_t hcalPht = 0;

			Double_t totalCal = 0;
			Double_t totalCalEvt = 0;
			Double_t calRatioEvt = 0;

			Double_t brojac_fotona_7 =0;
			Double_t brojac_fotona_15 =0;




			brojac++;
			vector<TLorentzVector> photons;
			vector<Double_t> photonsECal;
			vector<Double_t> photonsHCal;

			vector<TLorentzVector> particles;   // other than photons
			vector <Double_t> PtPhotons;
			vector <Double_t> ConeEnergyOfPhotons;
			vector <EVENT::ReconstructedParticle*> niz_fotoni;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
		//		std::vector<std::string> colNames = *evt -> getCollectionNames();

				EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("SelectedPandoraPFANewPFOs");/*("SelectedPandoraPFANewPFOs_Reprocess");*/
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");

			double_t Evis = 0;	//Ukupna energija po dogadjaju PandoraPFOsDefault_Reprocess, LooseSelectedPandoraPFANewPFOs_Reprocess,
			//SelectedPandoraPFANewPFOs_Reprocess, TightSelectedPandoraPFANewPFOs_Reprocess

			bool leptonFound = false;	//uslov da nemamo leptone
			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
			{

				 for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
						   {
						IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

						if (mcParticle->getGeneratorStatus() !=1) continue;

						if (mcParticle->getPDG()==12) {
							vl.neutrinoE = mcParticle->getEnergy();
						}
						if (mcParticle->getPDG()()==-12) {
							vl.antineutrinoE = mcParticle->getEnergy();
						}

						   }//end po MCParticles

				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici



				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(recParticle->getType());
			//	cout<<"particle PDG is: "<<particlePDG<<endl;
				std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) recParticle->getClusters();
				for ( std::vector<EVENT::Cluster*>::const_iterator iCluster=clusters.begin(); iCluster!=clusters.end(); ++iCluster)
				{

					if (*iCluster)
					{
						ecal += (*iCluster)->getSubdetectorEnergies()[0];
						hcal += (*iCluster)->getSubdetectorEnergies()[1];
					}

				}

				totalCalEvt = ecal+hcal;


				if(particlePDG == 22)	//rad sa fotonima (PDG=22)
				{


					theta = temp.Theta()*180/M_PI; //promenljiva koja nam daje Theta čestice
					//Double_t phi = temp.Phi();		//promenljiva koja nam daje Phi
					Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice

					PtPhotons.push_back(temp.Pt());
					allPhotons.push_back(temp);//sakupili smo sve fotone

					if (temp.Pt()>7){
						brojac_fotona_7++;
					}
					if (temp.Pt()>15){
						brojac_fotona_15++;
					}


					//uzimamo u obzir samo one fotone koji nam prodju uslove
					if(temp.Pt() >7 )//minpT1
					{

						counterPhotonsPt1++;

						TLorentzVector currentConeAxis = temp; //foton oko koga pravimo konus
						Double_t coneEnergy = 0; //Energija konusa oko fotona


						for (Int_t k = 0; k < recParticles->getNumberOfElements() ; k++) /*&& coneEnergy <= maxConeEnergy -zaustavlja fot loop kada predje maxE */
						{
							IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(k);
							const double *impuls = recParticle->getMomentum();
							double energija = recParticle->getEnergy();
							TLorentzVector otherParticle;
							otherParticle.SetPxPyPzE(impuls[0], impuls[1], impuls[2],energija);//četvorovektor drugih čestica
							if (currentConeAxis == otherParticle) continue; //da ne bi ubrojali foton oko kog pravimo konus
                            Double_t particleAngle = currentConeAxis.Angle(otherParticle.Vect())* 180 / TMath::Pi() ;//ugao izmedju fotona oko kog pravimo konus i čestice
							//anglePhotonParticle.Fill(particleAngle);//histogram uglova izmedju fotona i čestice
							//anglePhotonParticlewE.Fill(particleAngle, otherParticle.E());//histogram uglova izmedju fotona i čestice otežinjen energijom te čestice
							if (recParticle->getType()==22)
							{
								angle_Between_Photons =particleAngle;
							}

							//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
							if (particleAngle <=coneAngle)
							{
								coneEnergy += otherParticle.E();
							}


						}
					//	histoconeEnergy.Fill(coneEnergy); //histogram koji iscrtava energiju konusa
						ConeEnergyOfPhotons.push_back(coneEnergy);




						//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
						if (coneEnergy <= maxConeEnergy )
						{

								counterAllCuts++;//brojač za fotone koji prodju sve cutove

							/*	if (photons.size() == 2)
									break;
							*/

								std::vector<EVENT::Cluster*> clustersPht = (std::vector<EVENT::Cluster*>) recParticle->getClusters();
								//cout<<"rekonstrisana čestica je:"<<recParticle->getType()<<endl;

								 for ( std::vector<EVENT::Cluster*>::const_iterator kCluster=clustersPht.begin(); kCluster!=clustersPht.end(); ++kCluster)
								{

								if (*kCluster) {
									ecalPht += (*kCluster)->getSubdetectorEnergies()[0];
									hcalPht += (*kCluster)->getSubdetectorEnergies()[1];
								}

								 }
								 totalCal=ecalPht + hcalPht;



								photons.push_back(temp);
								photonsECal.push_back(ecal);
								photonsHCal.push_back(hcal);
								niz_fotoni.push_back(recParticle);
							//	histoconeEnergyFilter.Fill(coneEnergy);//histogram koji iscrtava energiju konusa



						}

					}


				} //pdg = 22
				//dogadjaj ne prolazi ako se u dogadjaju detektuje lepton
				else if (particlePDG != 22)
					{

					if(temp.Pt() > leptonPt)
						leptonFound = true;
				//	histoPtOtherParticles.Fill(temp.Pt());
					}
				if( temp.Pt() > 5)
				{
					Evis += e ;//ukupna enrgija po dogadjaju
				}



			}  // end of particle loop

			//ovde punimo histo fotona po evt
			histoPhotons15.Fill(photons.size());
			histoPhotons7.Fill(brojac_fotona_15);



			bool candFound = false;//
			int candCounter = 0;
			for (int i = 0; i < (int)allPhotons.size() - 1; ++i)
			{
				for (int j = i + 1; j <(int) allPhotons.size(); ++j)
				{

					TLorentzVector pair = allPhotons[i]+allPhotons[j];//sabira četvorovektore

					if (110 > pair.M() && pair.M() < 140)//uslov da bi par bio higs kandidat
					{

					 candFound = true;
					// candCounter++;
				    // histoCandidateMTest.Fill(pair.M());
				     //if (candFound==true) break;

					}
				}
			}

			if  (candFound == false) continue;



			sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());

		//	histoHighestPhotonPt.Fill(PtPhotons[0]);
		//	histo2ndHighestPhotonPt.Fill(PtPhotons[1]);
			vl.HighestPhotonPt = PtPhotons[0];
			vl.secondHighestPhotonPt = PtPhotons[1];



			if (ConeEnergyOfPhotons.size()==2)
			{
		//	histoconeEnergyfor2Photons.Fill(ConeEnergyOfPhotons[0]);
		//	histoconeEnergyfor2Photons.Fill(ConeEnergyOfPhotons[1]);
			}

			Int_t photonsSize = photons.size();
		//	histoNumberPhotonsbyEvent.Fill(photons.size());
			vl.NumberPhotonsbyEvent = photonsSize;


			Int_t Emiss = EnergyCenterMass - Evis; // missing energy
		//	histoMissingEnergy.Fill(Emiss);

			vl.MissingEnergy = Emiss;


			//cout << "lepton found: " <<leptonFound <<endl;
			if(leptonFound) continue; //ako je u finalnom stanju imamo lepton ili kvark (sem neutrina), događaj se preskače


			// At this point, all relevant photons have been collected.


		//	cout << "photon size is: " << photons.size()<<endl;
			//proveravamo da li u eventu ima detektovano više od jednog fotona koji bi bili kandidat za higsov bozon
			if (photons.size() == photonSize )
			{
				vector<CandidateData> candidates;
				CounterStrike++;

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
				//	numCutEvents++;//brojimo dogadjaje koji prodju uslove za energiju, invarijantnu masu, uglove,Pt
					eventsPrimaryCut++;
					TLorentzVector higgs = candidate->Higgs();
					Double_t Ehiggs = higgs.Energy();//energija kandidata


					Double_t Eremaining = Evis - Ehiggs;//vidljiva energija posle higgsa
					//		cout << "E remaining:"<< Eremaining<<endl;

					Double_t higgsInvM=higgs.M();
					Double_t higgsPt =higgs.Pt();
					TLorentzVector photon ;
					TLorentzVector photon2;
					//histoCandidateM.Fill(higgs.M());
					if (candidate->Photon1.Pt() > candidate->Photon2.Pt()) {
						photon = candidate->Photon1;
						photon2 = candidate->Photon2;
					}else{
						photon = candidate->Photon2;
						photon2 = candidate->Photon1;
					}

				    Double_t theta1 = photon.Theta()*180/M_PI;
				    Double_t theta2 = photon2.Theta()*180/M_PI;
				    Double_t pt1 = photon.Pt();
				    Double_t pt2 = photon2.Pt();

/*
						cout << "energija  prvog fotona je: " <<photon.E()<<endl;
						cout << "energija  drugog fotona je: " <<photon2.E()<<endl;
						cout << "px1 " <<photon.Pt()<<endl;
						cout << "px 2" <<photon2.Pt()<<endl;

						for (int p = 0; p < (int)allPhotons.size()-1; p++){
							if (allPhotons[p] == photon || allPhotons[p] == photon2) cout<<"našao sam isti foton\n";
						}*/



				    calRatioEvt = ecalEvt/totalCalEvt;



				    vl.CandidateM=higgsInvM;
				    vl.CandidateEnergy = higgs.E();
				    vl.CandidatePhi = higgs.Phi();
				    vl.CandidatePt = higgs.Pt();
				    vl.CandidateTheta = higgs.Theta();
				    vl.pLCandidate = higgs.Z();
				    vl.thetaCandidatePhotons1 = photon.Theta();
				    vl.thetaCandidatePhotons2 = photon2.Theta();
				    vl.pLCandidatePhotons1 = abs(photon.Z());
				    vl.pLCandidatePhotons2 = abs(photon2.Z());
				    vl.pTCandidatePhotons1 = photon.Pt();
				    vl.pTCandidatePhotons2 = photon2.Pt();
				    vl.energyCandidatePhotons1 = photon.E();
				    vl.energyCandidatePhotons2 = photon2.E();

				//    cout<<"theta pre boosta: "<< higgs.Theta()<<endl;


					bool eventTrue =false;
					bool eventISR =false;
					//
						IMPL::LCCollectionVec* mcPar = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticlesSkimmed");
						EVENT::LCRelation* link_to_rec_photon = 0;
	 					EVENT::MCParticle* pointer_to_mc_photon = 0;
	 					EVENT::ReconstructedParticle* pointer_to_rec_photon = 0;

	 				//	pointer_to_rec_photon =niz_fotoni[i];

	 					EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");
	 					for (int j = 0; j < links -> getNumberOfElements(); j++){

		 					EVENT::LCRelation* linki = (EVENT::LCRelation*) links -> getElementAt(j);
	 						EVENT::MCParticle* mcpj = (EVENT::MCParticle*) linki -> getTo();
							EVENT::ReconstructedParticle* rpi1 = (EVENT::ReconstructedParticle*) linki -> getFrom();
						//	cout <<"našao čestice sa istom energijom."<<endl;
							if (rpi1 ==0) continue;
							TLorentzVector temp (TVector3 (rpi1 -> getMomentum()), rpi1 -> getEnergy());
						//	if (rpi1->getType()==22 && rpi1->getEnergy() == niz_fotoni[0]->getEnergy())cout <<"našao čestice sa istom energijom."<<endl;
						//	cout<<"linkovana cestica je: "<<niz_fotoni[i]->getType()<<endl;

							int counter_ISR = 0;
							for (int k =0; (int)k< niz_fotoni.size();k++){
							if (rpi1 == niz_fotoni[k]){
							//	cout<<"linkovana cestica je: "<<rpi1->getType()<<endl;
								if (mcpj == 0 || mcpj->getPDG() !=22 || linki->getWeight()<0.985) continue;
								eventTrue=true;
								if(mcpj->getParents().size()== 0) continue;
								const EVENT::MCParticleVec & parent = mcpj -> getParents();

						//IF1
								//if (fabs(mcpj->getParents()[0]->getPDG())==22 &&  fabs(mcpj->getParents()[0]->getParents()[0]->getPDG())==11  && mcpj->getParents()[0]->getParents()[0]->getEnergy() > 1499.9 )
								if ( parent[0]->getPDG()==22 && parent[0]->getParents()[0]->getParents().size()==0 /*&& parent[0]->getParents()[0]->getParents()[0]->getEnergy()>1499*/)
								{
							//	cout<< "gen cestica ima koliko roditelja: "<<mcpj->getParents()[0]->getPDG()<<endl;
							//	cout << "energija roditelja je : "<< mcpj->getParents()[0]->getParents()[0]->getEnergy()<<endl;
						//			cout<< "baba je: "<< parent[0]->getParents()[0]->getPDG()<<", energija babe je: " <<parent[0]->getParents()[0]->getEnergy()<<endl;
							//		cout<<"_______________________________"<<endl;
									eventISR  = true;

									counter_ISR++;
							}//end of IF1

								bool realPhoton= false;

								// provera da li je foton pravi iz procesa (eeGAMMA), a ne isr ili neki slucajni
								if ( parent[0]->getPDG()==22 && parent[0]->getParents()[0]->getParents().size()!=0 && parent[0]->getParents()[0]->getParents()[0]->getEnergy()>1500 /* && parent[0]->getParents()[0]->getParents()[0]->getEnergy()>1499*/)
								{
						//		cout<< "da li postoji prababa: "<<parent[0]->getParents()[0]->getParents().size()<<endl;

								realPhoton=true;
								bool testIsr = false;
								bool testgamma1 = false;
								bool testgamma2 = false;

								if (realPhoton){
									real_photon_energy = temp.E();
									real_photon_pt = temp.Pt();
								    real_photon_theta = temp.Theta();

								    photonReal.Fill();


									if ( temp == photon ) testgamma1=true;
									if ( temp == photon2) testgamma2=true;


								}

								if (testgamma1){
									photon_x_energy = photon2.E();
									photon_x_pt = photon2.Pt();
								}

								if (testgamma2){
									photon_x_energy = photon.E();
									photon_x_pt = photon.Pt();
								}
								 photon_x.Fill();



							}//end of pravi foton

			 					if(eventISR ==false)continue;
			 				//	cout<<"isr bool: "<<eventISR<<", energija: "<<isr_photon_energy<<endl;
								isr_photon_energy = temp.E();
								isr_photon_pt = temp.Pt();
								isr_photon_theta = temp.Theta();

			 				if(counter_ISR !=1) continue;
			 				isr.Fill();

							}//end of rpi1 == niz_fotoni[k]
							}//end of for niz_fotoni.size

	 					}//end of links
	 				//	cout <<"redni broj događaja je: "<<brojac<<endl;
	 					if (eventTrue==false)continue;
	 					if(eventISR == true)continue;

	 					TVector3 boosttoparent = -(higgs.BoostVector());//prelazimo u sitem CM


	 						photon.Boost(boosttoparent);//prebacujemo fotone u sistem CM
	 						photon2.Boost(boosttoparent);
	 						TVector3 photon_3v=photon.Vect();
	 						TVector3 photon2_3v = photon2.Vect();
	 						TVector3 higgs_3v = higgs.Vect();



	 						Double_t anglePhoton1 = photon.Angle(higgs_3v);
	 						Double_t anglePhoton2 = photon2.Angle(higgs_3v);

	 					//	cout<<"higgs posle boosta: "<<higgs_3v.Theta()<<endl;

	 					//	cout<<cos(anglePhoton1)<< " ,"<<cos(anglePhoton2)<<endl;
	 					//	cout<<"________________________________________________"<<endl;

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
	 				//		eventList.Fill(); //KG puni se events

	 				//		cout <<"redni broj dogadjaja signala je: "<< CounterStrike<<endl;

	 						//	 						cout << "____________________________________________ " <<endl;



	 					eventList.Fill(); //dodato da bi se upisali samo evt koji nemaju ISR foton




				}//end of Candidates. ovo je visak. jer imamo samo 2 fotona
			}//end of photonSize ==2


		} // End of event loop


		totalCutEvents+=numCutEvents;//ukupan broj dogadjaja koji prodje cutove iz svih fajlova
		totalCutEvents_ecal+=numCutEvents_ecal;//ukupan broj dogadjaja koji prodje ecal preselekciju iz svih fajlova
		totalCutEvents_hcal+=numCutEvents_hcal;//ukupan broj dogadjaja koji prodje  hcal preseleckciju iz svih fajlova


		totalEvents+=brojDogadjaja; //ukupan broj dogadjaja iz svih fajlova
		primaryCutEvents+=eventsPrimaryCut;
		cout << "Broj event-a posle cut-a je: " << numCutEvents << endl;
		//cout << "Broj event-a posle ecall cut-a je: " << numCutEvents << endl;
		//cout << "Broj event-a posle hcall cut-a je: " << numCutEvents << endl;

		//cout << "Broj fotona posle primarnog cut-a  je: " << counterPhotonsPt1 << endl;
		//cout << "Broj fotona posle svih cut-ova  je: " << counterAllCuts << endl;
		cout << "Broj dogadjaja je: " << totalEvents << endl;





		lcReader->close();

	} // End of file loop



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	histoPhotons7.Write();
//	histoPhotons15.Write();
	eventList.Write();
	isr.Write();
	photonReal.Write();
	 photon_x.Write();






	rootFile.Close();



	return 0;
}



Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;
	UInt_t nFirstJob = 1;
	if(argc>iarg) nFirstJob = atoi(argv[iarg]); iarg++;
	UInt_t nLastJob = 535;
	if(argc>iarg) nLastJob   = atoi(argv[iarg]); iarg++;

	TString fName = "hvv_";//ee_rem_col_"
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "hgg.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
