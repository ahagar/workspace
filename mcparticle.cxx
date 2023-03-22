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
	#include "TGraph.h"
	#include "math.h"

	// LCIO includes
	#include <lcio.h>
	#include <IOIMPL/LCFactory.h>
	#include <IMPL/LCCollectionVec.h>
	#include <EVENT/MCParticle.h>
	#include <EVENT/ReconstructedParticle.h>
	#include <IMPL/CalorimeterHitImpl.h>
	#include <IMPL/MCParticleImpl.h>
	#include <IMPL/ReconstructedParticleImpl.h>
	#include <UTIL/LCRelationNavigator.h>
	#include <EVENT/LCRelation.h>
	#include <UTIL/LCTOOLS.h>
	#include <Exceptions.h>
	#include <IMPL/MCParticleImpl.h>
	#include <pre-generated/EVENT/LCRelation.h>
#endif

#include <vector>
#include <array>

using namespace std;
#include "stdlib.h"
#include "CandidateData.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include "varList.h"



using namespace std;

const Double_t minPt1 = 10;			   //minimalna energija fotona
const Double_t minPt2 = 15;			   //minimalna energija fotona
const Double_t minPt3 = 18;			   //minimalna energija fotona
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
	TH1F histoTheta ("histoTheta", "Theta; #theta_{#gamma}", 180, 0, 180);
	TH1F histoPhi ("histoPhi", "Phi; #phi_{#gamma}", 50, 0, 3.14);
	TH1F histoPhotonEnergy ("histoPhotonEnergy", "PhotonEnergy; E_{#gamma} (GeV)", 150, 0, 1500);
	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{#gamma#gamma} (GeV)", 30, 110, 140);
	TH1F histoCandidatePt ("CandidatePt", "CandidatePt; Pt_(#gamma#gamma) (GeV)", 150, 0, 1500);
	TH1F histoRemainingEnergy ("Remaining_Visible_Energy", "Remaining_Visible_Energy; E_{remaining}(GeV)", 150, 0, 3000);
	TH1F histoCandidateEnergy ("CandidateE", "CandidateE; E_{H} (GeV)",150, 0, 1550);
	TH1F histoCanditateTheta ("CandidateTheta", "CandidateTheta; #theta_{#gamma#gamma}", 90, 0, 180);
	TH1F histoCandidatePhi ("CandidatePhi", "CandidatePhi' #phi_{#gamma}", 50, 0, 3.14);
	TH1F histoBoost ("CandidateBoost", "CandidateBoost; #beta_{#gamma#gamma}", 180, 0, 1);
	TH1F histoZbirPt ("ZbirPt", "ZbirPt; Pt_{1} + Pt_{2}", 300, 0, 3000);
	TH1F histoHigherEnergyPhoton ("HigherEnergyPhoton", "HigherEnergyPhoton; E_{#gamma1}", 150, 0, 1500);
	TH1F histoLowerEnergyPhoton ("LowerEnergyPhoton", "LowerEnergyPhoton; E_{#gamma2}", 150, 0, 1500);
	TH1F histohelicityAngle ("helicityAngle", "helicityAngle; cos#theta", 100, 0, 1);
	TH1F histoarcCos ("arcCos", "arcCOs; arcCos#theta", 300, 0, 1);
	TH2F histogram ("test","title", 20, 0, 3.14, 20, 0, 3.14);
	TH1F histoangleBetweenPhotons ("angleBetweenPhotons", "angleBetweenPhotons; #alpha", 50, 0, 3.14);
	TH1F histoNumberPhotonsbyEvent ("Number_of_Photons_by_Event", "Number_of_Photons_by_Event", 10, 0, 10);
	TH1F anglePhotonParticle ("Angle Between photon and particle", "angle photon particle", 50, 0, 20);
	TH1F histoconeEnergy ("Cone_Energy", "Cone_Energy; E (GeV)", 100, 0, 100);
	TH1F anglePhotonParticlewE ("Angle_Between_photon_and_particle_wE", " photon_particle_angle_wE", 50, 0, 20);//otezinjeno sa Energijom
	TH1F histoconeEnergyFilter ("Cone Energy", "Cone Energy; E (GeV)", 100, 0, 10);
	TH1F histoHardPhotonsByEvent ("Hard Photons by Event", "Hard Photons by Event", 10, 0, 10);
	TH1F histoPtOtherParticles ("PtOtherParticles", "PtOtherParticles; Pt (GeV)", 100, 0, 20);
	TH1F histoMissingEnergy ("Missing_Energy", "Missing_Energy; E_{miss} (GeV)", 100, 0, 3000);
	TH1F histoTestPt ("Test Pt", "Test Pt; Pt (GeV)", 50, 0, 50);
	TH1F histoTheta1Photon ("histoTheta", "ThetaofcandidatePhoton; #theta_{#gamma}", 50, 0, 3.14);
	TH1F histoTheta2Photon ("histoTheta", "Theta; #theta_{#gamma}", 50, 0, 3.14);
    TH1F histoPhoton1CandidateTheta ("Thetaof1stPhotonofCandidate", "Theta of 1st Photon of Candidate; #theta_{#gamma}", 180, 0, 180 );
    TH1F histoPhoton2CandidateTheta ("Thetaof2ndPhotonofCandidate", "Theta of 2nd Photon of Candidate; #theta_{#gamma}", 180, 0, 180 );
	TH1F histoPtof1stPhotonofCandidate ("Ptof1stPhotonofCandidate", "Ptof1stPhotonofCandidate; Pt (GeV)", 150, 0, 150);
	TH1F histoPtof2ndPhotonofCandidate ("Ptof2ndPhotonofCandidate", "Ptof2ndPhotonofCandidate; Pt (GeV)", 150, 0, 150);
	TH1F histoHighestPhotonPt ("HighestPhotonPt", "HighestPhotonPt; Pt (GeV)", 250, 0, 500);//najveći pT od fotona koju prodju Pt cut
	TH1F histo2ndHighestPhotonPt ("2ndHighestPhotonPt", "2ndHighestPhotonPt; Pt (GeV)", 250, 0, 500);//drugi najveći Pt od fotona koji prodju pt cut
	TH1F histo2ndHighestPhotonPtZoomed ("2ndHighestPhotonPtZoomed", "2ndHighestPhotonPtZoomed; Pt (GeV)", 50, 0, 50);//drugi najveći Pt od fotona koji prodju pt cut
	TH1F histoNeutrinoEnergy ("NeutrinoE", "NeutrinoE; E_{neutrino} (GeV)",300, 0, 3000);
	TH1F histoVisibleEnergy ("Visible_Energy", "Visible_Energy; E_{vis} (GeV)",300, 0, 3000);
	TH1F CandidatePhotonEnergy2 ("EnergyPhoton2", "EnergyPhoton2; E_{#gamma2}", 80, 0, 800);
	TH1F CandidatePhotonEnergy1 ("EnergyPhoton1", "EnergyPhoton1; E_{#gamma1}", 80, 0, 800);
	TH1F histoHiggsPhotonsPt ("higgsPhotonsPt", "higgsPhotonsPt; Pt (GeV)", 1500, 0, 1500);//drugi najveći Pt od fotona koji prodju pt cut
	//TH1F histoDaughtersPt	("higgsPhotonsPt", "higgsPhotonsPt; Pt (GeV)", 150, 0, 150);
	TH1F ptOfPhotons ("ptOfPhotons", "ptOfPhotons; Pt_{#gamma} (GeV)", 1500, 0, 1500);
	TH1F ptOf2ndHiggsPhoton ("ptOf2ndHiggsPhoton", "ptOf2ndHiggsPhoton; p_{T#gamma} (GeV)", 1500, 0, 1500);

	TH1F histohiggsphotonstheta ("higgsphotonstheta", "ThetaofhiggsPhoton; #theta_{#gamma}", 180, 0, 180);


	TTree eventList("eventsSignal", "ILD event list");
	varListGoran vl; /* SL specific */
	eventList.Branch("invM", &(vl.CandidateM), "invM");
	eventList.Branch("CanE", &(vl.CandidateEnergy), "CanE");
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
	eventList.Branch("pTCandidatePhotons2", &(vl.pTCandidatePhotons2), "pTCandidatePhotons2");
	eventList.Branch("ptHiggsPhotons", &(vl.ptHiggsPhotons), "ptHiggsPhotons");
	eventList.Branch("thetaHiggsPhotons", &(vl.thetaHiggsPhotons), "thetaHiggsPhotons");
	eventList.Branch("energyHiggsPhotons", &(vl.energyHiggsPhotons), "energyHiggsPhotons");
	eventList.Branch("pTHiggsPhoton1", &(vl.pTHiggsPhoton1), "pTHiggsPhoton1");
	eventList.Branch("pTHiggsPhoton2", &(vl.pTHiggsPhoton2), "pTHiggsPhoton2");
	eventList.Branch("energyHiggsPhoton1", &(vl.energyHiggsPhoton1), "energyHiggsPhoton1");
	eventList.Branch("energyHiggsPhoton2", &(vl.energyHiggsPhoton2), "energyHiggsPhoton2");
	eventList.Branch("thetaHiggsPhoton1", &(vl.thetaHiggsPhoton1), "thetaHiggsPhoton1");
	eventList.Branch("thetaHiggsPhoton2", &(vl.thetaHiggsPhoton2), "thetaHiggsPhoton2");
	eventList.Branch("HiggsMass", &(vl.HiggsMass), "HiggsMass");

	Float_t recoPhoton1, recoPhoton2;

	TTree photonsReco("photonsReco", "photonsReco");
	varListPhotonsReco vp; /* SL specific */
	photonsReco.Branch("energyPhotonsReco", &(vp.energyPhotonsReco), "energyPhotonsReco");
	photonsReco.Branch("ptPhotonsReco", &(vp.ptPhotonsReco), "ptPhotonsReco");
	photonsReco.Branch("thetaPhotonsReco", &(vp.thetaPhotonsReco), "thetaPhotonsReco");

	TTree recoP("recoP", "recoP");

	recoP.Branch("recoPhoton1", &recoPhoton1, "recoPhoton1");
	recoP.Branch("recoPhoton2", &recoPhoton2, "recoPhoton2");


	Float_t pt_restPhotons;
	Float_t energy_restPhotons;
	TTree restPhotons("restPhotons", "restPhotons");
	restPhotons.Branch("pt_restPhotons", &pt_restPhotons, "pt_restPhotons");
	restPhotons.Branch("energy_restPhotons", &energy_restPhotons, "energy_restPhotons");




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
			vector<TLorentzVector> daughters;
			TLorentzVector higgsPhoton1;
			TLorentzVector higgsPhoton2;//

			TLorentzVector ptRecHiggsPhotons1;
			TLorentzVector ptRecHiggsPhotons2;

			vector<Double_t> ptOfHiggsPhotons;
			vector<TLorentzVector> ptRecHiggsPhotons;

			vector<TLorentzVector> ptRestPhotons;
			TLorentzVector ptRecHiggsPhotons3;
			TLorentzVector ptRecHiggsPhotons4;



			vector<Double_t> thetaOfHiggsPhotons;
			vector<TLorentzVector> higgsPhotons;


			vector <Double_t> PtPhotons;


				std::vector<std::string> colNames = *evt->getCollectionNames();
			/*std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < colNames.size(); ++i)
			{
			std::cout << colNames[i] << endl;
			}*/
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticlesSkimmed");/*("PandoraPFANewPFOs");("PandoraPFOCollection");*/
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*)evt->getCollection("PandoraPFANewPFOs");/*("PandoraPFOCollection");*/


			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)

			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(recParticle->getType());

				if (particlePDG !=22) continue;

				ptRestPhotons.push_back(temp);

			}

			double_t Evis = 0;	//Ukupna energija po dogadjaju
			double_t Eneutrino = 0;

			bool leptonFound = false;	//uslov da nemamo leptone
			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

				if (mcParticle->getGeneratorStatus() != 1) continue;

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici


				const double *p = mcParticle->getMomentum(); // impuls čestice
				double e = mcParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(mcParticle->getPDG());
			//	cout<<"particle PDG is: "<<particlePDG<<endl;


				if(particlePDG == 22)	//rad sa fotonima (PDG=22)
				{
					//cout << "da li je foton: "<<particlePDG<<endl;
					PtPhotons.push_back(temp.Pt());


					Double_t theta = temp.Theta()*180/M_PI; //promenljiva koja nam daje Theta čestice
					//Double_t phi = temp.Phi();		//promenljiva koja nam daje Phi
					Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice
					ptOfPhotons.Fill(temp.Pt());

					const EVENT::MCParticleVec & parents = mcParticle->getParents();
					if (parents.size() > 0)
					{
						IMPL::MCParticleImpl* parent = (IMPL::MCParticleImpl*)parents[0];

					/*	if (parent->getPDG() != 25){

							EVENT::LCRelation* link_to_rec_lepton1 = 0;
							EVENT::MCParticle* pointer_to_mc_lepton1 = 0;
							EVENT::ReconstructedParticle* pointer_to_rec_lepton1 = 0;
							EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");
							pointer_to_mc_lepton1 = mcParticle;

							for (int j = 0; j < links -> getNumberOfElements(); j++){
							EVENT::LCRelation* linki = (EVENT::LCRelation*) links -> getElementAt(j);
							EVENT::MCParticle* mcpj = (EVENT::MCParticle*) linki -> getTo();
							link_to_rec_lepton1 = 0;

								if (pointer_to_mc_lepton1 == mcpj)
								{
									link_to_rec_lepton1 = linki;
								}

								if (link_to_rec_lepton1){
									EVENT::ReconstructedParticle* rpi1 = (EVENT::ReconstructedParticle*) link_to_rec_lepton1 -> getFrom();
									pointer_to_rec_lepton1 = rpi1;


									if (rpi1 -> getType() == 22 ){
									TLorentzVector tempPhotons; //četvorovektor u koji sakupljamo informacije o svakoj čestici
									ptRestPhotons.push_back(tempPhotons);
									}

						}
						}
						}*/


						if (parent->getPDG() == 25)
						{
							ptOfHiggsPhotons.push_back(Pt);
							histohiggsphotonstheta.Fill(temp.Theta()*180/M_PI);
							histoHiggsPhotonsPt.Fill(Pt);
							vl.ptHiggsPhotons=temp.Pt();
							vl.energyHiggsPhotons=temp.E();
							vl.thetaHiggsPhotons = temp.Theta();
							higgsPhotons.push_back(temp);

		 					EVENT::LCRelation* link_to_rec_lepton1 = 0;
		 					EVENT::MCParticle* pointer_to_mc_lepton1 = 0;
		 					EVENT::ReconstructedParticle* pointer_to_rec_lepton1 = 0;
		 					EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");
		 					pointer_to_mc_lepton1 = mcParticle;

		 					for (int j = 0; j < links -> getNumberOfElements(); j++){
							EVENT::LCRelation* linki = (EVENT::LCRelation*) links -> getElementAt(j);
							EVENT::MCParticle* mcpj = (EVENT::MCParticle*) linki -> getTo();
							link_to_rec_lepton1 = 0;

								if (pointer_to_mc_lepton1 == mcpj)
								{
									link_to_rec_lepton1 = linki;
								}

								if (link_to_rec_lepton1){
									EVENT::ReconstructedParticle* rpi1 = (EVENT::ReconstructedParticle*) link_to_rec_lepton1 -> getFrom();
									pointer_to_rec_lepton1 = rpi1;


									if (rpi1 -> getType() == 22 ){
										TLorentzVector tempPhotons; //četvorovektor u koji sakupljamo informacije o svakoj čestici


										const double *pp = mcParticle->getMomentum(); // impuls čestice
										double ep = mcParticle->getEnergy();	//energija čestice
										tempPhotons.SetPxPyPzE(pp[0], pp[1], pp[2], ep);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
										ptRecHiggsPhotons.push_back(tempPhotons);

										vp.energyPhotonsReco = ep;
										vp.ptPhotonsReco = tempPhotons.Pt();
										vp.thetaPhotonsReco = tempPhotons.Theta()*180/M_PI;
										TLorentzVector currentConeAxis = tempPhotons; //foton oko koga pravimo konus
										Double_t coneEnergy = 0; //Energija konusa oko fotona


										for (Int_t k = 0; k < mcParticles->getNumberOfElements() ; k++) /*&& coneEnergy <= maxConeEnergy -zaustavlja fot loop kada predje maxE */
										{
											IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(k);
											if (mcParticle->getGeneratorStatus() != 1) continue;
											const double *impuls = mcParticle->getMomentum();
											double energija = mcParticle->getEnergy();
											TLorentzVector otherParticle;
											otherParticle.SetPxPyPzE(impuls[0], impuls[1], impuls[2],energija);//četvorovektor drugih čestica
											if (currentConeAxis == otherParticle) continue; //da ne bi ubrojali foton oko kog pravimo konus
											Double_t particleAngle = currentConeAxis.Angle(otherParticle.Vect())* 180 / TMath::Pi() ;//ugao izmedju fotona oko kog pravimo konus i čestice
											//anglePhotonParticle.Fill(particleAngle);//histogram uglova izmedju fotona i čestice
											//anglePhotonParticlewE.Fill(particleAngle, otherParticle.E());//histogram uglova izmedju fotona i čestice otežinjen energijom te čestice
										//	cout << "test ostale čestice: "<<particleAngle<<endl;


											//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
											if (particleAngle <=coneAngle)
											{
												coneEnergy += otherParticle.E();
												histoconeEnergy.Fill(coneEnergy); //histogram koji iscrtava energiju konusa
											//	cout << "test coneAngle: "<<coneEnergy<<endl;

											}

										}// end of for loop for cone
										vp.coneEnergyPhotonsReco = coneEnergy;
										photonsReco.Fill();


									}//end of if po rekontruisanoj cestici koja je foton

								}//end of link to rec lepton1


		 					}//end of link loop


						}
					}

				//	cout << "theta i pt: "<<theta<<"  "<<Pt<<endl;
					//cout << "Photons Size  : "<<PtPhotons.size()<<endl;



					//uzimamo u obzir samo one fotone koji nam prodju uslove
					if(temp.Pt() > minPt1    )
					{
					//	cout << "test1: "<< theta<<endl;

						//PtPhotons.push_back(temp.Pt());

						histoPhotonEnergy.Fill(temp.E());

						counterPhotonsPt1++;

						allPhotons.push_back(temp);//sakupili smo sve fotone
						TLorentzVector currentConeAxis = temp; //foton oko koga pravimo konus
						Double_t coneEnergy = 0; //Energija konusa oko fotona




						for (Int_t k = 0; k < mcParticles->getNumberOfElements() ; k++) /*&& coneEnergy <= maxConeEnergy -zaustavlja fot loop kada predje maxE */
						{
							IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(k);
							if (mcParticle->getGeneratorStatus() != 1) continue;
							const double *impuls = mcParticle->getMomentum();
							double energija = mcParticle->getEnergy();
							TLorentzVector otherParticle;
							otherParticle.SetPxPyPzE(impuls[0], impuls[1], impuls[2],energija);//četvorovektor drugih čestica
							if (currentConeAxis == otherParticle) continue; //da ne bi ubrojali foton oko kog pravimo konus
							Double_t particleAngle = currentConeAxis.Angle(otherParticle.Vect())* 180 / TMath::Pi() ;//ugao izmedju fotona oko kog pravimo konus i čestice
							//anglePhotonParticle.Fill(particleAngle);//histogram uglova izmedju fotona i čestice
							//anglePhotonParticlewE.Fill(particleAngle, otherParticle.E());//histogram uglova izmedju fotona i čestice otežinjen energijom te čestice
						//	cout << "test ostale čestice: "<<particleAngle<<endl;


							//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
							if (particleAngle <=coneAngle)
							{
								coneEnergy += otherParticle.E();
								histoconeEnergy.Fill(coneEnergy); //histogram koji iscrtava energiju konusa
							//	cout << "test coneAngle: "<<coneEnergy<<endl;

							}


						}
						histoconeEnergyFilter.Fill(coneEnergy);//histogram koji iscrtava energiju konusa

						//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
						if (coneEnergy <= maxConeEnergy)
							{
								counterAllCuts++;//brojač za fotone koji prodju sve cutove
								photons.push_back(temp);
							//	cout << "test coneEnergy: "<<endl;


							}

					}//end of minPt loop

				}//end of photon loop
				//dogadjaj ne prolazi ako se u dogadjaju detektuje lepton
				else if (particlePDG != 22)
					{

					if(temp.Pt() > leptonPt)
						leptonFound = true;
					histoPtOtherParticles.Fill(temp.Pt());
					}
				if (!(particlePDG ==12 || particlePDG ==14 || particlePDG ==16) && temp.Pt()>5)
				{
					Evis += e ;//ukupna enrgija po dogadjaju
				}

//				Evis += e ;//ukupna enrgija po dogadjaju


				if ((particlePDG ==12 || particlePDG ==14 || particlePDG ==16))
				{
					Eneutrino += e ;//ukupna enrgija po dogadjaju
					histoNeutrinoEnergy.Fill(Eneutrino);
				}


			}  // end of particle loop
			sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());
			histo2ndHighestPhotonPt.Fill(PtPhotons[1]);
			histo2ndHighestPhotonPtZoomed.Fill(PtPhotons[1]);

			sort (ptOfHiggsPhotons.begin(), ptOfHiggsPhotons.end(), greater<int>() );
			ptOf2ndHiggsPhoton.Fill(ptOfHiggsPhotons[1]);

			if (higgsPhotons[0].Pt()>higgsPhotons[1].Pt())
			{
				higgsPhoton1 = higgsPhotons[0];
				higgsPhoton2 = higgsPhotons[1];

			}else
				{
					higgsPhoton2 = higgsPhotons[0];
					higgsPhoton1 = higgsPhotons[1];


				}//
			if(ptRecHiggsPhotons.size()!=2) continue;
			if (ptRecHiggsPhotons[0].Pt()>ptRecHiggsPhotons[1].Pt())
			{
				ptRecHiggsPhotons1 = ptRecHiggsPhotons[0];
				ptRecHiggsPhotons2 = ptRecHiggsPhotons[1];

			}else
				{
				ptRecHiggsPhotons2 = higgsPhotons[0];
				ptRecHiggsPhotons1 = higgsPhotons[1];


				}
			ptRecHiggsPhotons3 =  ptRecHiggsPhotons1;
			ptRecHiggsPhotons4= ptRecHiggsPhotons2;

			recoPhoton1 = ptRecHiggsPhotons3.Pt();
			recoPhoton2 = ptRecHiggsPhotons4.Pt();
			recoP.Fill();

			vl.pTHiggsPhoton2 = higgsPhoton2.Pt();
			vl.pTHiggsPhoton1 = higgsPhoton1.Pt();

			vl.energyHiggsPhoton1 = higgsPhoton1.E();
			vl.energyHiggsPhoton2 = higgsPhoton2.E();
			vl.thetaHiggsPhoton1 = higgsPhoton1.Theta();
			vl.thetaHiggsPhoton2 = higgsPhoton2.Theta();

			eventList.Fill();






			Int_t Emiss = EnergyCenterMass - Evis; // missing energy
			histoMissingEnergy.Fill(Emiss);

			//cout << "lepton found: " <<leptonFound <<endl;
			if(leptonFound) continue; //ako je u finalnom stanju imamo lepton ili kvark (sem neutrina), događaj se preskače


			// At this point, all relevant photons have been collected.

			//cout << "photon size is: " << photons.size()<<endl;
			//proveravamo da li u eventu ima detektovano više od jednog fotona koji bi bili kandidat za higsov bozon
			if (photons.size()  == 2)
			{
				vector<CandidateData> candidates;
				histoHighestPhotonPt.Fill(PtPhotons[0]);

				histoNumberPhotonsbyEvent.Fill(photons.size());//histogram koji nam pokazuje broj fotona po događaju
				//histoHiggsPhotonsPt.Fill(PtPhotons[0]);
				//histoHiggsPhotonsPt.Fill(PtPhotons[1]);


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
					TLorentzVector higgs = candidate->Higgs();
					Double_t Ehiggs = higgs.Energy();//energija kandidata
					Double_t Eremaining = Evis - Ehiggs;//energija koja ostaje posle higsovog kandidata


					TVector3 boosttoparent = -(higgs.BoostVector());//prelazimo iz sistema CM u LAB sistem
					TLorentzVector photon = candidate->Photon1;
					TLorentzVector photon2 = candidate->Photon2;


					TLorentzVector higherPtphoton ;
					TLorentzVector lowerPtphoton ;
					if (photon.Pt() > photon2.Pt())
					{
						higherPtphoton = photon;
						lowerPtphoton = photon2;
					}else {
						higherPtphoton = photon2;
						lowerPtphoton = photon;
					}

				    Double_t theta1 = photon.Theta()*180/M_PI;
				    Double_t theta2 = photon2.Theta()*180/M_PI;




					histoPhoton1CandidateTheta.Fill(theta1);
					histoPhoton2CandidateTheta.Fill(theta2);



				    vl.CandidateM=higgs.M();
				    vl.CandidateEnergy = higgs.E();
				    vl.CandidatePhi = higgs.Phi();
				    vl.CandidatePt = higgs.Pt();
				    vl.CandidateTheta = higgs.Theta();
				    vl.pLCandidate = higgs.Z();
				    vl.thetaCandidatePhotons1 = higherPtphoton.Theta();
				    vl.thetaCandidatePhotons2 = lowerPtphoton.Theta();
				    vl.pLCandidatePhotons1 = higherPtphoton.Z();
				    vl.pLCandidatePhotons2 = lowerPtphoton.Z();
				    vl.pTCandidatePhotons1 = higherPtphoton.Pt();
				    vl.pTCandidatePhotons2 = lowerPtphoton.Pt();




					Double_t Px11 = candidate->Photon1.X();
					Double_t Py11 = candidate->Photon1.Y();
					Double_t Px22 = candidate->Photon2.X();
					Double_t Py22 = candidate->Photon2.Y();
					Double_t higherPtPhoton = TMath::Max(TMath::Sqrt(Px11*Px11 +Py11*Py11), TMath::Sqrt(Px22*Px22 +Py22*Py22));//Energija energičnijeg fotona koji je kandidat
					Double_t lowerPtPhoton = TMath::Min(TMath::Sqrt(Px11*Px11 +Py11*Py11), TMath::Sqrt(Px22*Px22 +Py22*Py22));//Energija manje energičnog fotona koji je kandidat

					CandidatePhotonEnergy1.Fill(candidate->Photon1.E());
					CandidatePhotonEnergy2.Fill(candidate->Photon2.E());

					histoPtof1stPhotonofCandidate.Fill(higherPtPhoton);
					histoPtof2ndPhotonofCandidate.Fill(lowerPtPhoton);



					//higgs.Boost(boosttoparent);
					photon.Boost(boosttoparent);//prebacujemo fotone u sistem CM
					photon2.Boost(boosttoparent);

					numCutEvents++;//brojimo dogadjaje koji prodju uslove za energiju, invarijantnu masu, uglove,Pt

					//TVector3 higgs3unit  = higgs.Vect();//.Unit();
					//TVector3 photon3unit  = photon.Vect();//.Unit();

					//Double_t numerator = photon.BoostVector() * higgs.BoostVector();
					//Double_t denominator = (photon.BoostVector().Mag())*(higgs.BoostVector().Mag());
					//Double_t temp = numerator/denominator;//računa kosinus izmedju dva fotona

					//Double_t helicity =  TMath::Abs( temp/* higgs3unit.Dot(photon3unit) */) ;
					//Double_t arcCos = TMath::ACos(helicity);

					Double_t higherEnergyPhotonE = TMath::Max(candidate->Photon1.E(), candidate->Photon2.E());//Energija energičnijeg fotona koji je kandidat
					Double_t lowerEnergyPhotonE = TMath::Min(candidate->Photon1.E(), candidate->Photon2.E());//Energija manje energičnog fotona koji je kandidat

					Double_t Px1 = candidate->Photon1.X();
					Double_t Py1 = candidate->Photon1.Y();
					Double_t Px2 = candidate->Photon2.X();
					Double_t Py2 = candidate->Photon2.Y();
					Double_t zbirPt = TMath::Sqrt(Px1*Px1 +Py1*Py1) + TMath::Sqrt(Px2*Px2 +Py2*Py2); // zbir transverzalnih impulsa
					Double_t Theta1 = candidate->Photon1.Theta()*180/M_PI;
					Double_t Theta2 = candidate->Photon2.Theta()*180/M_PI;
					Double_t Phi1 = candidate->Photon1.Phi();
					Double_t Phi2 = candidate->Photon2.Phi();
					Double_t CMTheta2 =photon2.Theta()*180/M_PI;
					Double_t CMTheta1 =photon.Theta()*180/M_PI;

				//	Double_t angleBetweenPotons = candidate->Photon1.Angle(candidate->Photon2.Vect());
				//	eventList.Fill();

					//histoPhotonEnergy.Fill(higherEnergyPhotonE, lowerEnergyPhotonE);
					histoTheta.Fill(Theta1);
					histoTheta.Fill(Theta2);
					histoPhi.Fill(Phi1);
					histoPhi.Fill(Phi2);
					histoCandidateM.Fill(higgs.M());
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
				//	histohelicityAngle.Fill(temp);
				//	histoarcCos.Fill(helicity);
					histogram.Fill(CMTheta1, CMTheta2);
					histoHardPhotonsByEvent.Fill(hardPhotonsByEvent);
					histoVisibleEnergy.Fill(Evis);



				}
			}


			cout <<"veličina niza svih fotona je: "<<ptRestPhotons.size()<<endl;
		for (int p=0; p < (int)ptRestPhotons.size(); p++){
			if (ptRestPhotons[p] == ptRecHiggsPhotons3) continue;
			if (ptRestPhotons[p] == ptRecHiggsPhotons4) continue;
			pt_restPhotons = ptRestPhotons[p].Pt();
			energy_restPhotons = ptRestPhotons[p].E();
			restPhotons.Fill();

		}


//cout <<"iznad event loopa_________________________________________________________"<<endl;
		} // End of event loop


		totalCutEvents+=numCutEvents;//ukupan broj dogadjaja koji prodje cutove iz svih fajlova
		totalEvents+=brojDogadjaja; //ukupan broj dogadjaja iz svih fajlova

		cout << "Broj event-a posle cut-a je: " << numCutEvents << endl;
		cout << "Broj fotona posle primarnog cut-a  je: " << counterPhotonsPt1 << endl;
		cout << "Broj fotona posle svih cut-ova  je: " << counterAllCuts << endl;


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
	file.open("MCproba.txt", ios_base::out | ios_base::app);

	file << minPt1 << "\t" << percentagePt1 << endl;
	file << minPt2 << "\t" << percentagePt2 << endl;
	file << minPt3 << "\t" << percentagePt3 << endl;
	file << minPt4 << "\t" << percentagePt4 << endl;


	file.close();

	ofstream results ;
	results.open("MCParticlesSkimmed.txt", ios_base::out | ios_base::app);
	results << "pT > " << minPt1 <<", "<< minInvMass<<" < M < " << maxInvMass << ", leptonPt > "<< leptonPt << ", cone Angle < "<< 2*coneAngle << ", efikasnost = "<<percentagePt1 <<endl;
	results.close();


	TGraph efikasnostproba ("MCproba.txt", "%lg %lg", "\t");
	efikasnostproba.GetXaxis()->SetTitle("Pt_{#gamma} (GeV)");
	efikasnostproba.GetYaxis()->SetTitle("Efficiency");
	efikasnostproba.SetTitle("Signal");



	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	eventList.Write();
	photonsReco.Write();
	recoP.Write();
	restPhotons.Write();

	histoCanditateTheta.Write();
	histoCandidatePt.Write();
	histoCandidateM.Write();
	histoCandidateEnergy.Write();
	histoCandidatePhi.Write();
	histoBoost.Write();
	histoRemainingEnergy.Write();
//	histohelicityAngle.Write();
	histoHigherEnergyPhoton.Write();
	histoLowerEnergyPhoton.Write();
	histoZbirPt.Write();
	histoPhotonEnergy.Write();
	histoTheta.Write();
	histoPhi.Write();
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
	histoVisibleEnergy.Write();
	histo2ndHighestPhotonPtZoomed.Write();
	ptOfPhotons.Write();
	histoHiggsPhotonsPt.Write();
	ptOf2ndHiggsPhoton.Write();
	histohiggsphotonstheta.Write();




	rootFile.Close();

//	gROOT->ProcessLine(".x /home/Goran/Programs/crtanje_histograma/test/test/CLICdpStyle.C");

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

	TString rfName = "mcparticle.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
