/***********************************************************
 *
 * 	Read reconstructed photons from a  .slcio file, and draws histograms for different variables
 *
 *  Author:Goran Kačaravić
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

using namespace std;


const double mH = 126.0;               //Higgs mass
const Double_t coneAngle = 92.5;		   //ugao konusa
const Double_t maxConeEnergy = 200000 ;	   //maksimalna energija konusa
const Double_t minInvM = 0;
const Double_t maxInvM = 140000000;
const Double_t minPhotonAngle = 0;
const Double_t maxPhotonAngle = 181;
const Double_t minPt = 0;


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



	//histogrami za različite kinematičke varijable
	TH1F histoTheta ("histoTheta", "Theta; #theta_{#gamma}", 50, 0, 3.14);
	TH1F histoPhi ("histoPhi", "Phi; #phi_{#gamma}", 50, 0, 3.14);
	TH1F energyOfPhotons ("energyOfPhotons", "PhotonEnergy; E_{#gamma} (GeV)", 150, 0, 1500);
	TH1F histoCandidateM ("CandidateInvariantM", "CandidateInvariantM; M_{H} (GeV)", 100, 110, 140);
	TH1F histoCandidatePt ("CandidatePt", "CandidatePt; Pt (GeV)", 300, 0, 1550);
	TH1F histoMissingEnergy ("MissingEnergy", "Missing Energy; miss E (GeV)", 300, 0, 3000);
	TH1F histoCandidateEnergy ("CandidateE", "CandidateE; E_{H} (GeV)",300, 0, 1550);
	TH1F histoCanditateTheta ("CandidateTheta", "CandidateTheta; #theta_{H}", 90, 0, 180);
	TH1F histoCandidatePhi ("CandidatePhi", "CandidatePhi' #phi_{#gamma}", 50, 0, 3.14);
	TH1F histoBoost ("CandidateBoost", "CandidateBoost; #beta_{H}", 180, 0, 1);
	TH1F histoZbirPt ("ZbirPt", "ZbirPt; Pt_{1} + Pt_{2}", 300, 0, 3000);
	TH1F histoHigherEnergyPhoton ("HigherEnergyPhoton", "HigherEnergyPhoton; E_{#gamma1}", 150, 0, 1500);
	TH1F histoLowerEnergyPhoton ("LowerEnergyPhoton", "LowerEnergyPhoton; E_{#gamma2}", 150, 0, 1500);
	TH1F histohelicityAngle ("helicityAngle", "helicityAngle; cos#theta", 100, 0, 1);
	TH1F histoarcCos ("arcCos", "arcCOs; arcCos#theta", 300, 0, 1);
	TH2F histogram ("test","title", 20, 0, 3.14, 20, 0, 3.14);
	TH1F histoangleBetweenPhotons ("angleBetweenPhotons", "angleBetweenPhotons; #alpha", 50, 0, 3.14);
	TH1F histoNumberPhotonsbyEvent ("NumberofPhotonsPerEvent", "NumberofPhotonsPerEvent", 5, 0, 5);
	TH1F anglePhotonParticle ("AngleBetweenphotonandparticle", "anglephotonparticle", 50, 0, 20);
	TH1F histoconeEnergy ("ConeEnergy", "ConeEnergy; E (GeV)", 100, 0, 100);
	TH1F anglePhotonParticlewE ("Angle Between photon and particle wE", " photon particle angle wE", 50, 0, 20);//otezinjeno sa Energijom
	TH1F histoconeEnergyFilter ("ConeEnergy", "ConeEnergy; E (GeV)", 100, 0, 10);
	TH1F histoHardPhotonsByEvent ("HardPhotonsbyEvent", "HardPhotonsbyEvent", 10, 0, 10);
	TH1F histoPtMax ("HihestPhotonPt", "MaxPhotonPt; Pt (GeV)", 300, 0, 1500);
	TH1F histoPt2Max ("2ndHighestPhotonPt", "2ndHighestPhotonPt; Pt (GeV)", 300, 0, 1500);
	TH1F histoUnknownParticle ("XInvariantM", "XInvariantM; M_{X} (GeV)", 100, 0, 10000);
	TH1F histoPt ("PhotonPt", "PhotonPt; Pt (GeV)", 300, 0, 1500);
	TH1F ptOfHiggsPhotons ("ptOfHiggsPhotons", "ptOfHiggsPhotons; Pt_{#gamma} (GeV)", 300, 0, 1500);
	TH1F pLOfHiggsPhotons ("pLOfHiggsPhotons", "pLOfHiggsPhotons; Pz_{#gamma} (GeV)", 300, 0, 1500);
// pt svih fotona	TH1F ptOfPhotons ("ptOfPhotons", "ptOfPhotons; Pt_{#gamma} (GeV)", 250, 0, 500);
	TH1F pLOfPhotons ("pLOfPhotons", "pLOfPhotons; Pz_{#gamma} (GeV)", 250, 0, 500);
	TH1F PtofCandidatePhotons ("PtofCandidatePhotons", "PtofCandidatePhotons; Pt_{#gamma} (GeV)", 100, 0, 1500);
	TH1F energyofCandidatePhotons ("energyofCandidatePhotons", "energyofCandidatePhotons; E_{#gamma} (GeV)", 100, 0, 1500);
	TH1F thetaOfCandidatePhotons ("thetaOfCandidatePhotons", "thetaOfCandidatePhotons; #theta_{#gamma}", 90, 0, 180);
	TH1F ptOfPhotons ("ptOfPhotons", "ptOfPhotons; Pt_{#gamma} (GeV)", 250, 0, 500);// pt higgsovih fotona





	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	//petlja koja iščitava .slcio podatke ukoliko ima više fajlova za jedan process
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{

		/*cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);
				try
				{
					lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));
				}
				catch(lcio::IOException &ex)
				{
					cout << ". Could not open.\n"; // Exception " << ex.what() << endl;
					continue;
				}*/

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
		cout << ". Reading.\n";

		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja
		cout << "Broj dogadjaja je  : " << brojDogadjaja << endl;
		Int_t numCutEvents = 0; 	//broj dogadjaja posle određenih cut-ova

		std::vector <TLorentzVector> allPhotons; //vektor koji prikuplja sve fotone
		Int_t brojac =0;
		Int_t counterPrimaryCuts = 0;
		Int_t counterAllCuts = 0;
		Int_t hardPhotons = 0;
		//Int_t counterX = 0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
			Int_t hardPhotonsByEvent = 0;
			brojac++;
			vector<TLorentzVector> photons;
			vector<TLorentzVector> particles;   // other than photons
			vector <Double_t> PtPhotons;

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");//("MCParticle");

			double_t Etot = 0;	//Ukupna energija po dogadjaju

			bool leptonFound = false;	//uslov da nemamo leptone
			//Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
			for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
				/*const double *p = mcParticle->getMomentum(); // impuls čestice
				double e = mcParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				if (mcParticle->getPDG() == 94)
					{
						histoUnknownParticle.Fill(temp.M());
						counterX++;
					}
*/
				// Interesuju nas samo finalni fotoni
				if (mcParticle->getGeneratorStatus() != 1) continue;

				const double *p = mcParticle->getMomentum(); // impuls čestice
				double e = mcParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

				Int_t particlePDG = fabs(mcParticle->getPDG());

				if(particlePDG == 22)	//rad sa fotonima (PDG=22)
				{


					Double_t theta = temp.Theta(); //promenljiva koja nam daje Theta čestice
					//Double_t phi = temp.Phi();		//promenljiva koja nam daje Phi
					double Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice

					//uzimamo u obzir samo one fotone koji nam prodju uslove
					if(theta > minPhotonAngle * TMath::Pi() / 180  && theta < maxPhotonAngle  * TMath::Pi() / 180  && temp.Pt() > 0)
					{

						TLorentzVector temp1; //četvorovektor u koji sakupljamo informacije o svakoj čestici
						const double *p1 = mcParticle->getMomentum(); // impuls čestice
						double e1 = mcParticle->getEnergy();	//energija čestice
						temp1.SetPxPyPzE(p1[0], p1[1], p1[2], e1);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

						counterPrimaryCuts++;
						allPhotons.push_back(temp);//sakupili smo sve fotone
						TLorentzVector currentConeAxis = temp; //foton oko koga pravimo konus
						Double_t coneEnergy = 0; //Energija konusa oko fotona

						PtPhotons.push_back(temp.Pt());
						sort(PtPhotons.begin(), PtPhotons.end(), greater<int>());
						energyOfPhotons.Fill(temp.E());

						histoPtMax.Fill(PtPhotons[0]);
						histoPt2Max.Fill(PtPhotons[1]);
						histoPt.Fill(temp1.Pt());
					//	ptOfPhotons.Fill(Pt);
						pLOfPhotons.Fill(temp.Pz());


						for (Int_t l = 0; l < (int)mcParticle->getParents().size() ; l++)
						{
							IMPL::MCParticleImpl* parent = (IMPL::MCParticleImpl*) mcParticle->getParents()[l];
							double const *parentTempP = parent->getMomentum();//impuls roditelja
							double parentTempE = parent->getEnergy();//energija roditelja

							TLorentzVector parentP;//četvorovektor roditelja fotona
							parentP.SetPxPyPzE(parentTempP[0], parentTempP[1], parentTempP[2], parentTempE);


							/*Int_t parentPDG = fabs(parent->getParent(l)->getPDG());//proverava koja čestica je roditelj
							if (parentPDG == abs(25))//proverava da li je Higgs
							{
								ptOfHiggsPhotons.Fill(Pt);
								pLOfHiggsPhotons.Fill(temp.Pz());

							}*/
						}


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
							anglePhotonParticle.Fill(particleAngle);//histogram uglova izmedju fotona i čestice
							anglePhotonParticlewE.Fill(particleAngle, otherParticle.E());//histogram uglova izmedju fotona i čestice otežinjen energijom te čestice


							//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
							if (particleAngle <=coneAngle)
							{
								coneEnergy += otherParticle.E();
								histoconeEnergy.Fill(coneEnergy); //histogram koji iscrtava energiju konusa
							}


						}
						//uslov u kom proveravamo da li je foton izolovan ili je deo jet-a
						if (coneEnergy <= maxConeEnergy)
							{
								counterAllCuts++;//brojač za fotone koji prodju sve cutove
								photons.push_back(temp);
								histoconeEnergyFilter.Fill(coneEnergy);//histogram koji iscrtava energiju konusa

								histoNumberPhotonsbyEvent.Fill(photons.size());//histogram koji nam pokazuje broj fotona po događaju

								//petlja u kojoj proveravamo da li je elektron ili pozitron iz snopa parent fotona (indrirektno preko Pt)



							}

					}

				}


				//dogadjaj ne prolazi ako se u dogadjaju detektuje lepton
				else if ((1 <= particlePDG && particlePDG <= 15) && !(particlePDG == 12 || particlePDG == 14))
				{
					// napraviti bez veta na neutrina!
					if(8 * TMath::Pi() / 180 <= temp.Theta() && temp.Theta() <= 172 * TMath::Pi() / 180)
						leptonFound = true;
				}
				Etot += e ;//ukupna enrgija po dogadjaju

			}  // end of particle loop

			if(leptonFound) continue; //ako je u finalnom stanju imamo lepton ili kvark (sem neutrina), događaj se preskače

			// At this point, all relevant photons have been collected.
			vector<TLorentzVector>::iterator it;
			for (it = photons.begin(); it < photons.end(); ++it)
			{
				//histoTheta.Fill(it->Theta());
				//histoPhi->Fill(it->Phi());
				//energyOfPhotons->Fill(it->Energy());

			}

			//proveravamo da li u eventu ima detektovano više od jednog fotona koji bi bili kandidat za higsov bozon
			if (photons.size() > 1)
			{
				vector<CandidateData> candidates;

				for (int i = 0; i < (int)photons.size() - 1; ++i)
				{
					for (int j = i + 1; j < photons.size(); ++j)
					{
						CandidateData current(photons[i], photons[j]);//uzima dva fotona

						TLorentzVector pair = current.Higgs();//sabira četvorovektore

						if (minInvM < pair.M() && pair.M() < maxInvM)//uslov da bi par bio higs kandidat
						{
							candidates.push_back(current);
						}
					}
				}

				//proverava koji kandidat ima najbližu vrednost higsovom bozonu
				vector<CandidateData>::iterator candidate = min_element(candidates.begin(), candidates.end(), IsBetterHiggsCandidate);
				if (candidate != candidates.end())

				{
					TLorentzVector higgs = candidate->Higgs();
					Double_t Ehiggs = higgs.Energy();//energija kandidata
					Double_t Emiss = Etot - Ehiggs;//energija neutrina


					TVector3 boosttoparent = -(higgs.BoostVector());//prelazimo iz sistema CM u LAB sistem
					TLorentzVector photon = candidate->Photon1;
					TLorentzVector photon2 = candidate->Photon2;

					PtofCandidatePhotons.Fill(photon.Pt());
					PtofCandidatePhotons.Fill(photon2.Pt());

					ptOfPhotons.Fill(photon.Pt());
					ptOfPhotons.Fill(photon2.Pt());


					energyofCandidatePhotons.Fill(photon.E());
					energyofCandidatePhotons.Fill(photon2.E());

					thetaOfCandidatePhotons.Fill(photon.Theta()*180/M_PI);
					thetaOfCandidatePhotons.Fill(photon2.Theta()*180/M_PI);
					//higgs.Boost(boosttoparent);
					photon.Boost(boosttoparent);//prebacujemo fotone u sistem CM
					photon2.Boost(boosttoparent);

					numCutEvents++;//brojimo dogadjaje koji prodju uslove za energiju, invarijantnu masu, uglove,Pt

					//TVector3 higgs3unit  = higgs.Vect();//.Unit();
					//TVector3 photon3unit  = photon.Vect();//.Unit();

					/*Double_t numerator = photon.BoostVector() * higgs.BoostVector();
					Double_t denominator = (photon.BoostVector().Mag())*(higgs.BoostVector().Mag());
					Double_t temp = numerator/denominator;//računa kosinus izmedju dva fotona

					Double_t helicity =  TMath::Abs( temp higgs3unit.Dot(photon3unit) ) ;
					Double_t arcCos = TMath::ACos(helicity);*/

					Double_t higherEnergyPhotonE = TMath::Max(candidate->Photon1.E(), candidate->Photon2.E());//Energija energičnijeg fotona koji je kandidat
					Double_t lowerEnergyPhotonE = TMath::Min(candidate->Photon1.E(), candidate->Photon2.E());//Energija manje energičnog fotona koji je kandidat

					Double_t Px1 = candidate->Photon1.X();
					Double_t Py1 = candidate->Photon1.Y();
					Double_t Px2 = candidate->Photon2.X();
					Double_t Py2 = candidate->Photon2.Y();
					Double_t zbirPt = TMath::Sqrt(Px1*Px1 +Py1*Py1) + TMath::Sqrt(Px2*Px2 +Py2*Py2); // zbir transverzalnih impulsa
					Double_t Theta1 = candidate->Photon1.Theta();
					Double_t Theta2 = candidate->Photon2.Theta();
					Double_t Phi1 = candidate->Photon1.Phi();
					Double_t Phi2 = candidate->Photon2.Phi();
					Double_t CMTheta2 =photon2.Theta();
					Double_t CMTheta1 =photon.Theta();

					Double_t angleBetweenPotons = candidate->Photon1.Angle(candidate->Photon2.Vect());

				//	energyOfPhotons.Fill(higherEnergyPhotonE, lowerEnergyPhotonE);
					histoTheta.Fill(Theta1);
					histoTheta.Fill(Theta2);
					histoPhi.Fill(Phi1);
					histoPhi.Fill(Phi2);
					histoCandidateM.Fill(higgs.M());
					histoCandidatePt.Fill(higgs.Pt());
					histoCandidateEnergy.Fill(higgs.Energy());
					histoCanditateTheta.Fill(higgs.Theta()*180/M_PI);
					histoCandidatePhi.Fill(higgs.Phi());
					histoMissingEnergy.Fill(Emiss) ;
					histoBoost.Fill(higgs.BoostVector().Mag());
					histoHigherEnergyPhoton.Fill(higherEnergyPhotonE);
					histoLowerEnergyPhoton.Fill(lowerEnergyPhotonE);
					histoZbirPt.Fill(zbirPt);
				//	histohelicityAngle.Fill(temp);
				//	histoarcCos.Fill(helicity);
					histogram.Fill(CMTheta1, CMTheta2);
					histoHardPhotonsByEvent.Fill(hardPhotonsByEvent);



					//}

				}
			}

		} // End of event loop
		cout << "Broj event-a posle cut-a je: " << numCutEvents << endl;
		cout << "Broj fotona posle primarnog cut-a  je: " << counterPrimaryCuts << endl;
		cout << "Broj fotona posle svih cut-ova  je: " << counterAllCuts << endl;
	//	cout << "Broj čestica PDG = 94 je: " << counterX << endl;
		cout << "Broj dogadjaja: " << brojDogadjaja << endl;


		Double_t percentage = ((Double_t) numCutEvents/brojDogadjaja)*100;
		cout << "procenat dogadjaja koji prodju cutove je: " << percentage << endl;


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

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	//eventList.Write();

	histoCanditateTheta.Write();
	histoCandidatePt.Write();
	histoCandidateM.Write();
	histoCandidateEnergy.Write();
	histoCandidatePhi.Write();
	histoBoost.Write();
	histoMissingEnergy.Write();
	histohelicityAngle.Write();
	histoHigherEnergyPhoton.Write();
	histoLowerEnergyPhoton.Write();
	histoZbirPt.Write();
	energyOfPhotons.Write();
	histoTheta.Write();
	histoPhi.Write();
	histoangleBetweenPhotons.Write();
	histoPt2Max.Write();
	histoPt.Write();
	histoPtMax.Write();
	ptOfHiggsPhotons.Write();
	pLOfHiggsPhotons.Write();
	ptOfPhotons.Write();
	pLOfPhotons.Write();
	PtofCandidatePhotons.Write();
	energyofCandidatePhotons.Write();
	thetaOfCandidatePhotons.Write();


	rootFile.Close();


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
	histoMissingEnergy.Draw();
	c4.Print("histoMissingEnergy.pdf");

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
	gStyle->SetOptStat(111111);
	c9.SetCanvasSize(1000,650);
	c9.Divide(1,1,0.01,0.01);
	c9.cd(1);
	//c9.GetPad(1)->SetLogy();
	histoHigherEnergyPhoton.Draw();
	c9.Print("histoHigherEnergyPhoton.pdf");

	TCanvas c10;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c10.SetCanvasSize(1000,650);
	c10.Divide(1,1,0.01,0.01);
	c10.cd(1);
	//c10.GetPad(1)->SetLogy();
	histoLowerEnergyPhoton.Draw();
	c10.Print("histoLowerEnergyPhoton.pdf");

	TCanvas c11;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
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
	c12.GetPad(1)->SetLogy();
	histohelicityAngle.Draw();
	c12.Print("histohelicityAngle.pdf");

	TCanvas c13;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c13.SetCanvasSize(1000,650);
	c13.Divide(1,1,0.01,0.01);
	c13.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoconeEnergy.Draw();
	c13.Print("ConeEnergy.pdf");

	TCanvas c14;
	gStyle->SetOptStat(111111);
	gStyle->SetPalette( 1 );
	c14.SetCanvasSize(1000,650);
	c14.Divide(1,1,0.01,0.01);
	c14.cd(1);
	//c13.GetPad(1)->SetLogy();
	energyOfPhotons.Draw();
	c14.Print("PhotonEnergy.pdf");

	TCanvas c15;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c15.SetCanvasSize(1000,650);
	c15.Divide(1,1,0.01,0.01);
	c15.cd(1);
	//c15.GetPad(1)->SetLogy();
	histoPhi.Draw();
	c15.Print("PhotonPhi.pdf");

	TCanvas c16;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c16.SetCanvasSize(1000,650);
	c16.Divide(1,1,0.01,0.01);
	c16.cd(1);
	//c15.GetPad(1)->SetLogy();
//	histogram.Draw("colz");
	//c16.Print("ISR.pdf");

	TCanvas c17;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c17.SetCanvasSize(1000,650);
	c17.Divide(1,1,0.01,0.01);
	c17.cd(1);
	//c15.GetPad(1)->SetLogy();
//	histoangleBetweenPhotons.Draw();
	//c17.Print("angleBetweenPhotons.pdf");

	TCanvas c18;
	gStyle->SetPalette( 1 );
	c18.SetCanvasSize(1000,650);
	c18.Divide(1,1,0.01,0.01);
	c18.cd(1);
	//c15.GetPad(1)->SetLogy();
	gStyle->SetOptStat(111111);
	histoNumberPhotonsbyEvent.Draw();
	c18.Print("NumberPhotonsbyEvnt.pdf");

	TCanvas c19;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c19.SetCanvasSize(1000,650);
	c19.Divide(1,1,0.01,0.01);
	c19.cd(1);
	//c19.GetPad(1)->SetLogy();
//	anglePhotonParticle.Draw();
	//c19.Print("PhotonParticleAngle.pdf");

	TCanvas c20;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c20.SetCanvasSize(1000,650);
	c20.Divide(1,1,0.01,0.01);
	c20.cd(1);
	//c20.GetPad(1)->SetLogy();
//	anglePhotonParticlewE.Draw();
	//c20.Print("PhotonParticleAngleWE.pdf");

	TCanvas c21;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c21.SetCanvasSize(1000,650);
	c21.Divide(1,1,0.01,0.01);
	c21.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoconeEnergyFilter.Draw();
	c21.Print("ConeEnergyFilter.pdf");

	TCanvas c22;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c22.SetCanvasSize(1000,650);
	c22.Divide(1,1,0.01,0.01);
	c22.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoHardPhotonsByEvent.Draw();
	c22.Print("HardPhotonsByEvent.pdf");

	TCanvas c23;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c23.SetCanvasSize(1000,650);
	c23.Divide(1,1,0.01,0.01);
	c23.cd(1);
//	c23.GetPad(1)->SetLogy();
	histoPtMax.Draw();
	c23.Print("histoMaxPt.pdf");



	TCanvas c24;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c24.SetCanvasSize(1000,650);
	c24.Divide(1,1,0.01,0.01);
	c24.cd(1);
//	c24.GetPad(1)->SetLogy();
	histoPt2Max.Draw();
	c24.Print("histo2ndMax.pdf");

	TCanvas c25;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c25.SetCanvasSize(1000,650);
	c25.Divide(1,1,0.01,0.01);
	c25.cd(1);
	//c13.GetPad(1)->SetLogy();
	//histoUnknownParticle.Draw();
	//c25.Print("XInvariantM.pdf");

	TCanvas c26;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c26.SetCanvasSize(1000,650);
	c26.Divide(1,1,0.01,0.01);
	c26.cd(1);
	//c13.GetPad(1)->SetLogy();
	histoPt.Draw();
	c26.Print("PhotonsPt.pdf");


	TCanvas c27;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c27.SetCanvasSize(1000,650);
	c27.Divide(1,1,0.01,0.01);
	c27.cd(1);
	//c13.GetPad(1)->SetLogy();
	ptOfPhotons.Draw();
	c27.Print("ptOfPhotons.pdf");

	TCanvas c28;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c28.SetCanvasSize(1000,650);
	c28.Divide(1,1,0.01,0.01);
	c28.cd(1);
	//c13.GetPad(1)->SetLogy();
	pLOfPhotons.Draw();
	c28.Print("pLOfPhotons.pdf");

	TCanvas c29;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c29.SetCanvasSize(1000,650);
	c29.Divide(1,1,0.01,0.01);
	c29.cd(1);
	//c13.GetPad(1)->SetLogy();
	ptOfHiggsPhotons.Draw();
	c29.Print("ptOfHiggsPhotons.pdf");

	TCanvas c30;
	gStyle->SetPalette( 1 );
	gStyle->SetOptStat(111111);
	c30.SetCanvasSize(1000,650);
	c30.Divide(1,1,0.01,0.01);
	c30.cd(1);
	//c13.GetPad(1)->SetLogy();
	pLOfHiggsPhotons.Draw();
	c30.Print("pLOfHiggsPhotons.pdf");

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

	return 0;
}



Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;
	UInt_t nFirstJob = 1;
	if(argc>iarg) nFirstJob = atoi(argv[iarg]); iarg++;
	UInt_t nLastJob = 10;
	if(argc>iarg) nLastJob   = atoi(argv[iarg]); iarg++;

	TString fName = "name";
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "pregencut.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
