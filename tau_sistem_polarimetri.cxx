/*
		Polarimetri za 350 GeV (Jinx)
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Napravljeno: 7530. šumopad 28.
		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		Autori: KG, Gordana, Goran
*/

#ifndef __CINT__
	#include "TROOT.h"
	#include "TFile.h"
	#include "Riostream.h"
	#include "TFile.h"
	#include "TH1.h"
	#include "TH2.h"
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
	#include "math.h"
	#include <algorithm>
	#include "TMatrixD.h"

// LCIO includes
	#include "lcio.h"
	#include <IOIMPL/LCFactory.h>
	#include <IMPL/LCCollectionVec.h>
	#include <EVENT/MCParticle.h>
	#include "EVENT/ReconstructedParticle.h"
	#include <IMPL/CalorimeterHitImpl.h>
	#include <IMPL/MCParticleImpl.h>
	#include <IMPL/ReconstructedParticleImpl.h>
	#include <UTIL/LCRelationNavigator.h>
	#include "EVENT/LCRelation.h"
	#include <UTIL/LCTOOLS.h>
	#include <UTIL/PIDHandler.h>
	#include <Exceptions.h>
	#include "TMatrixD.h"
#endif


#include "stdlib.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <array>

using namespace std;

// --------------------------------------------------------------------------------------------------------------------

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
	#ifdef __CINT__
		gSystem -> Load("${LCIO}/lib/liblcio.so");
		gSystem -> Load("${LCIO}/lib/liblcioDict.so");
	#endif



		Float_t P_Tauon_minus;
		Float_t P_Tauon_plus;
		Float_t P_Pion_minus;
		Float_t P_Pion_plus;
		Float_t tauMinus_azimut;
		Float_t tauPlus_azimut;
		Float_t pionMinus_azimut;
		Float_t pionMinus_polar;
		Float_t higsMinus_azimut;
		Float_t higsPlus_azimut;
		Float_t tauMinus_polar;
		Float_t tauPlus_polar;
		Float_t pionPlus_azimut;
		Float_t pionPlus_polar;
		Float_t higsMinus_polar;
		Float_t higsPlus_polar;
		Float_t pionMinus_Higgs;
		Float_t deltaPfi;

	// TFile rootFile("Polarimetar_Jinx_MC.root","RECREATE", "Polarimetar_Jinx_MC", 1);
	TTree JinxTree ("JinxTree", "Generator particle tree");

	TTree polarimeterTree ("polarimeterTree", "Generator particle tree");

	polarimeterTree.Branch("tauMinus_azimut", &tauMinus_azimut, "tauMinus_azimut");
	polarimeterTree.Branch("tauPlus_azimut", &tauPlus_azimut, "tauPlus_azimut");
		polarimeterTree.Branch("tauMinus_polar", &tauMinus_polar, "tauMinus_polar");
		polarimeterTree.Branch("tauPlus_polar", &tauPlus_polar, "tauPlus_polar");
	polarimeterTree.Branch("higsMinus_azimut", &higsMinus_azimut, "higsMinus_azimut");
	polarimeterTree.Branch("higsPlus_azimut", &higsPlus_azimut, "higsPlus_azimut");
	polarimeterTree.Branch("higsMinus_polar", &higsMinus_polar, "higsMinus_polar");
	polarimeterTree.Branch("higsPlus_polar", &higsPlus_polar, "higsPlus_polar");
	polarimeterTree.Branch("pionMinus_azimut", &pionMinus_azimut, "pionMinus_azimut");
	polarimeterTree.Branch("pionMinus_polar", &pionMinus_polar, "pionMinus_polar");
	polarimeterTree.Branch("pionPlus_azimut", &pionPlus_azimut, "pionPlus_azimut");
	polarimeterTree.Branch("pionPlus_polar", &pionPlus_polar, "pionPlus_polar");
	polarimeterTree.Branch("deltaPfi", &deltaPfi, "deltaPfi");

//	polarimeterTree.Branch("pionMinus_Higgs", &pionMinus_Higgs, "pionMinus_Higgs");

	// Int_t N_piona_plus = 0;
	// Int_t N_piona_minus = 0;
	// Float_t Polarimetar_6_minus;
	// Float_t Polarimetar_6_plus;
	Double_t Fi_Plus, Fi_Minus, Delta_Fi;

	TH1F histogram_Fi_Minus( "histogram_Fi_minus", " Ugao #it{#Phi}^{-}; #it{#Phi}^{-} ", 45, 0., 3.14);
	TH1F histogram_Fi_Plus( "histogram_Fi_Plus", " Ugao #it{#Phi}^{+}; #it{#Phi}^{+} ", 45, 0., 3.14);
	TH1F histogram_DeltaFi( "histogram_DeltaFi", " Ugao #Delta_{#it{#Phi}}; #Delta_{#it{#Phi}}", 45, -3.14, 3.14);
	TH1F histogram_TauonMinusko( "histogram_TauonMinusko", " Masa tauona minus; #it{m} (GeV) ", 40, 0, 4 );
	TH1F histogram_TauonPlusko( "histogram_TauonPlusko", " Masa tauona plus; #it{m} (GeV) ", 40, 0, 4 );
 	TH1F histogram_Higs_Tauoni( "histogram_Higs_Tauoni", " Masa Higsa od tauona; #it{m} (GeV) ", 10, 120, 130 );
 	TH1F histogram_Higs_Pioni( "histogram_Higs_Pioni", " Masa Higsa od produkta raspada; #it{m} (GeV) ", 10, 120, 130 );
	TH1F histogram_Ugao_tauon_foton_minus3("histogram_Ugao_tauon_foton_minus3", " Ugao izmedju tauona minus i fotona; #it{#theta} ", 100, 0, 3.14 );
	TH1F histogram_Ugao_tauon_foton_plus3("histogram_Ugao_tauon_foton_plus3", " Ugao izmedju tauona plus i fotona; #it{#theta} ", 100, 0, 3.14 );
	TH1F histogram_azimut("testPhi", " testPhi; #phi (rad) ", 90, 0, 3.14 );

	JinxTree.Branch("P_Tauon_minus", &P_Tauon_minus, "P_Tauon_minus");
	JinxTree.Branch("P_Tauon_plus", &P_Tauon_plus, "P_Tauon_plus");
	JinxTree.Branch("P_Pion_minus", &P_Pion_minus, "P_Pion_minus");
	JinxTree.Branch("P_Pion_plus", &P_Pion_plus, "P_Pion_plus");
	// JinxTree.Branch("Polarimetar_6_minus", &Polarimetar_6_minus, "Polarimetar_6_minus");
	// JinxTree.Branch("Polarimetar_6_plus", &Polarimetar_6_plus, "Polarimetar_6_plus");
	// JinxTree.Branch("Delta_Fi", &Delta_Fi, "Delta_Fi");

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader();
	TString fName = fn;
	stringstream fNameStream;

	// Int_t N_event_ukupno = 0;
	// Int_t N_isolep = 0;
	// Int_t N_signal = 0;
	// const Double_t mH_theory = 125.0;
	// bool right_evt = false;

	// N_event = 0;

	Int_t N_pion_minus_od_tauona = 0;
	Int_t N_pion_plus_od_tauona = 0;
	Int_t N_signjala = 0;
	Int_t N_tau_minus_from_photon = 0;
	Int_t N_tau_plus_from_photon = 0;
	Int_t N_evt_with_photon = 0;


	// Petlja koja iščitava .slcio fajlove
	for(UInt_t iJob = nFirstJob; iJob <= nLastJob; iJob++)
	{
    cout << "            " << endl;
    cout << "============================================================== " << endl;
		cout << "Otvara se " << Form("%s%i.slcio", fName.Data(), iJob);

		try
		{
			lcReader -> open(Form("%s%i.slcio", fName.Data(), iJob));
		}

		catch(lcio::IOException &ex)
		{
			cout << ". Ne mere se otvorit.\n";
			continue;
		}

		cout << ". Učitavanje.\n";

		Int_t N_event = 0;

		// Petlja po dogadjajima
		EVENT::LCEvent* evt = 0;

		while( (evt = lcReader -> readNextEvent()) != 0)
		{
			N_event++;
			bool B_tauon_potomak_od_tauona_minus = false;
			bool B_foton_potomak_od_tauona_minus = false;

			bool B_tauon_potomak_od_tauona_plus = false;
			bool B_foton_potomak_od_tauona_plus = false;


			std::vector<std::string> colNames = *evt -> getCollectionNames();

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticlesSkimmed");
			// IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFANewPFOs");
			// IMPL::LCCollectionVec* colJet = (IMPL::LCCollectionVec*) evt -> getCollection("twoRefJetsZep");
			// IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_2");
			// IMPL::LCCollectionVec* jets4 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_4");
			// IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected");

			TLorentzVector Higs_pocetni;
			TLorentzVector Tauon_minus;
			TLorentzVector Tauon_plus;
			vector <TLorentzVector> Pion_minus;
			vector <TLorentzVector> Pion_plus;

			/*Int_t N_tauona_plus = 0;
			Int_t N_tauona_minus = 0;

			Int_t N_tauona_plus_ukupno = 0;
			Int_t N_tauona_minus_ukupno = 0;*/

			TLorentzVector TLVpionMinus, TLVtauonskiNeutrino, TauonMinusko;
			TLorentzVector TLVpionPlus, TLVtauonskiAntineutrino, TauonPlusko;

			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				// cout << "usao u for po cesticama" << endl;

				if( mcParticle -> getPDG() == 25 && (mcParticle -> getDaughters()[0] -> getPDG() == 15 && mcParticle -> getDaughters()[1] -> getPDG() == -15) )
				{
			//		cout << "   " << endl;
			//		cout << "------------" << endl;
			//		cout << N_event << ". događaj" << endl;
				//	cout << "------------" << endl;

					EVENT::MCParticle* tauonMinus;
					EVENT::MCParticle* tauonPlus;

					EVENT::MCParticle* tauonMinusTemp;
					EVENT::MCParticle* tauonPlusTemp;

//				EVENT::MCParticle* pionMinus;
//				EVENT::MCParticle* pionPlus;

/*				EVENT::MCParticle* mionMinus;
					EVENT::MCParticle* mionPlus;
					EVENT::MCParticle* mionskiAntineutrino;
					EVENT::MCParticle* mionskiNeutrino; */

					const double *p = mcParticle -> getMomentum();
					double e = mcParticle -> getEnergy();

					Higs_pocetni.SetPxPyPzE( p[0], p[1], p[2], e );

					tauonMinus = mcParticle -> getDaughters()[0];
					tauonPlus = mcParticle -> getDaughters()[1];

					tauonMinusTemp = mcParticle -> getDaughters()[0];
					tauonPlusTemp = mcParticle -> getDaughters()[1];

					if (tauonMinus->getDaughters().size() > 0 && tauonMinus->getDaughters()[0]->getPDG()==15) tauonMinus = tauonMinus->getDaughters()[0];
					if (tauonPlus->getDaughters().size() > 0 && tauonPlus->getDaughters()[0]->getPDG()==-15) tauonPlus = tauonPlus->getDaughters()[0];

					bool B_pionMinus = false;
					bool B_pionPlus = false;


					if ( tauonMinusTemp -> getDaughters().size() == 1 )
					{
						tauonMinus = tauonMinusTemp -> getDaughters()[0] -> getDaughters()[0];
						tauonPlus = tauonMinusTemp -> getDaughters()[0] -> getDaughters()[1];
				//		cout <<"94->15->[0]->getPDG: "<< tauonMinus->getDaughters()[0]->getPDG()<<endl;
				//		cout <<"94->-15->[0]->getPDG: "<< tauonPlus->getDaughters()[0]->getPDG()<<endl;
						if (tauonMinus->getDaughters()[0]->getPDG()==15) tauonMinus = tauonMinus->getDaughters()[0];
						if (tauonPlus->getDaughters()[0]->getPDG()==-15) tauonPlus = tauonPlus->getDaughters()[0];


					}

					if ( tauonPlusTemp -> getDaughters().size() == 1 )
					{
						tauonMinus = tauonPlusTemp -> getDaughters()[0] -> getDaughters()[0];
						tauonPlus = tauonPlusTemp -> getDaughters()[0] -> getDaughters()[1];
						if (tauonMinus->getDaughters()[0]->getPDG()==15) tauonMinus = tauonMinus->getDaughters()[0];
						if (tauonPlus->getDaughters()[0]->getPDG()==-15) tauonPlus = tauonPlus->getDaughters()[0];
					}

// -------------------------------------------------------------------------------------

					const double *pTauMinus = tauonMinus -> getMomentum();
					double eTauMinus = tauonMinus -> getEnergy();

					Tauon_minus.SetPxPyPzE(pTauMinus[0],pTauMinus[1],pTauMinus[2],eTauMinus);

					TLorentzVector TLVtauon_potomak_od_tauona_minus, TLVfoton_potomak_od_tauona_minus;
					TLorentzVector TLVtauon_potomak_od_tauona_plus, TLVfoton_potomak_od_tauona_plus;

					for ( int i = 0; i < (int) tauonMinus -> getDaughters().size(); i++ )
					{
						Int_t PDG = tauonMinus -> getDaughters()[i] -> getPDG();

						const double *p = tauonMinus -> getDaughters()[i] -> getMomentum();
						double e = tauonMinus -> getDaughters()[i] -> getEnergy();

						TLorentzVector temp;
						temp.SetPxPyPzE(p[0],p[1],p[2],e);

						if ( PDG == 16 )
						{
							TLVtauonskiNeutrino = temp;
						}

						if ( PDG == -211 )
						{
							B_pionMinus = true;

							TLVpionMinus = temp;

							N_pion_minus_od_tauona++;
						}

						if ( PDG == 15 )
						{
							B_tauon_potomak_od_tauona_minus = true;
							TLVtauon_potomak_od_tauona_minus = temp;
					//		TLVpionMinus = temp;
							N_tau_minus_from_photon++;
						}

						if ( PDG == 22 )
						{
							B_foton_potomak_od_tauona_minus = true;
							TLVfoton_potomak_od_tauona_minus = temp;
						}

					}

					if ( B_pionMinus  )
					{
						TauonMinusko = TLVpionMinus + TLVtauonskiNeutrino;

						histogram_TauonMinusko.Fill( TauonMinusko.M() );
					}

					if ( B_tauon_potomak_od_tauona_minus && B_foton_potomak_od_tauona_minus )
					{
						TVector3 tauon_potomak_od_tauona_minus3(TLVtauon_potomak_od_tauona_minus.Px(),TLVtauon_potomak_od_tauona_minus.Py(),TLVtauon_potomak_od_tauona_minus.Pz()) ;
						TVector3 foton_potomak_od_tauona_minus3(TLVfoton_potomak_od_tauona_minus.Px(),TLVfoton_potomak_od_tauona_minus.Py(),TLVfoton_potomak_od_tauona_minus.Pz());

						histogram_Ugao_tauon_foton_minus3.Fill( tauon_potomak_od_tauona_minus3.Angle(foton_potomak_od_tauona_minus3) );
					}

// -------------------------------------------------------------------------------------

					const double *pTauPlus = tauonPlus -> getMomentum();
					double eTauPlus = tauonPlus -> getEnergy();

					Tauon_plus.SetPxPyPzE(pTauPlus[0],pTauPlus[1],pTauPlus[2],eTauPlus);

					for ( int i = 0; i < (int) tauonPlus -> getDaughters().size(); i++ )
					{
						Int_t PDG = tauonPlus -> getDaughters()[i] -> getPDG();

						const double *p = tauonPlus -> getDaughters()[i] -> getMomentum();
						double e = tauonPlus -> getDaughters()[i] -> getEnergy();

						TLorentzVector temp;
						temp.SetPxPyPzE(p[0],p[1],p[2],e);

						if ( PDG == -16 )
						{
							TLVtauonskiAntineutrino = temp;
						}

						if ( PDG == 211 )
						{
							B_pionPlus = true;

							TLVpionPlus = temp;

							N_pion_plus_od_tauona++;
						}

						if ( PDG == -15 )
						{
							B_tauon_potomak_od_tauona_plus = true;
							TLVtauon_potomak_od_tauona_plus = temp;
					//		TLVpionPlus = temp;
							N_tau_plus_from_photon++;
						}

						if ( PDG == 22 )
						{
							B_foton_potomak_od_tauona_plus = true;
							TLVfoton_potomak_od_tauona_plus = temp;
						}
					}

					if ( B_pionPlus )
					{
						TauonPlusko = TLVpionPlus + TLVtauonskiAntineutrino;

						histogram_TauonPlusko.Fill( TauonPlusko.M() );
					}

					if ( B_tauon_potomak_od_tauona_plus && B_foton_potomak_od_tauona_plus )
					{
						TVector3 tauon_potomak_od_tauona_plus3(TLVtauon_potomak_od_tauona_plus.Px(),TLVtauon_potomak_od_tauona_plus.Py(),TLVtauon_potomak_od_tauona_plus.Pz()) ;
						TVector3 foton_potomak_od_tauona_plus3(TLVfoton_potomak_od_tauona_plus.Px(),TLVfoton_potomak_od_tauona_plus.Py(),TLVfoton_potomak_od_tauona_plus.Pz());

						histogram_Ugao_tauon_foton_plus3.Fill( tauon_potomak_od_tauona_plus3.Angle(foton_potomak_od_tauona_plus3) );
					}

				// histogram_Higs_Tauoni.Fill( TauonMinusko.M() + TauonPlusko.M() );
				// histogram_Higs_Pioni.Fill( TLVpionMinus.M() + TLVtauonskiNeutrino.M() + TLVpionPlus.M() + TLVtauonskiAntineutrino.M() );

				} // Kraj IF za tauone 15 ili -15

			} // Kraj FOR petlje po broju elemenata
			if (B_tauon_potomak_od_tauona_minus && B_tauon_potomak_od_tauona_plus ) N_evt_with_photon++;

			if ( TLVpionPlus.P() > 0 && TLVpionMinus.P() > 0 )
			{
				N_signjala++;

				P_Pion_minus = TLVpionMinus.P();
				P_Pion_plus = TLVpionPlus.P();

				TLorentzVector Higsonja = TauonMinusko + TauonPlusko;

				histogram_Higs_Tauoni.Fill( (TauonMinusko + TauonPlusko).M() );
				histogram_Higs_Pioni.Fill( (TLVpionMinus + TLVtauonskiNeutrino + TLVpionPlus + TLVtauonskiAntineutrino).M() );

				TVector3 BoostToTauonMinus = -( Tauon_minus.BoostVector() );
				TVector3 BoostToTauonPlus = -( Tauon_plus.BoostVector() );

				TLorentzVector TLVHigs_tau_minus = Higs_pocetni;
				TLorentzVector TLVHigs_tau_plus = Higs_pocetni;

				TLVHigs_tau_minus.Boost(BoostToTauonMinus);
				TLVHigs_tau_plus.Boost(BoostToTauonPlus);

				TLVpionMinus.Boost(BoostToTauonMinus);
				TLVpionPlus.Boost(BoostToTauonPlus);
				Tauon_minus.Boost(BoostToTauonMinus);
				Tauon_plus.Boost(BoostToTauonPlus);
				if (Tauon_minus.Px() > 0.00001 ) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;
				if ( Tauon_minus.Py() > 0.00001 ) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;
				if (Tauon_minus.Pz() > 0.00001) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;


				TVector3 pionMinus3(TLVpionMinus.Px(),TLVpionMinus.Py(),TLVpionMinus.Pz()) ;
				TVector3 pionPlus3(TLVpionPlus.Px(),TLVpionPlus.Py(),TLVpionPlus.Pz());
				TVector3 tauMinus3(Tauon_minus.Px(),Tauon_minus.Py(),Tauon_minus.Pz()) ;
				TVector3 tauPlus3(Tauon_plus.Px(),Tauon_plus.Py(),Tauon_plus.Pz()) ;


				TLVtauonskiNeutrino.Boost(BoostToTauonMinus);
				TLVtauonskiAntineutrino.Boost(BoostToTauonPlus);

				TVector3 Higs_tauon_minus3( TLVHigs_tau_minus.Px(),TLVHigs_tau_minus.Py(),TLVHigs_tau_minus.Pz() );
				TVector3 Higs_tauon_plus3( TLVHigs_tau_plus.Px(),TLVHigs_tau_plus.Py(),TLVHigs_tau_plus.Pz() );

				Fi_Minus = Higs_tauon_minus3.Angle( pionMinus3 );
				Fi_Plus = Higs_tauon_plus3.Angle( pionPlus3 );

				histogram_Fi_Minus.Fill(Fi_Minus);
				histogram_Fi_Plus.Fill(Fi_Plus);

				Delta_Fi = Fi_Plus - Fi_Minus;

	//			cout <<"Higgs px: "<< Higs_tauon_minus3.X()<<", py: "<< Higs_tauon_minus3.Y()<<", pz: "<<Higs_tauon_minus3.Z()<<endl;
				//rotacija u prostoru
				//   TRotation r ();        // r initialized as identity
				 //  XYZVector u (0,1,0);


				TVector3 x(1,0,0);
				TVector3 y(0,1,0);
				TVector3 z(0,0,1);


				//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-****-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
				//KG Vršimo transformaciju koordinatnih sistema na takav način da Higgsov pravac postaje x-osa






				//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-****-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


	//			Higs_tauon_minus3.Print();

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-****-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
				//KG Vršimo transformaciju koordinatnih sistema na takav način da Higgsov pravac postaje x-osa
				//KG y osa je z x Higgs, dok novu z osu funkcija sama definiše


		/*		TVector3 directionMinusTest = Higs_tauon_minus3;
				directionMinusTest *= 1./directionMinusTest.Mag(); // KG minus zato što je pravac x ose suprotan pravcu kretanja Higsa
				TVector3 rotatedYAxisMinusTest = z.Cross(directionMinusTest);
				rotatedYAxisMinusTest *= 1.0/rotatedYAxisMinusTest.Mag();
				TRotation rMinusTest;
				rMinusTest.SetXAxis(directionMinusTest, rotatedYAxisMinusTest);
				rMinusTest.Invert();

				TVector3 testVec1 = Higs_tauon_minus3;

				testVec1.Transform(rMinusTest);

				testVec1.Print();*/


		/*		TVector3 directionPlus = Higs_tauon_plus3;
				directionPlus *= 1./directionPlus.Mag();

				TVector3 rotatedYAxisPlus = z.Cross(directionPlus);
				rotatedYAxisPlus *= 1.0/rotatedYAxisPlus.Mag();

				TRotation rPlus;
				rPlus.SetXAxis(directionPlus, rotatedYAxisPlus);
				rPlus.Invert();

				Higs_tauon_plus3.Transform(rPlus);
				pionPlus3.Transform(rPlus);
				tauPlus3.Transform(rPlus);


				TVector3 directionMinus = Higs_tauon_minus3;
				directionMinus *= -1./directionMinus.Mag(); // KG minus zato što je pravac x ose suprotan pravcu kretanja Higsa


				TVector3 rotatedYAxisMinus = z.Cross(directionMinus);
				rotatedYAxisMinus *= 1.0/rotatedYAxisMinus.Mag();

				TRotation rMinus;
				rMinus.SetXAxis(directionMinus, rotatedYAxisMinus);
				rMinus.Invert();

				Higs_tauon_minus3.Transform(rMinus);
				pionMinus3.Transform(rMinus);
				tauMinus3.Transform(rMinus);*/

//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-****-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

			//	pionMinus3.Print();
			//	pionPlus3.Print();
			//	Higs_tauon_minus3.Print();
			//	Higs_tauon_plus3.Print();

				cout <<"Higgs u tau minus:" <<endl;
				Higs_tauon_minus3.Print();
				cout <<"Higgs u tau plus:" <<endl;
				Higs_tauon_plus3.Print();


				Double_t pMinusPhi = pionMinus3.Phi();
				if (pMinusPhi < 0) pMinusPhi = 2*M_PI-abs(pMinusPhi);


				Double_t pPlusPhi = pionPlus3.Phi();
//
				if (pPlusPhi < 0) pPlusPhi = 2*M_PI-abs(pPlusPhi);
			//	if (Higs_tauon_minus3.Phi() < 0)  Higs_tauon_minus3.Print();


				pionMinus_azimut = pMinusPhi;
				pionMinus_polar = pionMinus3.Theta();
				pionPlus_azimut = pPlusPhi;
				pionPlus_polar = pionPlus3.Theta();
				higsMinus_azimut = Higs_tauon_minus3.Phi();
				higsMinus_polar = Higs_tauon_minus3.Theta();
				higsPlus_azimut = Higs_tauon_plus3.Phi();
				higsPlus_polar = Higs_tauon_plus3.Theta();
				tauMinus_polar = tauMinus3.Theta();
				tauPlus_polar = tauPlus3.Theta();
				tauMinus_azimut = tauMinus3.Phi();
				tauPlus_azimut = tauPlus3.Phi();
				deltaPfi = pPlusPhi - pMinusPhi;



				polarimeterTree.Fill();
			}

		} // Kraj WHILE petlje

		JinxTree.Fill();

		lcReader -> close();

	} // Kraj FOR petlje koja iščitava .slcio fajlove

	cout << "Ukupan broj signjala " << N_signjala << endl;

	cout << "Ukupan broj piona minus od tauona " << N_pion_minus_od_tauona << endl;
	cout << "Ukupan broj piona plus od tauona " << N_pion_plus_od_tauona << endl;
	//N_tau_minus_from_photon
//	cout << "Ukupan broj tau minus sa fotonom: " << N_tau_minus_from_photon << endl;
//	cout << "Ukupan broj tau plus sa fotonom: " << N_tau_plus_from_photon << endl;
//	cout << "Ukupan broj dogadjaja sa fotonom: " << N_evt_with_photon << endl;


//	JinxTree.Fill();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	histogram_Fi_Minus.Write();
	histogram_Fi_Plus.Write();

	histogram_DeltaFi.Write();
	histogram_TauonMinusko.Write();
	histogram_TauonPlusko.Write();
	histogram_Higs_Tauoni.Write();
	histogram_Higs_Pioni.Write();
	histogram_Ugao_tauon_foton_minus3.Write();
	histogram_Ugao_tauon_foton_plus3.Write();

	JinxTree.Write();
	polarimeterTree.Write();

	// evtTree.Write();
	rootFile.Write();
	rootFile.Close();

	return 0;

} // Kraj funkcije

// ------------------------------------------------------------------

Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;

	UInt_t nFirstJob = 1;

	if(argc > iarg) {nFirstJob = atoi(argv[iarg]); iarg++;}

	UInt_t nLastJob = 50;

	if(argc > iarg) {nLastJob  = atoi(argv[iarg]); iarg++;}

	TString fName = "test_";
	if(argc > iarg) fName = argv[iarg]; //iarg++;

	TString rfName = "Polarimetar_Jinx_MC.root";
	if(argc > iarg) rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
