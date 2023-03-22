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

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn,
		const char * rfn) {
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
	Float_t pionPlus_azimut, higs_azimut;
	Float_t pionPlus_polar;
	Float_t higsMinus_polar, higs_polar;
	Float_t higsPlus_polar;
	Float_t pionMinus_Higgs;
	Float_t deltaPfi, deltaPfi1;

	Float_t Higgs3_x, Higgs3_y, Higgs3_z, HSHigs_pt, HigsSys_Mass;
	Float_t tauMinus3_x, tauMinus3_y, tauMinus3_z, HS_tauMinus_pt,
			HS_tauMinus_Mass;
	Float_t tauPlus3_x, tauPlus3_y, tauPlus3_z, HS_tauPlus_pt, HS_tauPlus_Mass;
	Float_t pionPlus3_x, pionPlus3_y, pionPlus3_z, HS_pionPlus_pt,
			HS_pionPlus_Mass;
	Float_t pionMinus3_x, pionMinus3_y, pionMinus3_z, HS_pionMinus_pt,
			HS_pionMinus_Mass;

	// TFile rootFile("Polarimetar_Jinx_MC.root","RECREATE", "Polarimetar_Jinx_MC", 1);
	TTree JinxTree("JinxTree", "Generator particle tree");

	TTree polarimeterTree("polarimeterTree", "Generator particle tree");

	polarimeterTree.Branch("tauMinus_azimut", &tauMinus_azimut,
			"tauMinus_azimut");
	polarimeterTree.Branch("tauMinus_polar", &tauMinus_polar, "tauMinus_polar");

	polarimeterTree.Branch("tauPlus_azimut", &tauPlus_azimut, "tauPlus_azimut");
	polarimeterTree.Branch("tauPlus_polar", &tauPlus_polar, "tauPlus_polar");

	polarimeterTree.Branch("higs_azimut", &higs_azimut, "higs_azimut");
	polarimeterTree.Branch("higs_polar", &higs_polar, "higs_polar");

	polarimeterTree.Branch("pionMinus_azimut", &pionMinus_azimut,
			"pionMinus_azimut");
	polarimeterTree.Branch("pionMinus_polar", &pionMinus_polar,
			"pionMinus_polar");

	polarimeterTree.Branch("pionPlus_azimut", &pionPlus_azimut,
			"pionPlus_azimut");
	polarimeterTree.Branch("pionPlus_polar", &pionPlus_polar, "pionPlus_polar");

	polarimeterTree.Branch("deltaPfi", &deltaPfi, "deltaPfi");
	polarimeterTree.Branch("deltaPfi1", &deltaPfi1, "deltaPfi1");

	TTree distroTree("distroTree", "Generator particle tree");
	distroTree.Branch("HS_tauMinus_Mass", &HS_tauMinus_Mass,
			"HS_tauMinus_Mass");
	distroTree.Branch("HSHigs_pt", &HSHigs_pt, "HSHigs_pt");
	distroTree.Branch("Higgs3_x", &Higgs3_x, "Higgs3_x");
	distroTree.Branch("Higgs3_y", &Higgs3_y, "Higgs3_y");
	distroTree.Branch("Higgs3_z", &Higgs3_z, "Higgs3_z");

	distroTree.Branch("HS_tauMinus_Mass", &HS_tauMinus_Mass,
			"HS_tauMinus_Mass");
	distroTree.Branch("HS_tauMinus_pt", &HS_tauMinus_pt, "HS_tauMinus_pt");
	distroTree.Branch("tauMinus3_x", &tauMinus3_x, "tauMinus3_x");
	distroTree.Branch("tauMinus3_y", &tauMinus3_y, "tauMinus3_y");
	distroTree.Branch("tauMinus3_z", &tauMinus3_z, "tauMinus3_z");

	distroTree.Branch("HS_tauPlus_Mass", &HS_tauPlus_Mass, "HS_tauPlus_Mass");
	distroTree.Branch("HS_tauPlus_pt", &HS_tauPlus_pt, "HS_tauPlus_pt");
	distroTree.Branch("tauPlus3_x", &tauPlus3_x, "tauPlus3_x");
	distroTree.Branch("tauPlus3_y", &tauPlus3_y, "tauPlus3_y");
	distroTree.Branch("tauPlus3_z", &tauPlus3_z, "tauPlus3_z");

	distroTree.Branch("HS_pionMinus_Mass", &HS_pionMinus_Mass,
			"HS_pionMinus_Mass");
	distroTree.Branch("HS_pionMinus_pt", &HS_pionMinus_pt, "HS_pionMinus_pt");
	distroTree.Branch("pionMinus3_x", &pionMinus3_x, "pionMinus3_x");
	distroTree.Branch("pionMinus3_y", &pionMinus3_y, "pionMinus3_y");
	distroTree.Branch("pionMinus3_z", &pionMinus3_z, "pionMinus3_z");

	distroTree.Branch("HS_pionMinus_Mass", &HS_pionMinus_Mass,
			"HS_pionMinus_Mass");
	distroTree.Branch("HS_pionMinus_pt", &HS_pionMinus_pt, "HS_pionMinus_pt");
	distroTree.Branch("pionPlus3_x", &pionPlus3_x, "pionPlus3_x");
	distroTree.Branch("pionPlus3_y", &pionPlus3_y, "pionPlus3_y");
	distroTree.Branch("pionPlus3_z", &pionPlus3_z, "pionPlus3_z");

	float labHigs_Mass, labHigs_Polar, labHigs_Azi, labHigs_pt, labHigs_X,
			labHigs_Y, labHigs_Z, labHigs_E;
	float labTauM_Mass, labTauM_Polar, labTauM_Azi, labTauM_pt, labTauM_X,
			labTauM_Y, labTauM_Z, labTauM_E;
	float labTauP_Mass, labTauP_Polar, labTauP_Azi, labTauP_pt, labTauP_X,
			labTauP_Y, labTauP_Z, labTauP_E;
	float labpionM_Mass, labpionM_Polar, labpionM_Azi, labpionM_pt, labpionM_X,
			labpionM_Y, labpionM_Z, labpionM_E;
	float labpionP_Mass, labpionP_Polar, labpionP_Azi, labpionP_pt, labpionP_X,
			labpionP_Y, labpionP_Z, labpionP_E;

	TTree labTree("labTree", "Generator particle tree");
	labTree.Branch("labHigs_Mass", &labHigs_Mass, "labHigs_Mass");
	labTree.Branch("labHigs_Polar", &labHigs_Polar, "labHigs_Polar");
	labTree.Branch("labHigs_Azi", &labHigs_Azi, "labHigs_Azi");
	labTree.Branch("labHigs_E", &labHigs_E, "labHigs_E");
	labTree.Branch("labHigs_pt", &labHigs_pt, "labHigs_pt");
	labTree.Branch("labHigs_X", &labHigs_X, "labHigs_X");
	labTree.Branch("labHigs_Y", &labHigs_Y, "labHigs_Y");
	labTree.Branch("labHigs_Z", &labHigs_Z, "labHigs_Z");
	// tau minus
	labTree.Branch("labTauM_Mass", &labTauM_Mass, "labTauM_Mass");
	labTree.Branch("labTauM_Polar", &labTauM_Polar, "labTauM_Polar");
	labTree.Branch("labTauM_Azi", &labTauM_Azi, "labTauM_Azi");
	labTree.Branch("labTauM_E", &labTauM_E, "labTauM_E");
	labTree.Branch("labTauM_pt", &labTauM_pt, "labTauM_pt");
	labTree.Branch("labTauM_X", &labTauM_X, "labTauM_X");
	labTree.Branch("labTauM_Y", &labTauM_Y, "labTauM_Y");
	labTree.Branch("labTauM_Z", &labTauM_Z, "labTauM_Z");

//tau plus
	labTree.Branch("labTauP_Mass", &labTauP_Mass, "labTauP_Mass");
	labTree.Branch("labTauP_Polar", &labTauP_Polar, "labTauP_Polar");
	labTree.Branch("labTauP_Azi", &labTauP_Azi, "labTauP_Azi");
	labTree.Branch("labTauP_E", &labTauP_E, "labTauP_E");
	labTree.Branch("labTauP_pt", &labTauP_pt, "labTauP_pt");
	labTree.Branch("labTauP_X", &labTauP_X, "labTauP_X");
	labTree.Branch("labTauP_Y", &labTauP_Y, "labTauP_Y");
	labTree.Branch("labTauP_Z", &labTauP_Z, "labTauP_Z");

//pion minus
	labTree.Branch("labpionM_Mass", &labpionM_Mass, "labpionM_Mass");
	labTree.Branch("labpionM_Polar", &labpionM_Polar, "labpionM_Polar");
	labTree.Branch("labpionM_Azi", &labpionM_Azi, "labpionM_Azi");
	labTree.Branch("labpionM_pt", &labpionM_pt, "labpionM_pt");
	labTree.Branch("labpionM_E", &labpionM_E, "labpionM_E");
	labTree.Branch("labpionM_X", &labpionM_X, "labpionM_X");
	labTree.Branch("labpionM_Y", &labpionM_Y, "labpionM_Y");
	labTree.Branch("labpionM_Z", &labpionM_Z, "labpionM_Z");
	//pion plus
	labTree.Branch("labpionp_Mass", &labpionP_Mass, "labpionP_Mass");
	labTree.Branch("labpionP_Polar", &labpionP_Polar, "labpionP_Polar");
	labTree.Branch("labpionP_Azi", &labpionP_Azi, "labpionP_Azi");
	labTree.Branch("labpionP_E", &labpionP_E, "labpionP_E");
	labTree.Branch("labpionP_pt", &labpionP_pt, "labpionP_pt");
	labTree.Branch("labpionP_X", &labpionP_X, "labpionP_X");
	labTree.Branch("labpionP_Y", &labpionP_Y, "labpionP_Y");
	labTree.Branch("labpionP_Z", &labpionP_Z, "labpionP_Z");

	// Higsov sistem reference

	float HS_Higs_Mass, HS_Higs_Polar, HS_Higs_Azi, HS_Higs_pt, HS_Higs_X,
			HS_Higs_Y, HS_Higs_Z, HS_Higs_E;
	float HS_TauM_Mass, HS_TauM_Polar, HS_TauM_Azi, HS_TauM_pt, HS_TauM_X,
			HS_TauM_Y, HS_TauM_Z, HS_TauM_E;
	float HS_TauP_Mass, HS_TauP_Polar, HS_TauP_Azi, HS_TauP_pt, HS_TauP_X,
			HS_TauP_Y, HS_TauP_Z, HS_TauP_E;
	float HS_pionM_Mass, HS_pionM_Polar, HS_pionM_Azi, HS_pionM_pt, HS_pionM_X,
			HS_pionM_Y, HS_pionM_Z, HS_pionM_E;
	float HS_pionP_Mass, HS_pionP_Polar, HS_pionP_Azi, HS_pionP_pt, HS_pionP_X,
			HS_pionP_Y, HS_pionP_Z, HS_pionP_E;
	float angleTaus, angleTausOriginal;

	TTree HSTree("HSTree", "Generator particle tree");
	HSTree.Branch("HS_Higs_Mass", &HS_Higs_Mass, "HS_Higs_Mass");
	HSTree.Branch("HS_Higs_Polar", &HS_Higs_Polar, "HS_Higs_Polar");
	HSTree.Branch("HS_Higs_Azi", &HS_Higs_Azi, "HS_Higs_Azi");
	HSTree.Branch("HS_Higs_E", &HS_Higs_E, "HS_Higs_E");
	HSTree.Branch("HS_Higs_pt", &HS_Higs_pt, "HS_Higs_pt");
	HSTree.Branch("HS_Higs_X", &HS_Higs_X, "HS_Higs_X");
	HSTree.Branch("HS_Higs_Y", &HS_Higs_Y, "HS_Higs_Y");
	HSTree.Branch("HS_Higs_Z", &HS_Higs_Z, "HS_Higs_Z");
	// tau minus
	HSTree.Branch("HS_TauM_Mass", &HS_TauM_Mass, "HS_TauM_Mass");
	HSTree.Branch("HS_TauM_Polar", &HS_TauM_Polar, "HS_TauM_Polar");
	HSTree.Branch("HS_TauM_Azi", &HS_TauM_Azi, "HS_TauM_Azi");
	HSTree.Branch("HS_TauM_E", &HS_TauM_E, "HS_TauM_E");
	HSTree.Branch("HS_TauM_pt", &HS_TauM_pt, "HS_TauM_pt");
	HSTree.Branch("HS_TauM_X", &HS_TauM_X, "HS_TauM_X");
	HSTree.Branch("HS_TauM_Y", &HS_TauM_Y, "HS_TauM_Y");
	HSTree.Branch("HS_TauM_Z", &HS_TauM_Z, "HS_TauM_Z");

//tau plus
	HSTree.Branch("HS_TauP_Mass", &HS_TauP_Mass, "HS_TauP_Mass");
	HSTree.Branch("HS_TauP_Polar", &HS_TauP_Polar, "HS_TauP_Polar");
	HSTree.Branch("HS_TauP_Azi", &HS_TauP_Azi, "HS_TauP_Azi");
	HSTree.Branch("HS_TauP_E", &HS_TauP_E, "HS_TauP_E");
	HSTree.Branch("HS_TauP_pt", &HS_TauP_pt, "HS_TauP_pt");
	HSTree.Branch("HS_TauP_X", &HS_TauP_X, "HS_TauP_X");
	HSTree.Branch("HS_TauP_Y", &HS_TauP_Y, "HS_TauP_Y");
	HSTree.Branch("HS_TauP_Z", &HS_TauP_Z, "HS_TauP_Z");

//pion minus
	HSTree.Branch("HS_pionM_Mass", &HS_pionM_Mass, "HS_pionM_Mass");
	HSTree.Branch("HS_pionM_Polar", &HS_pionM_Polar, "HS_pionM_Polar");
	HSTree.Branch("HS_pionM_Azi", &HS_pionM_Azi, "HS_pionM_Azi");
	HSTree.Branch("HS_pionM_E", &HS_pionM_E, "HS_pionM_E");
	HSTree.Branch("HS_pionM_pt", &HS_pionM_pt, "HS_pionM_pt");
	HSTree.Branch("HS_pionM_X", &HS_pionM_X, "HS_pionM_X");
	HSTree.Branch("HS_pionM_Y", &HS_pionM_Y, "HS_pionM_Y");
	HSTree.Branch("HS_pionM_Z", &HS_pionM_Z, "HS_pionM_Z");
	//pion plus
	HSTree.Branch("HS_pionP_Mass", &HS_pionP_Mass, "HS_pionP_Mass");
	HSTree.Branch("HS_pionP_Polar", &HS_pionP_Polar, "HS_pionP_Polar");
	HSTree.Branch("HS_pionP_Azi", &HS_pionP_Azi, "HS_pionP_Azi");
	HSTree.Branch("HS_pionP_E", &HS_pionP_E, "HS_pionP_E");
	HSTree.Branch("HS_pionP_pt", &HS_pionP_pt, "HS_pionP_pt");
	HSTree.Branch("HS_pionP_X", &HS_pionP_X, "HS_pionP_X");
	HSTree.Branch("HS_pionP_Y", &HS_pionP_Y, "HS_pionP_Y");
	HSTree.Branch("HS_pionP_Z", &HS_pionP_Z, "HS_pionP_Z");
	HSTree.Branch("angleTaus", &angleTaus, "angleTaus");
	HSTree.Branch("angleTausOriginal", &angleTausOriginal, "angleTausOriginal");

	//boostovan sistem gde je tau minus (0, 0, M) vektor

	float boosted_Higs_Mass, boosted_Higs_Polar, boosted_Higs_Azi,
			boosted_Higs_pt, boosted_Higs_X, boosted_Higs_Y, boosted_Higs_Z;
	float boosted_TauM_Mass, boosted_TauM_Polar, boosted_TauM_Azi,
			boosted_TauM_pt, boosted_TauM_X, boosted_TauM_Y, boosted_TauM_Z;
	float boosted_TauP_Mass, boosted_TauP_Polar, boosted_TauP_Azi,
			boosted_TauP_pt, boosted_TauP_X, boosted_TauP_Y, boosted_TauP_Z;
	float boosted_pionM_Mass, boosted_pionM_Polar, boosted_pionM_Azi,
			boosted_pionM_pt, boosted_pionM_X, boosted_pionM_Y, boosted_pionM_Z;
	float boosted_pionP_Mass, boosted_pionP_Polar, boosted_pionP_Azi,
			boosted_pionP_pt, boosted_pionP_X, boosted_pionP_Y, boosted_pionP_Z;

	TTree boostedTree("boostedTree", "Generator particle tree");
	boostedTree.Branch("boosted_Higs_Mass", &boosted_Higs_Mass,
			"boosted_Higs_Mass");
	boostedTree.Branch("boosted_Higs_Polar", &boosted_Higs_Polar,
			"boosted_Higs_Polar");
	boostedTree.Branch("boosted_Higs_Azi", &boosted_Higs_Azi,
			"boosted_Higs_Azi");
	boostedTree.Branch("boosted_Higs_pt", &boosted_Higs_pt, "boosted_Higs_pt");
	boostedTree.Branch("boosted_Higs_X", &boosted_Higs_X, "boosted_Higs_X");
	boostedTree.Branch("boosted_Higs_Y", &boosted_Higs_Y, "boosted_Higs_Y");
	boostedTree.Branch("boosted_Higs_Z", &boosted_Higs_Z, "boosted_Higs_Z");
	// tau minus
	boostedTree.Branch("boosted_TauM_Mass", &boosted_TauM_Mass,
			"boosted_TauM_Mass");
	boostedTree.Branch("boosted_TauM_Polar", &boosted_TauM_Polar,
			"boosted_TauM_Polar");
	boostedTree.Branch("boosted_TauM_Azi", &boosted_TauM_Azi,
			"boosted_TauM_Azi");
	boostedTree.Branch("boosted_TauM_pt", &boosted_TauM_pt, "boosted_TauM_pt");
	boostedTree.Branch("boosted_TauM_X", &boosted_TauM_X, "boosted_TauM_X");
	boostedTree.Branch("boosted_TauM_Y", &boosted_TauM_Y, "boosted_TauM_Y");
	boostedTree.Branch("boosted_TauM_Z", &boosted_TauM_Z, "boosted_TauM_Z");

//tau plus
	boostedTree.Branch("boosted_TauP_Mass", &boosted_TauP_Mass,
			"boosted_TauP_Mass");
	boostedTree.Branch("boosted_TauP_Polar", &boosted_TauP_Polar,
			"boosted_TauP_Polar");
	boostedTree.Branch("boosted_TauP_Azi", &boosted_TauP_Azi,
			"boosted_TauP_Azi");
	boostedTree.Branch("boosted_TauP_pt", &boosted_TauP_pt, "boosted_TauP_pt");
	boostedTree.Branch("boosted_TauP_X", &boosted_TauP_X, "boosted_TauP_X");
	boostedTree.Branch("boosted_TauP_Y", &boosted_TauP_Y, "boosted_TauP_Y");
	boostedTree.Branch("boosted_TauP_Z", &boosted_TauP_Z, "boosted_TauP_Z");

//pion minus
	boostedTree.Branch("boosted_pionM_Mass", &boosted_pionM_Mass,
			"boosted_pionM_Mass");
	boostedTree.Branch("boosted_pionM_Polar", &boosted_pionM_Polar,
			"boosted_pionM_Polar");
	boostedTree.Branch("boosted_pionM_Azi", &boosted_pionM_Azi,
			"boosted_pionM_Azi");
	boostedTree.Branch("boosted_pionM_pt", &boosted_pionM_pt,
			"boosted_pionM_pt");
	boostedTree.Branch("boosted_pionM_X", &boosted_pionM_X, "boosted_pionM_X");
	boostedTree.Branch("boosted_pionM_Y", &boosted_pionM_Y, "boosted_pionM_Y");
	boostedTree.Branch("boosted_pionM_Z", &boosted_pionM_Z, "boosted_pionM_Z");
	//pion plus
	boostedTree.Branch("boosted_pionP_Mass", &boosted_pionP_Mass,
			"boosted_pionP_Mass");
	boostedTree.Branch("boosted_pionP_Polar", &boosted_pionP_Polar,
			"boosted_pionP_Polar");
	boostedTree.Branch("boosted_pionP_Azi", &boosted_pionP_Azi,
			"boosted_pionP_Azi");
	boostedTree.Branch("boosted_pionP_pt", &boosted_pionP_pt,
			"boosted_pionP_pt");
	boostedTree.Branch("boosted_pionP_X", &boosted_pionP_X, "boosted_pionP_X");
	boostedTree.Branch("boosted_pionP_Y", &boosted_pionP_Y, "boosted_pionP_Y");
	boostedTree.Branch("boosted_pionP_Z", &boosted_pionP_Z, "boosted_pionP_Z");

//	polarimeterTree.Branch("pionMinus_Higgs", &pionMinus_Higgs, "pionMinus_Higgs");

// Int_t N_piona_plus = 0;
// Int_t N_piona_minus = 0;
// Float_t Polarimetar_6_minus;
// Float_t Polarimetar_6_plus;

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
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
	for (UInt_t iJob = nFirstJob; iJob <= nLastJob; iJob++) {
		cout << "            " << endl;
		cout
				<< "============================================================== "
				<< endl;
		cout << "Otvara se " << Form("%s%i.slcio", fName.Data(), iJob);

		try {
			lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));
		}

		catch (lcio::IOException &ex) {
			cout << ". Ne mere se otvorit.\n";
			continue;
		}

		cout << ". Učitavanje.\n";

		Int_t N_event = 0;

		// Petlja po dogadjajima
		EVENT::LCEvent* evt = 0;

		while ((evt = lcReader->readNextEvent()) != 0) {
			N_event++;
			bool B_tauon_potomak_od_tauona_minus = false;
			bool B_foton_potomak_od_tauona_minus = false;

			bool B_tauon_potomak_od_tauona_plus = false;
			bool B_foton_potomak_od_tauona_plus = false;

			std::vector<std::string> colNames = *evt->getCollectionNames();

			IMPL::LCCollectionVec* mcParticles =
					(IMPL::LCCollectionVec*) evt->getCollection(
							"MCParticlesSkimmed");
			// IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFANewPFOs");
			// IMPL::LCCollectionVec* colJet = (IMPL::LCCollectionVec*) evt -> getCollection("twoRefJetsZep");
			// IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_2");
			// IMPL::LCCollectionVec* jets4 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_4");
			// IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected");

			TLorentzVector Higs_pocetni;
			TLorentzVector Tauon_minus;
			TLorentzVector Tauon_plus;
			vector<TLorentzVector> Pion_minus;
			vector<TLorentzVector> Pion_plus;

			/*Int_t N_tauona_plus = 0;
			 Int_t N_tauona_minus = 0;

			 Int_t N_tauona_plus_ukupno = 0;
			 Int_t N_tauona_minus_ukupno = 0;*/

			TLorentzVector TLVpionMinus, TLVtauonskiNeutrino, TauonMinusko;
			TLorentzVector TLVpionPlus, TLVtauonskiAntineutrino, TauonPlusko;

			for (Int_t i = 0; i < mcParticles->getNumberOfElements(); i++) {
				IMPL::MCParticleImpl* mcParticle =
						(IMPL::MCParticleImpl*) mcParticles->getElementAt(i);
				// cout << "usao u for po cesticama" << endl;

				if (mcParticle->getPDG() == 25
						&& (mcParticle->getDaughters()[0]->getPDG() == 15
								&& mcParticle->getDaughters()[1]->getPDG()
										== -15)) {
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

					tauonMinus = mcParticle->getDaughters()[0];
					tauonPlus = mcParticle->getDaughters()[1];

					tauonMinusTemp = mcParticle->getDaughters()[0];
					tauonPlusTemp = mcParticle->getDaughters()[1];

					const double *p = mcParticle->getMomentum();
					double e = mcParticle->getEnergy();

					Higs_pocetni.SetPxPyPzE(p[0], p[1], p[2], e);

					// ************************************************************
					//provera da li su pocetni tauoni back to back u higsovom sistemu reference

					/*		const double *pM = tauonMinus->getMomentum();
					 double eM = tauonMinus -> getEnergy();

					 const double *pP = tauonPlus->getMomentum();
					 double eP = tauonPlus -> getEnergy();

					 TVector3 BoostToHiggs1 = -( Higs_pocetni.BoostVector() );
					 Higs_pocetni.Boost(BoostToHiggs1);
					 cout <<"higgs boost: "<<endl;
					 Higs_pocetni.Print();

					 TLorentzVector tauM ;
					 TLorentzVector tauP ;
					 tauM.SetPxPyPzE(pM[0],pM[1],pM[2], eM);
					 tauP.SetPxPyPzE(pP[0],pP[1],pP[2], eP);

					 tauM.Boost(BoostToHiggs1);

					 tauP.Boost(BoostToHiggs1);

					 TVector3 tauM3 (tauM.Px(), tauM.Py(), tauM.Pz());
					 TVector3 tauP3 (tauP.Px(), tauP.Py(), tauP.Pz());
					 cout <<"ugao izmedju tau čestica nakon boost a u Higsov sistem: "<< tauM3.Angle(tauP3)<<endl;*/

					// ************************************************************

					if (tauonMinus->getDaughters().size() > 0
							&& tauonMinus->getDaughters()[0]->getPDG() == 15)
						tauonMinus = tauonMinus->getDaughters()[0];
					if (tauonPlus->getDaughters().size() > 0
							&& tauonPlus->getDaughters()[0]->getPDG() == -15)
						tauonPlus = tauonPlus->getDaughters()[0];

					bool B_pionMinus = false;
					bool B_pionPlus = false;

					if (tauonMinusTemp->getDaughters().size() == 1) {
						tauonMinus =
								tauonMinusTemp->getDaughters()[0]->getDaughters()[0];
						tauonPlus =
								tauonMinusTemp->getDaughters()[0]->getDaughters()[1];
						//		cout <<"94->15->[0]->getPDG: "<< tauonMinus->getDaughters()[0]->getPDG()<<endl;
						//		cout <<"94->-15->[0]->getPDG: "<< tauonPlus->getDaughters()[0]->getPDG()<<endl;
						if (tauonMinus->getDaughters()[0]->getPDG() == 15)
							tauonMinus = tauonMinus->getDaughters()[0];
						if (tauonPlus->getDaughters()[0]->getPDG() == -15)
							tauonPlus = tauonPlus->getDaughters()[0];

					}

					if (tauonPlusTemp->getDaughters().size() == 1) {
						tauonMinus =
								tauonPlusTemp->getDaughters()[0]->getDaughters()[0];
						tauonPlus =
								tauonPlusTemp->getDaughters()[0]->getDaughters()[1];
						if (tauonMinus->getDaughters()[0]->getPDG() == 15)
							tauonMinus = tauonMinus->getDaughters()[0];
						if (tauonPlus->getDaughters()[0]->getPDG() == -15)
							tauonPlus = tauonPlus->getDaughters()[0];
					}

// -------------------------------------------------------------------------------------

					const double *pTauMinus = tauonMinus->getMomentum();
					double eTauMinus = tauonMinus->getEnergy();

					Tauon_minus.SetPxPyPzE(pTauMinus[0], pTauMinus[1],
							pTauMinus[2], eTauMinus);

					TLorentzVector TLVtauon_potomak_od_tauona_minus,
							TLVfoton_potomak_od_tauona_minus;
					TLorentzVector TLVtauon_potomak_od_tauona_plus,
							TLVfoton_potomak_od_tauona_plus;

					for (int i = 0; i < (int) tauonMinus->getDaughters().size();
							i++) {
						Int_t PDG = tauonMinus->getDaughters()[i]->getPDG();

						const double *p =
								tauonMinus->getDaughters()[i]->getMomentum();
						double e = tauonMinus->getDaughters()[i]->getEnergy();

						TLorentzVector temp;
						temp.SetPxPyPzE(p[0], p[1], p[2], e);

						if (PDG == 16) {
							TLVtauonskiNeutrino = temp;
						}

						if (PDG == -211) {
							B_pionMinus = true;

							TLVpionMinus = temp;

							N_pion_minus_od_tauona++;
						}

						if (PDG == 15) {
							B_tauon_potomak_od_tauona_minus = true;
							TLVtauon_potomak_od_tauona_minus = temp;
							//		TLVpionMinus = temp;
							N_tau_minus_from_photon++;
						}

						if (PDG == 22) {
							B_foton_potomak_od_tauona_minus = true;
							TLVfoton_potomak_od_tauona_minus = temp;
						}

					}

					if (B_pionMinus) {
						TauonMinusko = TLVpionMinus + TLVtauonskiNeutrino;

					}

					if (B_tauon_potomak_od_tauona_minus
							&& B_foton_potomak_od_tauona_minus) {
						TVector3 tauon_potomak_od_tauona_minus3(
								TLVtauon_potomak_od_tauona_minus.Px(),
								TLVtauon_potomak_od_tauona_minus.Py(),
								TLVtauon_potomak_od_tauona_minus.Pz());
						TVector3 foton_potomak_od_tauona_minus3(
								TLVfoton_potomak_od_tauona_minus.Px(),
								TLVfoton_potomak_od_tauona_minus.Py(),
								TLVfoton_potomak_od_tauona_minus.Pz());

					}

// -------------------------------------------------------------------------------------

					const double *pTauPlus = tauonPlus->getMomentum();
					double eTauPlus = tauonPlus->getEnergy();

					Tauon_plus.SetPxPyPzE(pTauPlus[0], pTauPlus[1], pTauPlus[2],
							eTauPlus);

					for (int i = 0; i < (int) tauonPlus->getDaughters().size();
							i++) {
						Int_t PDG = tauonPlus->getDaughters()[i]->getPDG();

						const double *p =
								tauonPlus->getDaughters()[i]->getMomentum();
						double e = tauonPlus->getDaughters()[i]->getEnergy();

						TLorentzVector temp;
						temp.SetPxPyPzE(p[0], p[1], p[2], e);

						if (PDG == -16) {
							TLVtauonskiAntineutrino = temp;
						}

						if (PDG == 211) {
							B_pionPlus = true;

							TLVpionPlus = temp;

							N_pion_plus_od_tauona++;
						}

						if (PDG == -15) {
							B_tauon_potomak_od_tauona_plus = true;
							TLVtauon_potomak_od_tauona_plus = temp;
							//		TLVpionPlus = temp;
							N_tau_plus_from_photon++;
						}

						if (PDG == 22) {
							B_foton_potomak_od_tauona_plus = true;
							TLVfoton_potomak_od_tauona_plus = temp;
						}
					}

					if (B_pionPlus) {
						TauonPlusko = TLVpionPlus + TLVtauonskiAntineutrino;

					}

					if (B_tauon_potomak_od_tauona_plus
							&& B_foton_potomak_od_tauona_plus) {
						TVector3 tauon_potomak_od_tauona_plus3(
								TLVtauon_potomak_od_tauona_plus.Px(),
								TLVtauon_potomak_od_tauona_plus.Py(),
								TLVtauon_potomak_od_tauona_plus.Pz());
						TVector3 foton_potomak_od_tauona_plus3(
								TLVfoton_potomak_od_tauona_plus.Px(),
								TLVfoton_potomak_od_tauona_plus.Py(),
								TLVfoton_potomak_od_tauona_plus.Pz());

					}

					// histogram_Higs_Tauoni.Fill( TauonMinusko.M() + TauonPlusko.M() );
					// histogram_Higs_Pioni.Fill( TLVpionMinus.M() + TLVtauonskiNeutrino.M() + TLVpionPlus.M() + TLVtauonskiAntineutrino.M() );

				} // Kraj IF za tauone 15 ili -15

			} // Kraj FOR petlje po broju elemenata
			if (B_tauon_potomak_od_tauona_minus
					&& B_tauon_potomak_od_tauona_plus)
				N_evt_with_photon++;

			if (TLVpionPlus.P() > 0 && TLVpionMinus.P() > 0) {
				N_signjala++;

				P_Pion_minus = TLVpionMinus.P();
				P_Pion_plus = TLVpionPlus.P();

				TLorentzVector TLVHigs = Tauon_minus + Tauon_plus;
				TLorentzVector Tminus = Tauon_minus;
				TLorentzVector Tplus = Tauon_plus;



				TVector3 Higs3(TLVHigs.Px(), TLVHigs.Py(), TLVHigs.Pz());

				labHigs_Mass = TLVHigs.M();
				labHigs_Polar = TLVHigs.Theta();
				labHigs_Azi = TLVHigs.Phi();
				labHigs_E = TLVHigs.E();
				labHigs_pt = TLVHigs.Pt();
				labHigs_X = TLVHigs.X();
				labHigs_Y = TLVHigs.Y();
				labHigs_Z = TLVHigs.Z();

				labTauM_Mass = Tauon_minus.M();
				labTauM_Polar = Tauon_minus.Theta();
				labTauM_Azi = Tauon_minus.Phi();
				labTauM_E = Tauon_minus.E();
				labTauM_pt = Tauon_minus.Pt();
				labTauM_X = Tauon_minus.X();
				labTauM_Y = Tauon_minus.Y();
				labTauM_Z = Tauon_minus.Z();

				labTauP_Mass = Tauon_plus.M();
				labTauP_Polar = Tauon_plus.Theta();
				labTauP_Azi = Tauon_plus.Phi();
				labTauP_E = Tauon_plus.E();
				labTauP_pt = Tauon_plus.Pt();
				labTauP_X = Tauon_plus.X();
				labTauP_Y = Tauon_plus.Y();
				labTauP_Z = Tauon_plus.Z();

				labpionM_Mass = TLVpionMinus.M();
				labpionM_Polar = TLVpionMinus.Theta();
				labpionM_Azi = TLVpionMinus.Phi();
				labpionM_E = TLVpionMinus.E();
				labpionM_pt = TLVpionMinus.Pt();
				labpionM_X = TLVpionMinus.X();
				labpionM_Y = TLVpionMinus.Y();
				labpionM_Z = TLVpionMinus.Z();

				labpionP_Mass = TLVpionMinus.M();
				labpionP_Polar = TLVpionMinus.Theta();
				labpionP_Azi = TLVpionMinus.P();
				labpionP_E = TLVpionPlus.E();
				labpionP_pt = TLVpionMinus.Pt();
				labpionP_X = TLVpionMinus.X();
				labpionP_Y = TLVpionMinus.Y();
				labpionP_Z = TLVpionMinus.Z();

				//pocinje boost ***************************************************************************

				TVector3 BoostToOriginal = -(Higs_pocetni.BoostVector()); // KG  vektor boost u Higsov sistem reference
				Tminus.Boost(BoostToOriginal);
				Tplus.Boost(BoostToOriginal);

				TVector3 tauPlusBoost(Tplus.Px(), Tplus.Py(),Tplus.Pz());
				TVector3 tauMinusBoost(Tminus.Px(), Tminus.Py(),Tminus.Pz());

				angleTausOriginal = tauMinusBoost.Angle(tauPlusBoost);





				TVector3 BoostToHiggs = -(TLVHigs.BoostVector()); // KG  vektor boost u Higsov sistem reference

				TLVHigs.Boost(BoostToHiggs);
				TLVpionMinus.Boost(BoostToHiggs); //
				TLVpionPlus.Boost(BoostToHiggs);
				Tauon_minus.Boost(BoostToHiggs);
				Tauon_plus.Boost(BoostToHiggs);
				//	cout<<"Higgs 4V nakon boost u Higgs sistem: "<<endl;	TLVHigs.Print();

				//

				HS_Higs_Mass = TLVHigs.M();
				HS_Higs_Polar = TLVHigs.Theta();
				HS_Higs_Azi = TLVHigs.Phi();
				HS_Higs_E = TLVHigs.E();
				HS_Higs_pt = TLVHigs.Pt();
				HS_Higs_X = TLVHigs.X();
				HS_Higs_Y = TLVHigs.Y();
				HS_Higs_Z = TLVHigs.Z();

				HS_TauM_Mass = Tauon_minus.M();
				HS_TauM_Polar = Tauon_minus.Theta();
				HS_TauM_Azi = Tauon_minus.Phi();
				HS_TauM_E = Tauon_minus.E();
				HS_TauM_pt = Tauon_minus.Pt();
				HS_TauM_X = Tauon_minus.X();
				HS_TauM_Y = Tauon_minus.Y();
				HS_TauM_Z = Tauon_minus.Z();

				HS_TauP_Mass = Tauon_plus.M();
				HS_TauP_Polar = Tauon_plus.Theta();
				HS_TauP_Azi = Tauon_plus.Phi();
				HS_TauP_E = Tauon_plus.E();
				HS_TauP_pt = Tauon_plus.Pt();
				HS_TauP_X = Tauon_plus.X();
				HS_TauP_Y = Tauon_plus.Y();
				HS_TauP_Z = Tauon_plus.Z();

				HS_pionM_Mass = TLVpionMinus.M();
				HS_pionM_Polar = TLVpionMinus.Theta();
				HS_pionM_Azi = TLVpionMinus.Phi();
				HS_pionM_E = TLVpionMinus.E();
				HS_pionM_pt = TLVpionMinus.Pt();
				HS_pionM_X = TLVpionMinus.X();
				HS_pionM_Y = TLVpionMinus.Y();
				HS_pionM_Z = TLVpionMinus.Z();

				HS_pionP_Mass = TLVpionPlus.M();
				HS_pionP_Polar = TLVpionPlus.Theta();
				HS_pionP_Azi = TLVpionPlus.P();
				HS_pionP_E = TLVpionPlus.E();
				HS_pionP_pt = TLVpionPlus.Pt();
				HS_pionP_X = TLVpionPlus.X();
				HS_pionP_Y = TLVpionPlus.Y();
				HS_pionP_Z = TLVpionPlus.Z();
			//	cout << "Higs vektor: "   << endl;
		//		TLVHigs.Print();


				/*	if (Tauon_minus.Px() > 0.00001 ) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;
				 if ( Tauon_minus.Py() > 0.00001 ) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;
				 if (Tauon_minus.Pz() > 0.00001) cout << "Px = " <<Tauon_minus.Px()<<", Py= "<<Tauon_minus.Py()<<", pz = "<<Tauon_minus.Pz()<<endl;*/

				TVector3 pionMinus3(TLVpionMinus.Px(), TLVpionMinus.Py(),
						TLVpionMinus.Pz());
				TVector3 pionPlus3(TLVpionPlus.Px(), TLVpionPlus.Py(),
						TLVpionPlus.Pz());
				TVector3 tauMinus3(Tauon_minus.Px(), Tauon_minus.Py(),
						Tauon_minus.Pz());
				TVector3 tauPlus3(Tauon_plus.Px(), Tauon_plus.Py(),
						Tauon_plus.Pz());
				Higs3.SetXYZ(TLVHigs.Px(), TLVHigs.Py(), TLVHigs.Pz());

			angleTaus = tauPlus3.Angle(tauMinus3);
				//	cout <<"ugao izmedju tau čestica nakon boost a u Higsov sistem: "<< Tauon_minus.Angle(tauPlus3)<<endl;

				//	cout<<"Higgs nakon boostovanja u Higgs sistem: "<<endl;	Higs3.Print();
				cout	<< "ugao izmedju tau i piona pre rotacije ali nakon boosta: "
						<< pionMinus3.Angle(tauMinus3) << endl;

				TVector3 x(1, 0, 0);
				TVector3 y(0, 1, 0);
				TVector3 z(0, 0, 1);

				TVector3 directionMinus = tauMinus3;
				directionMinus *= 1. / directionMinus.Mag(); // KG pravac Tau minus je z osa

						/*			TVector3 testy = z.Cross(directionMinus);
						 double testAngle = TMath::ACos(z.Angle(directionMinus));
						 TRotation rTest;
						 rTest.AngleAxis(testAngle, directionMinus);

						 tauMinus3.Transform(rTest);

						 cout <<"tau minus ose  x= "<<tauMinus3.X()<<", y = " <<tauMinus3.Y()<<", z = "<<tauMinus3.Z()<<endl;*/

				//		TVector3 rotatedYAxisMinus = x.Cross(directionMinus);// test
				//	TVector3 rotatedYAxisMinus = directionMinus.Cross(x);
				TVector3 rotatedXAxisMinus = y.Cross(directionMinus);

				//	rotatedYAxisMinus *= 1.0/rotatedYAxisMinus.Mag();

			/*	TRotation rMinus;
				rMinus.SetZAxis(directionMinus, rotatedXAxisMinus);
				rMinus.Invert();*/

				TRotation r;
				TVector3 axis = z.Cross(tauMinus3);
				double angle_rot = tauMinus3.Theta();
				r.Rotate(angle_rot, axis);

				//	Higs3.Transform(rMinus);
				pionMinus3.Transform(r);
				pionPlus3.Transform(r);
				tauMinus3.Transform(r);
				tauPlus3.Transform(r);

				//	cout<<rotatedYAxisMinus.Z()<<endl;
				cout << "ugao izmedju tau i piona nakon rotacije: "
						<< pionMinus3.Angle(tauMinus3) << endl;

				// rotirati x y osu matricom rotacije koja je napravljana

				//cout<<"Higgs nakon rotacije: "<<endl;	Higs3.Print();



				cout << "tau minus: " << endl;
				tauMinus3.Print();
				cout << "tau plus: " << endl;
				tauPlus3.Print();
				/*		cout<<"pion minus: "<<endl;	pionMinus3.Print();
				 cout<<"pion plus: "<<endl;	pionPlus3.Print();*/
				cout << "____________________________________________" << endl;

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
				Double_t pMinusPhi = pionMinus3.Phi();
				if (pMinusPhi < 0)
					pMinusPhi = 2 * M_PI - abs(pMinusPhi);

				Double_t pPlusPhi = pionPlus3.Phi();			//
				if (pPlusPhi < 0)
					pPlusPhi = 2 * M_PI - abs(pPlusPhi);

				Double_t pMinusPhi1 = pionMinus3.Phi();
				Double_t pPlusPhi1 = pionPlus3.Phi();			//
				deltaPfi1 = pPlusPhi1 - pMinusPhi1;

				//	if (Higs_tauon_minus3.Phi() < 0)  Higs_tauon_minus3.Print();

				pionMinus_azimut = pMinusPhi;
				pionMinus_polar = pionMinus3.Theta();
				pionPlus_azimut = pPlusPhi;
				pionPlus_polar = pionPlus3.Theta();
				higs_azimut = Higs3.Phi();
				higs_polar = Higs3.Theta();
				tauMinus_polar = tauMinus3.Theta();
				tauPlus_polar = tauPlus3.Theta();
				tauMinus_azimut = tauMinus3.Phi();
				tauPlus_azimut = tauPlus3.Phi();
				deltaPfi = pPlusPhi - pMinusPhi;

				tauMinus3_x = tauMinus3.X();
				tauMinus3_y = tauMinus3.Y();
				tauMinus3_z = tauMinus3.Z();
				tauPlus3_x = tauPlus3.X();
				tauPlus3_y = tauPlus3.Y();
				tauPlus3_z = tauPlus3.Z();
				pionPlus3_x = pionPlus3.X();
				pionPlus3_y = pionPlus3.Y();
				pionPlus3_z = pionPlus3.Z();
				pionMinus3_x = pionMinus3.X();
				pionMinus3_y = pionMinus3.Y();
				pionMinus3_z = pionMinus3.Z();
				Higgs3_x = Higs3.X();
				Higgs3_y = Higs3.Y();
				Higgs3_z = Higs3.Z();

				distroTree.Fill();
				polarimeterTree.Fill();
				labTree.Fill();
				HSTree.Fill();
			}

		} // Kraj WHILE petlje

		JinxTree.Fill();

		lcReader->close();

	} // Kraj FOR petlje koja iščitava .slcio fajlove

	cout << "Ukupan broj signjala " << N_signjala << endl;

	cout << "Ukupan broj piona minus od tauona " << N_pion_minus_od_tauona
			<< endl;
	cout << "Ukupan broj piona plus od tauona " << N_pion_plus_od_tauona
			<< endl;
	//N_tau_minus_from_photon
//	cout << "Ukupan broj tau minus sa fotonom: " << N_tau_minus_from_photon << endl;
//	cout << "Ukupan broj tau plus sa fotonom: " << N_tau_plus_from_photon << endl;
//	cout << "Ukupan broj dogadjaja sa fotonom: " << N_evt_with_photon << endl;

//	JinxTree.Fill();

	TString tfName(rfn);
	if (!tfName.EndsWith(".root"))
		tfName.Append(".root");
	TFile rootFile(tfName.Data(), "RECREATE");

	//JinxTree.Write();
	polarimeterTree.Write();
	distroTree.Write();
	labTree.Write();
	HSTree.Write();

	// evtTree.Write();
	rootFile.Write();
	rootFile.Close();

	return 0;

} // Kraj funkcije

// ------------------------------------------------------------------

Int_t main(int argc, char* argv[]) {
	Int_t iarg = 1;

	UInt_t nFirstJob = 1;

	if (argc > iarg) {
		nFirstJob = atoi(argv[iarg]);
		iarg++;
	}

	UInt_t nLastJob = 50;

	if (argc > iarg) {
		nLastJob = atoi(argv[iarg]);
		iarg++;
	}

	TString fName = "test_";
	if (argc > iarg)
		fName = argv[iarg]; //iarg++;

	TString rfName = "Polarimetar_Jinx_MC.root";
	if (argc > iarg)
		rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
