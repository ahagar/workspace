// ee -> H -> ee
// Napravljen dana: 26.03.2018.
// Autor: Natasa
// Modifikovao MIRKO dana 7.01.2021.

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

#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <array>

using namespace std;

//     --------------------    RecoMCTruth Link    --------------------

Int_t slcio2appTree(UInt_t nPrviFajl, UInt_t nZadnjiFajl, const char* fn, const char* rfn)
{
	#ifdef __CINT__
		gSystem -> Load("${LCIO}/lib/liblcio.so");
		gSystem -> Load("${LCIO}/lib/liblcioDict.so");
	#endif

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	TTree leptonTree ("leptonTree", "Generator particle tree");
	TTree ptTree ("ptTree", "Generator particle tree");
	TTree pfoTree ("pfoTree", "Generator particle tree");
	TTree mcTree ("mcTree", "Generator particle tree"); // definicija drveta gde upisujemo histograme
	TTree elektronskoTree("elektronskoTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za sve elektrone i pozitrone
	TTree elektronsko6i7Tree("elektronsko6i7Tree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami samo za elektron 6 i pozitron 7
	TTree elektronskoRekonstrTree("elektronskoRekonstrTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za rekonstruisane elektrone i pozitrone
//	TTree pozitronskoRekonstrTree("pozitronskoRekonstrTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za rekonstruisane pozitrone
//	elektronskoRekonstrTree.AddFriend("pozitronskoRekonstrTree", "");
	leptonTree.AddFriend("ptTree", "");
	leptonTree.AddFriend("pfoTree", "");
	TTree testTree ("testTree", "testTree"); // definicija drveta gde upisujemo histograme


	Float_t odnos_ecal_ukupnical_leptona;
	Float_t odnos_ecal_ukupnical_leptona_cutted = 0;
	Float_t d0_recparticle, z0_recparticle, r03d_recparticle;
	Float_t d0_recparticle_cutted = 0, z0_recparticle_cutted = 0, r03d_recparticle_cutted = 0;
	Float_t Energija_konusa = 0, Energija_leptonskog_traga = 0;
	Float_t Energija_konusa_cutted = 0, Energija_leptonskog_traga_cutted = 0;
	Float_t lep_pfoj_angle;
	Float_t pt_leptona, p_leptona;
	Float_t pt_leptona_cutted = 0;
	Float_t n_lep_evt, n_kriva;
	Float_t n_lep_evt_cutted = 0;
	Float_t n_lep_evt_d0_cutted = 0, n_lep_evt_z0_cutted = 0, n_lep_evt_r03d_cutted = 0;
	Float_t n_lep_evt_ratio_cutted = 0, n_lep_evt_pt_cutted = 0;
	Float_t brojac_prezivelih_leptona = 0;
	Float_t n_track, n_track_pt;
//	Float_t n_ro_track;
//	Float_t n_do_track;
//	Float_t n_zo_track;
//	Float_t n_rdzo_track;
	Float_t n_svipfo;
	Float_t pt_pfo, p_pfo;
	Float_t pt_pfo_cutted = 0;
	Float_t odnos_ecal_ukupnical_pfo;
	Float_t odnos_ecal_ukupnical_pfo_cutted = 0;
	Float_t d0_pfo, z0_pfo, r03d_pfo;
	Float_t d0_pfo_cutted = 0, z0_pfo_cutted = 0, r03d_pfo_cutted = 0;
	Float_t Energija_konusa_pfo, pfoTrackEn;
	Float_t pfoTrackEn_cutted = 0;
	Float_t coneAngle_pfo;

	// Promenljive za KRAJNJI elektron/pozitron u cesticnom lancu
	Float_t brojac_preostalih_KG_elektrona = 0;
	Float_t pt_KG_elektrona = 0, p_KG_elektrona = 0;
	Float_t Energija_KG_elektrona = 0;
	Float_t Energija_KG_elektrona_cutted = 0;
	Float_t Theta_KG_elektrona = 0;
	Float_t odnos_Ecal_TOTALcal_KG_elektrona = 0;
	Float_t d0_KG_elektrona = 0, z0_KG_elektrona = 0, r0_KG_elektrona = 0;
	Float_t Energija_konusa_KG_elektrona = 0;

	// Promenljive za POCETNI elektron 6
//	Float_t brojac_preostalih_KG_elektrona_6 = 0;
	Float_t pt_KG_elektrona_6 = 0, p_KG_elektrona_6 = 0;
	Float_t Energija_KG_elektrona_6 = 0;
	Float_t Theta_KG_elektrona_6 = 0;
//	Float_t odnos_Ecal_TOTALcal_KG_elektrona_6 = 0;
	Float_t d0_KG_elektrona_6 = 0, z0_KG_elektrona_6 = 0, r0_KG_elektrona_6 = 0;
	Float_t mee_gen= 0;

	// Promenljive za POCETNI pozitron 7
//	Float_t brojac_preostalih_KG_pozitrona_7 = 0;
	Float_t pt_KG_pozitrona_7 = 0, p_KG_pozitrona_7 = 0;
	Float_t Energija_KG_pozitrona_7 = 0;
	Float_t Theta_KG_pozitrona_7 = 0;
//	Float_t odnos_Ecal_TOTALcal_KG_pozitrona_7 = 0;
	Float_t d0_KG_pozitrona_7 = 0, z0_KG_pozitrona_7 = 0, r0_KG_pozitrona_7 = 0;
//	Float_t Energija_konusa_KG_pozitrona_7 = 0;

	// Promenljive za REKONSTRUISANI elektron
//	Float_t brojac_preostalih_KG_elektrona = 0;
	Float_t pt_KG_elektrona_rekonstr = 0, p_KG_elektrona_rekonstr = 0;
	Float_t Energija_KG_elektrona_rekonstr = 0;
	Float_t Energija_KG_elektrona_rekonstr_cutted = 0;
	Float_t Theta_KG_elektrona_rekonstr = 0;
	Float_t odnos_Ecal_TOTALcal_KG_elektrona_rekonstr = 0;
	Float_t d0_KG_elektrona_rekonstr = 0, z0_KG_elektrona_rekonstr = 0, r0_KG_elektrona_rekonstr = 0;
	Float_t Energija_konusa_KG_elektrona_rekonstr = 0;
	Float_t Brojac_za_dvojku = 0;

	// Promenljive za REKONSTRUISANI pozitron
//	Float_t brojac_preostalih_KG_pozitrona = 0;
	Float_t pt_KG_pozitrona_rekonstr = 0, p_KG_pozitrona_rekonstr = 0;
	Float_t Energija_KG_pozitrona_rekonstr = 0;
	Float_t Energija_KG_pozitrona_rekonstr_cutted = 0;
	Float_t Theta_KG_pozitrona_rekonstr = 0;
	Float_t odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr = 0;
	Float_t d0_KG_pozitrona_rekonstr = 0, z0_KG_pozitrona_rekonstr = 0, r0_KG_pozitrona_rekonstr = 0;
	Float_t Energija_konusa_KG_pozitrona_rekonstr = 0;
//	Float_t Brojac_za_dvojku = 0;

	// Promenljive za PFO rekonstruisanih elektrona/pozitrona
//	Float_t brojac_preostalih_KG_elektrona = 0;
	Float_t pt_KG_pfo = 0, p_KG_pfo = 0;
	Float_t Energija_KG_pfo = 0;
	Float_t Theta_KG_pfo = 0;
	Float_t odnos_Ecal_TOTALcal_KG_pfo = 0;
	Float_t d0_KG_pfo = 0, z0_KG_pfo = 0, r0_KG_pfo = 0;
	Float_t Energija_konusa_KG_pfo = 0;

	Float_t mHg = 0, mZ1g = 0, mee = 0, mZ2g =0, Eee = 0, thetaZ1 = 0, thetaZ2 = 0, pxZ1=0, pxZ2=0, pyZ1=0, pyZ2=0, pzZ1=0, pzZ2=0;
	Float_t thetaH = 0, pxH = 0, pyH = 0 , pzH = 0, Ees = 0 , Eps = 0, Ee6 = 0, Ep7 = 0, Eestarac = 0, Epstarac = 0;

	// Grane drveta za KRAJNJI generisani elektron/pozitron
	elektronskoTree.Branch("pt_KG_elektrona", &pt_KG_elektrona, "pt_KG_elektrona");
	elektronskoTree.Branch("p_KG_elektrona", &p_KG_elektrona, "p_KG_elektrona");
	elektronskoTree.Branch("Theta_KG_elektrona", &Theta_KG_elektrona, "Theta_KG_elektrona");
	elektronskoTree.Branch("Energija_KG_elektrona", &Energija_KG_elektrona, "Energija_KG_elektrona");
	elektronskoTree.Branch("odnos_Ecal_TOTALcal_KG_elektrona", &odnos_Ecal_TOTALcal_KG_elektrona, "odnos_Ecal_TOTALcal_KG_elektrona");

	testTree.Branch("mee", &mee, "mee");
	testTree.Branch("Eee", &Eee, "Eee");
	testTree.Branch("mHg", &mHg, "mHg");
	testTree.Branch("mZ1g", &mZ1g, "mZ1g");
	testTree.Branch("mZ2g", &mZ2g, "mZ2g");
	testTree.Branch("thetaZ1", &thetaZ1, "thetaZ1");
	testTree.Branch("thetaZ2", &thetaZ2, "thetaZ2");
	testTree.Branch("pxZ1", &pxZ1, "pxZ1");
	testTree.Branch("pyZ1", &pyZ1, "pyZ1");
	testTree.Branch("pzZ1", &pzZ1, "pzZ1");
	testTree.Branch("pxZ2", &pxZ2, "pxZ2");
	testTree.Branch("pyZ2", &pyZ2, "pyZ2");
	testTree.Branch("pzZ2", &pzZ2, "pzZ2");
	testTree.Branch("pxH", &pxH, "pxH");
	testTree.Branch("pyH", &pyH, "pyH");
	testTree.Branch("pzH", &pzH, "pzH");
	testTree.Branch("thetaH", &thetaH, "thetaH");
	testTree.Branch("Ees", &Ees, "Ees");
	testTree.Branch("Eps", &Eps, "Eps");
	testTree.Branch("Ee6", &Ee6, "Ee6");
	testTree.Branch("Ep7", &Ep7, "Ep7");
	testTree.Branch("Eestarac", &Eestarac, "Eestarac");
	testTree.Branch("Epstarac", &Epstarac, "Epstarac");





	// Grane drveta za POCETNI generisani elektron 6
	elektronsko6i7Tree.Branch("mee_gen", &mee_gen, "mee_gen");
	elektronsko6i7Tree.Branch("pt_KG_elektrona_6", &pt_KG_elektrona_6, "pt_KG_elektrona_6");
	elektronsko6i7Tree.Branch("p_KG_elektrona_6", &p_KG_elektrona_6, "p_KG_elektrona_6");
	elektronsko6i7Tree.Branch("Theta_KG_elektrona_6", &Theta_KG_elektrona_6, "Theta_KG_elektrona_6");
//	elektronsko6i7Tree.Branch("odnos_Ecal_TOTALcal_KG_elektrona_6", &odnos_Ecal_TOTALcal_KG_elektrona_6, "odnos_Ecal_TOTALcal_KG_elektrona_6");
	elektronsko6i7Tree.Branch("d0_KG_elektrona_6", &d0_KG_elektrona_6, "d0_KG_elektrona_6");
	elektronsko6i7Tree.Branch("z0_KG_elektrona_6", &z0_KG_elektrona_6, "z0_KG_elektrona_6");
	elektronsko6i7Tree.Branch("r0_KG_elektrona_6", &r0_KG_elektrona_6, "r0_KG_elektrona_6");
	elektronsko6i7Tree.Branch("Energija_KG_elektrona_6", &Energija_KG_elektrona_6, "Energija_KG_elektrona_6");
//	elektronsko6i7Tree.Branch("Energija_konusa_KG_elektrona_6", &Energija_konusa_KG_elektrona_6, "Energija_konusa_KG_elektrona_6");

	// Grane drveta za POCETNI generisani pozitron 7
	elektronsko6i7Tree.Branch("pt_KG_pozitrona_7", &pt_KG_pozitrona_7, "pt_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("p_KG_pozitrona_7", &p_KG_pozitrona_7, "p_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("Theta_KG_pozitrona_7", &Theta_KG_pozitrona_7, "Theta_KG_pozitrona_7");
//	elektronsko6i7Tree.Branch("odnos_Ecal_TOTALcal_KG_pozitrona_7", &odnos_Ecal_TOTALcal_KG_pozitrona_7, "odnos_Ecal_TOTALcal_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("d0_KG_pozitrona_7", &d0_KG_pozitrona_7, "d0_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("z0_KG_pozitrona_7", &z0_KG_pozitrona_7, "z0_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("r0_KG_pozitrona_7", &r0_KG_pozitrona_7, "r0_KG_pozitrona_7");
	elektronsko6i7Tree.Branch("Energija_KG_pozitrona_7", &Energija_KG_pozitrona_7, "Energija_KG_pozitrona_7");
//	elektronsko6i7Tree.Branch("Energija_konusa_KG_pozitrona_7", &Energija_konusa_KG_pozitrona_7, "Energija_konusa_KG_pozitrona_7");

	// Grane drveta za REKONSTRUISANI elektron
	elektronskoRekonstrTree.Branch("pt_KG_elektrona_rekonstr", &pt_KG_elektrona_rekonstr, "pt_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("p_KG_elektrona_rekonstr", &p_KG_elektrona_rekonstr, "p_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("Theta_KG_elektrona_rekonstr", &Theta_KG_elektrona_rekonstr, "Theta_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("d0_KG_elektrona_rekonstr", &d0_KG_elektrona_rekonstr, "d0_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("z0_KG_elektrona_rekonstr", &z0_KG_elektrona_rekonstr, "z0_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("r0_KG_elektrona_rekonstr", &r0_KG_elektrona_rekonstr, "r0_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("Energija_KG_elektrona_rekonstr", &Energija_KG_elektrona_rekonstr, "Energija_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("Energija_konusa_KG_elektrona_rekonstr", &Energija_konusa_KG_elektrona_rekonstr, "Energija_konusa_KG_elektrona_rekonstr");
	elektronskoRekonstrTree.Branch("odnos_Ecal_TOTALcal_KG_elektrona_rekonstr", &odnos_Ecal_TOTALcal_KG_elektrona_rekonstr, "odnos_Ecal_TOTALcal_KG_elektrona_rekonstr");

//	elektronskoRekonstrTree.Branch("Brojac_za_dvojku",&Brojac_za_dvojku,"Brojac_za_dvojku");

	// Grane drveta za REKONSTRUISANI pozitron
	elektronskoRekonstrTree.Branch("pt_KG_pozitrona_rekonstr", &pt_KG_pozitrona_rekonstr, "pt_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("p_KG_pozitrona_rekonstr", &p_KG_pozitrona_rekonstr, "p_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("Theta_KG_pozitrona_rekonstr", &Theta_KG_pozitrona_rekonstr, "Theta_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("d0_KG_pozitrona_rekonstr", &d0_KG_pozitrona_rekonstr, "d0_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("z0_KG_pozitrona_rekonstr", &z0_KG_pozitrona_rekonstr, "z0_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("r0_KG_pozitrona_rekonstr", &r0_KG_pozitrona_rekonstr, "r0_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("Energija_KG_pozitrona_rekonstr", &Energija_KG_pozitrona_rekonstr, "Energija_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("Energija_konusa_KG_pozitrona_rekonstr", &Energija_konusa_KG_pozitrona_rekonstr, "Energija_konusa_KG_pozitrona_rekonstr");
	elektronskoRekonstrTree.Branch("odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr", &odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr, "odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr");

//	elektronskoRekonstrTree.Branch("Brojac_za_dvojku",&Brojac_za_dvojku,"Brojac_za_dvojku");

	// Grane drveta za sve PFOs
	pfoTree.Branch("pt_KG_pfo", &pt_KG_pfo, "pt_KG_pfo");
//	pfoTree.Branch("pt_KG_pfo_cutted", &pt_KG_pfo_cutted, "pt_KG_pfo_cutted");
//	pfoTree.Branch("p_KG_pfo", &p_KG_pfo, "p_KG_pfo");
	pfoTree.Branch("Theta_KG_pfo", &Theta_KG_pfo, "Theta_KG_pfo");
	pfoTree.Branch("d0_KG_pfo", &d0_KG_pfo, "d0_KG_pfo");
//	pfoTree.Branch("d0_pfo_cutted", &d0_pfo_cutted, "d0_pfo_cutted");
	pfoTree.Branch("z0_KG_pfo", &z0_KG_pfo, "z0_KG_pfo");
//	pfoTree.Branch("z0_pfo_cutted", &z0_pfo_cutted, "z0_pfo_cutted");
	pfoTree.Branch("r0_KG_pfo", &r0_KG_pfo, "r0_KG_pfo");
//	pfoTree.Branch("r0_KG_pfo_cutted", &r0_KG_pfo_cutted, "r0_KG_pfo_cutted");
	pfoTree.Branch("Energija_KG_pfo", &Energija_KG_pfo, "Energija_KG_pfo");
//	pfoTree.Branch("Energija_KG_pfo_cutted", &Energija_KG_pfo_cutted, "Energija_KG_pfo_cutted");
//	pfoTree.Branch("Energija_konusa_KG_pfo", &Energija_konusa_KG_pfo, "Energija_konusa_KG_pfo");
	pfoTree.Branch("odnos_Ecal_TOTALcal_KG_pfo", &odnos_Ecal_TOTALcal_KG_pfo, "odnos_Ecal_TOTALcal_KG_pfo");
//	pfoTree.Branch("odnos_Ecal_TOTALcal_KG_pfo_cutted", &odnos_Ecal_TOTALcal_KG_pfo_cutted, "odnos_Ecal_TOTALcal_KG_pfo_cutted");

/*	leptonTree.Branch("Energija_konusa", &Energija_konusa, "Energija_konusa");
	leptonTree.Branch("Energija_konusa_cutted", &Energija_konusa_cutted, "Energija_konusa_cutted");
	leptonTree.Branch("Energija_leptonskog_traga", &Energija_leptonskog_traga, "Energija_leptonskog_traga");
	leptonTree.Branch("Energija_leptonskog_traga_cutted", &Energija_leptonskog_traga_cutted, "Energija_leptonskog_traga_cutted");
	leptonTree.Branch("odnos_ecal_ukupnical_leptona", &odnos_ecal_ukupnical_leptona, "odnos_ecal_ukupnical_leptona");
	leptonTree.Branch("odnos_ecal_ukupnical_leptona_cutted", &odnos_ecal_ukupnical_leptona_cutted, "odnos_ecal_ukupnical_leptona_cutted");
	leptonTree.Branch("r03d_recparticle", &r03d_recparticle, "r03d_recparticle");
	leptonTree.Branch("r03d_recparticle_cutted", &r03d_recparticle_cutted, "r03d_recparticle_cutted");
	leptonTree.Branch("d0_recparticle", &d0_recparticle, "d0_recparticle");
	leptonTree.Branch("d0_recparticle_cutted", &d0_recparticle_cutted, "d0_recparticle_cutted");
	leptonTree.Branch("z0_recparticle", &z0_recparticle, "z0_recparticle");
	leptonTree.Branch("z0_recparticle_cutted", &z0_recparticle_cutted, "z0_recparticle_cutted");
	leptonTree.Branch("lep_pfoj_angle", &lep_pfoj_angle, "lep_pfoj_angle");
	leptonTree.Branch("pt_leptona", &pt_leptona, "pt_leptona");
	leptonTree.Branch("pt_leptona_cutted", &pt_leptona_cutted, "pt_leptona_cutted");
	leptonTree.Branch("p_leptona", &p_leptona, "p_leptona");

	ptTree.Branch("n_lep_evt", &n_lep_evt, "n_lep_evt");
	ptTree.Branch("n_kriva", &n_kriva, "n_kriva");
	ptTree.Branch("n_svipfo", &n_svipfo, "n_svipfo");      */

	TH1F histo_niso_lep ("histo_niso_lep", " ; n_{izolep}", 50, 0, 50);
	TH1F histo_mcvect_size ("histo_mcvect_size", " ; size", 100, 0, 10);
	TH2F histoTrackEnergija_konusa ("histoTrackEnergija_konusa", " ; #it{E}_{track} (GeV); #it{E}_{cone} (GeV)", 100, 0, 750, 100, 0, 200);
	TH2F histoTrackEnergija_konusa_cutted ("histoTrackEnergija_konusa_cutted", " ; #it{E}_{track} (GeV); #it{E}_{cone} (GeV)", 100, 0, 750, 100, 0, 200);
	// TH2F histoTrackEnergija_konusa_pfo ("histoTrackEnergija_konusa_pfo", " ; Track energy (GeV); Cone energy (Gev)", 200, 0, 400, 80, 0, 80);

	// Brojacke promenljive za KRAJNJI elektron/pozitron
	Int_t N_ukupnih_dogadjaja = 0;
	Int_t N_KG_elektrona = 0;
	Int_t N_KG_kvarkova = 0;
	Int_t N_MC_elektronskih_signala = 0;
	Int_t N_rekonstruisanih_MC_elektronskih_signala = 0;
	Int_t N_truthlinkovanih_MC_elektronskih_signala = 0;

	// Brojacka promenljiva za POCETNI elektron 6 i pozitron 7
//	Int_t N_ukupnih_dogadjaja_6i7 = 0;
	Int_t N_KG_elektrona_6 = 0;
	Int_t N_KG_pozitrona_7 = 0;
//	Int_t N_MC_elektronskih_signala_6i7 = 0;
//	Int_t N_rekonstruisanih_MC_elektronskih_signala_6i7 = 0;
//	Int_t N_truthlinkovanih_MC_elektronskih_signala_6i7 = 0;

/*	Int_t nevent_uk = 0;
	Int_t nsignala_mc = 0;
	Int_t n_recomc_link_lep = 0;
	Int_t n_truthlep_kriva = 0;
	Int_t n_ukpfo = 0;              */

	Int_t N_za_dva = 0;
	int brDogadjajazaKrivu = 0;
	//int event_number;


	// Petlja koja iscitava .slcio fajlove
	for(UInt_t iFajl = nPrviFajl; iFajl <= nZadnjiFajl; iFajl++)
	{
		cout << "Otvara se " << Form("%s%i.slcio", fName.Data(), iFajl);

		try
		{
			lcReader -> open(Form("%s%i.slcio", fName.Data(), iFajl));
		}

		catch(lcio::IOException &ex)
		{
			cout << ". Ne moze se otvoriti.\n";
			continue;
		}

		cout << ". Ucitavanje.\n";

		// KG: ovde su Kragujevcani definisali broj dogadjaja u petlji
		Int_t N_dogadjaja = 0;

//		Int_t nevent = 0;

		EVENT::LCEvent* evt = 0;

		// WHILE petlja po dogadjajima
		while( (evt = lcReader -> readNextEvent()) != 0)
		{
			// KG: ovde su Kragujevcani definisali unutrasnje brojace
			N_dogadjaja++;
			N_ukupnih_dogadjaja++;
			TLorentzVector temp_e; // definisan privremeni četvorovektor u koji se zapisuju impuls i energija svakog pozitrona
			TLorentzVector temp_p;

			Int_t N_leptona_energija_ugao_po_dogadjaju = 0;
			vector <TLorentzVector> leptoni;
			vector<TLorentzVector> esps;


		    	bool elektronac = false;
			bool pozitronac = false;

		    	bool elektron2 = false;
			bool pozitron2 = false;

		    	bool elektron2cutted = false;
			bool pozitron2cutted = false;

		    	bool mcElektron = false;
			bool mcPozitron = false;

			double ekp = 0 ;
			double eke = 0;
/*			nevent++;
			nevent_uk ++;
			Int_t n_eventlok = 0;
			n_lep_evt = 0;
			n_track = 0;
			n_track_pt = 0;
			n_kriva = 0;
			n_svipfo = 0; */

			vector <EVENT::ReconstructedParticle*> niz_linklep, niz_pfos;

			// KG: ovde su Kragujevcani definisali potrebne nizove kvadri-vektora
			vector <TLorentzVector> niz_KG_elektrona, niz_KG_pozitrona, niz_KG_kvarkova;

//			vector <TLorentzVector> nizl, nizq;
			TLorentzVector estarac, pstarac, es, ps;

			std::vector<std::string> colNames = *evt -> getCollectionNames();

			EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticlesSkimmed");
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFOs");
//			IMPL::LCCollectionVec* pfo = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFANewPFOs");
//			IMPL::LCCollectionVec* iso_lep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected_new");
//			IMPL::LCCollectionVec* iso_lep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected_new_new_new");
//			IMPL::LCCollectionVec* iso_lep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected_2017_test");
//			IMPL::LCCollectionVec* iso_lep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected");
//			IMPL::LCCollectionVec* pfo_ptcut = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFOs_new");
//			IMPL::LCCollectionVec* bez_isolep = (IMPL::LCCollectionVec*) evt -> getCollection("WithoutIsoLep_Selected_new");
//			IMPL::LCCollectionVec* jets = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_2_2017_test");
//			IMPL::LCCollectionVec* jets = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_2");

			bool signal = false;
			vector <EVENT::ReconstructedParticle*> reclep_vector; // Definisan niz u koji se zapisuju rekonstrusani elektroni/pozitroni

			// Definisani broj elektrona, broj pozitrona i broj kvarkova
//			Int_t N_KG_kvarkova = 0;

/*			Int_t Broj_leptona = 0;
			Int_t Broj_kvarkova = 0; */

//			Float_t Brojac_za_dvicu = 0;
			Float_t e_number;


			for (Int_t i = 0; i < recParticles->getNumberOfElements() ; i++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles->getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
				Int_t particlePDG = recParticle->getType();

				if (particlePDG ==11 && e > 100) es = temp;
				if (particlePDG ==-11 && e > 100) ps = temp;

			}
		//	if (e_number == 2) event_number++;

			// FOR petlja u kojoj se traze signali kroz sve redove jedne kolekcije MCParticlesSkimmmed
			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				// ako je i = 2 to je pocetni elektron koji se sudara, ako je i = 3 to je pocetni pozitron, i trebaju nam njegove promenljive



			//	if (mcParticle -> getGeneratorStatus() != 1) continue;

				TLorentzVector temp; // privremeni četvorovektor u koji sakupljamo informacije o svakoj čestici

				const double *p = mcParticle -> getMomentum(); // impuls čestice
				double e = mcParticle -> getEnergy();          // energija čestice

				temp.SetPxPyPzE(p[0], p[1], p[2], e);          // zapisujemo vrednosti energije i impulsa u privremeni četvorovektor

				int PDG = mcParticle -> getPDG();

				const EVENT::MCParticleVec & parent = mcParticle -> getParents();
				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
				
			//	if (PDG == 25 && !(daughter[0]->getPDG() == 5 || daughter[0]->getPDG() == -5))  ;

				//cout <<"ispred sam \n";

				if (PDG == 25 && parent[0]->getPDG() == 11)
						/*&& parent[0]->getPDG()==11 && parent[0]->getParents().size()==0*/  {

						TLorentzVector estarac1 (TVector3 (parent[0] -> getMomentum()), parent[0] -> getEnergy());
						estarac = estarac1;
						TLorentzVector pstarac1 (TVector3 (parent[1] -> getMomentum()), parent[1] -> getEnergy());
						pstarac = pstarac1;
					//	cout <<"estarac1 : "<< estarac1.E();

				/*	TLorentzVector estarac1 (TVector3 (parent[0] -> getMomentum()), parent[0] -> getEnergy());
					TLorentzVector pstarac1 (TVector3 (parent[1] -> getMomentum()), parent[1] -> getEnergy());
					pstarac = pstarac1;*/



				}

				if (fabs(PDG) == 11 && temp.E() > 100 && mcParticle -> getGeneratorStatus() == 1)
				{
					if (temp.Theta()* 180/M_PI > 8  && temp.Theta()  * 180/M_PI < 172) N_leptona_energija_ugao_po_dogadjaju++;

				//	cout << " Redni broj dogadjaja: " << N_dogadjaja << endl;
				//	cout << " broj leptona u evt: " <<N_leptona_energija_ugao_po_dogadjaju<< endl;
					pt_KG_elektrona = temp.Pt();
					Energija_KG_elektrona = temp.E();
					Theta_KG_elektrona = temp.Theta() * 180/M_PI;
 					N_KG_elektrona++;

 					EVENT::LCRelation* link_to_rec_lepton1 = 0;
 					EVENT::MCParticle* pointer_to_mc_lepton1 = 0;
 					EVENT::ReconstructedParticle* pointer_to_rec_lepton1 = 0;
 					EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");

 					pointer_to_mc_lepton1 = mcParticle;

 					for (int j = 0; j < links -> getNumberOfElements(); j++)
 					{
	 					EVENT::LCRelation* linki = (EVENT::LCRelation*) links -> getElementAt(j);
 						EVENT::MCParticle* mcpj = (EVENT::MCParticle*) linki -> getTo();
					    	link_to_rec_lepton1 = 0;

					    	if (pointer_to_mc_lepton1 == mcpj)
 					    	{
 					    		link_to_rec_lepton1 = linki;
 					    	}


						if (link_to_rec_lepton1)
					    	{
 					    		EVENT::ReconstructedParticle* rpi1 = (EVENT::ReconstructedParticle*) link_to_rec_lepton1 -> getFrom();
 					    		pointer_to_rec_lepton1 = rpi1;
							// IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles -> getElementAt(i);

							if (rpi1 -> getType() == 11 && rpi1 -> getEnergy() > 100)
							{
								elektronac = true;

								//cout << " ******* Energija elektrona " << rpi1 -> getEnergy() << endl;

								// cout << "Dobijen linkovano-rekonstruisani elektron energije vece od 100 GeV-a" << rpi1 -> getType() << endl;

								TLorentzVector rekonstrE_temp; // definisan privremeni četvorovektor u koji se zapisuju impuls i energija svakog elektrona


								const double *p = rpi1 -> getMomentum(); // impuls svakog elektrona
								double e = rpi1 -> getEnergy();          // energija svakog elektrona
 							//	cout << " ******* Energija elektrona " << rpi1 -> getEnergy() << endl;

								rekonstrE_temp.SetPxPyPzE(p[0], p[1], p[2], e); // zapisuju se impuls i energija svakog elektrona u privremeni četvorovektor

								pt_KG_elektrona_rekonstr = rekonstrE_temp.Pt();
								Energija_KG_elektrona_rekonstr = rekonstrE_temp.E();
								Theta_KG_elektrona_rekonstr = rekonstrE_temp.Theta() * 180/M_PI;

								reclep_vector.push_back(rpi1); // ubacivanje u niz elektrona signala

								const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) rpi1 -> getTracks();

								if (trkvec.size() > 0)
								{
									d0_KG_elektrona_rekonstr = fabs (trkvec[0] -> getD0());
									z0_KG_elektrona_rekonstr = fabs (trkvec[0] -> getZ0());
									r0_KG_elektrona_rekonstr = sqrt (d0_KG_elektrona_rekonstr * d0_KG_elektrona_rekonstr + z0_KG_elektrona_rekonstr * z0_KG_elektrona_rekonstr);
								}

								// if (d0_KG_elektrona_rekonstr >= 0.04 || z0_KG_elektrona_rekonstr >= 0.1 || Energija_KG_elektrona_rekonstr <= 5.0) continue;
								// odbacuju se elektroni sa d0 vecim od 0.04 mm, sa z0 vecim od 0.1 mm i sa energijom manjom od 5 GeV

								std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) rpi1 -> getClusters();

								Float_t Ecal_leptona = 0;
								Float_t Hcal_leptona = 0;

								// energija elektrona deponovana u kalorimetrima
								for (std::vector<EVENT::Cluster*>::const_iterator iCluster = clusters.begin(); iCluster != clusters.end(); ++iCluster)
								{
									Ecal_leptona += (*iCluster) -> getSubdetectorEnergies()[0];
									Hcal_leptona += (*iCluster) -> getSubdetectorEnergies()[1];
								}

								odnos_Ecal_TOTALcal_KG_elektrona_rekonstr = Ecal_leptona / (Ecal_leptona + Hcal_leptona);

								// if (odnos_Ecal_TOTALcal_KG_elektrona_rekonstr < 0.9) continue;   // odbacuju se elektroni sa Rcal manjim od 0.9

								// Brojac_za_dvicu++;

								// Brojac_za_dvojku = Brojac_za_dvicu;

								Float_t Energija_konusa = 0;

								for (Int_t k = 0; k < recParticles -> getNumberOfElements(); k++)
								{
									EVENT::ReconstructedParticle* pfok = (EVENT::ReconstructedParticle*) recParticles -> getElementAt(k);

									TLorentzVector pfok_LV (TVector3 (pfok -> getMomentum()), pfok -> getEnergy());

									if (rekonstrE_temp == pfok_LV) continue;

									lep_pfoj_angle = rekonstrE_temp.Angle(pfok_LV.Vect()) * 180/M_PI;

									if (cos(lep_pfoj_angle) > 0.995 )
									{
										Energija_konusa += pfok_LV.Energy();
									}

								} // Kraj brojne petlje po PFO (ParticleFlowObjects) za konus

								Energija_konusa_KG_elektrona_rekonstr = Energija_konusa;
								eke = Energija_konusa;

							//	if (d0_KG_elektrona_rekonstr < 0.04 && z0_KG_elektrona_rekonstr < 0.1 && Energija_KG_elektrona_rekonstr > 100.0 && odnos_Ecal_TOTALcal_KG_elektrona_rekonstr >= 0.94){
								if (d0_KG_elektrona_rekonstr < 0.1 && z0_KG_elektrona_rekonstr < 1 && Energija_KG_elektrona_rekonstr > 100.0 && odnos_Ecal_TOTALcal_KG_elektrona_rekonstr >= 0.95){
								temp_e = rekonstrE_temp;
								leptoni.push_back(temp_e);

							         elektron2 = true;
								if((pow(Energija_konusa,2) < 900 * (0.01 * Energija_KG_elektrona_rekonstr + 0.01))) elektron2cutted = true;
											}

								//KG histoTrackEnergija_konusa.Fill(temp_e.E(), Energija_konusa);

								 // KG formula za izolacionu krivu, menjati var po potrebi



							} // end za IF getType() == 11 && rpi1 -> getEnergy() > 100

							if (rpi1 -> getType() == -11 && rpi1 -> getEnergy() > 100)
							{
								pozitronac = true;

								//cout << " ******* Energija pozitrona " << rpi1 -> getEnergy() << endl;

								// cout << "Dobijen linkovano-rekonstruisani pozitron energije vece od 100 GeV-a" << rpi1 -> getType() << endl;
								TLorentzVector rekonstr_temp; // definisan privremeni četvorovektor u koji se zapisuju impuls i energija svakog pozitrona
								 // definisan privremeni četvorovektor u koji se zapisuju impuls i energija svakog pozitrona
								const double *p = rpi1 -> getMomentum(); // impuls svakog pozitrona
								double e = rpi1 -> getEnergy();          // energija svakog pozitrona

								rekonstr_temp.SetPxPyPzE(p[0], p[1], p[2], e); // zapisuju se impuls i energija svakog pozitrona u privremeni četvorovektor

								pt_KG_pozitrona_rekonstr = rekonstr_temp.Pt();
								Energija_KG_pozitrona_rekonstr = rekonstr_temp.E();
								Theta_KG_pozitrona_rekonstr = rekonstr_temp.Theta() * 180/M_PI;

								reclep_vector.push_back(rpi1); // ubacivanje u niz leptona signala

								const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) rpi1 -> getTracks();

								if (trkvec.size() > 0)
								{
									d0_KG_pozitrona_rekonstr = fabs (trkvec[0] -> getD0());
									z0_KG_pozitrona_rekonstr = fabs (trkvec[0] -> getZ0());
									r0_KG_pozitrona_rekonstr = sqrt (d0_KG_pozitrona_rekonstr * d0_KG_pozitrona_rekonstr + z0_KG_pozitrona_rekonstr * z0_KG_pozitrona_rekonstr);
								}

								// if (d0_KG_pozitrona_rekonstr >= 0.04 || z0_KG_pozitrona_rekonstr >= 0.1 || Energija_KG_pozitrona_rekonstr <= 5.0) continue;
								// odbacuju se pozitroni sa d0 vecim od 0.04 mm, sa z0 vecim od 0.1 mm i sa energijom manjom od 5 GeV

								std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) rpi1 -> getClusters();

								Float_t Ecal_leptona = 0;
								Float_t Hcal_leptona = 0;

								// energija pozitrona deponovana u kalorimetrima
								for (std::vector<EVENT::Cluster*>::const_iterator iCluster = clusters.begin(); iCluster != clusters.end(); ++iCluster)
								{
									Ecal_leptona += (*iCluster) -> getSubdetectorEnergies()[0];
									Hcal_leptona += (*iCluster) -> getSubdetectorEnergies()[1];
								}

								odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr = Ecal_leptona / (Ecal_leptona + Hcal_leptona);


								// if (odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr < 0.9) continue;   // odbacuju se pozitroni sa Rcal manjim od 0.9

//								Brojac_za_dvicu++;

//								Brojac_za_dvojku = Brojac_za_dvicu;

								Float_t Energija_konusa = 0;

								for (Int_t k = 0; k < recParticles -> getNumberOfElements(); k++)
								{
									EVENT::ReconstructedParticle* pfok = (EVENT::ReconstructedParticle*) recParticles -> getElementAt(k);

									TLorentzVector pfok_LV (TVector3 (pfok -> getMomentum()), pfok -> getEnergy());

									if (rekonstr_temp == pfok_LV) continue;

									lep_pfoj_angle = rekonstr_temp.Angle(pfok_LV.Vect()) * 180/M_PI;

									if (cos(lep_pfoj_angle) > 0.995 )
									{
										Energija_konusa += pfok_LV.Energy();
									}

								} // Kraj brojne petlje po PFO (ParticleFlowObjects) za konus

								Energija_konusa_KG_pozitrona_rekonstr = Energija_konusa;
								ekp = Energija_konusa;

							//	if (d0_KG_pozitrona_rekonstr < 0.04 && z0_KG_pozitrona_rekonstr < 0.1 && Energija_KG_pozitrona_rekonstr > 100.0 && odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr >= 0.94)							{
								if (d0_KG_pozitrona_rekonstr < 0.1 && z0_KG_pozitrona_rekonstr < 1 && Energija_KG_pozitrona_rekonstr > 100.0 && odnos_Ecal_TOTALcal_KG_pozitrona_rekonstr >= 0.95)	{
								temp_p = rekonstr_temp;
								leptoni.push_back(temp_p);

								pozitron2 = true;
								if((pow(ekp,2) < 900 * (0.01 * Energija_KG_pozitrona_rekonstr + 0.01))) pozitron2cutted=true;
								}



							} // getType() == -11 && rpi1 -> getEnergy() > 100


							//	elektronskoRekonstrTree.Fill();


 						  } // end 332


 				  	} // end for link loop
						// cout << " ------- Energija_konusa_KG_elektrona_rekonstr: " <<elektron2  <<endl;
			/*		if (leptoni.size()==2{ //elektron2 == true && pozitron2 == true)
								brDogadjajazaKrivu++;
								histoTrackEnergija_konusa.Fill(temp_e.E(), eke);
								histoTrackEnergija_konusa.Fill(temp_p.E(), ekp);
								    }

					if (elektron2cutted == true && pozitron2cutted == true) {

								histoTrackEnergija_konusa_cutted.Fill(temp_e.E(), eke);
								histoTrackEnergija_konusa_cutted.Fill(temp_p.E(), ekp);
												}*/

					elektronskoTree.Fill();
				} // kraj FOR petlje za abs(PDF) = 11, E_traga > 100 GeV, izlazni kod = 1

				// IF uslov za Higsov bozon koji se raspada na ma koja dva kvarka

				// KG: ispituje se da li Higs ima elektron/pozitron za roditelja
				if (mcParticle -> getPDG() == 25 && fabs(mcParticle -> getParents()[0] -> getPDG()) == 11 )
				{
					// KG traže se potomci Higsa, koristi se za traženje signala u liniji 1617
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					const EVENT::MCParticleVec & parent = mcParticle -> getParents(); // Traži se roditelj Higsa
					const EVENT::MCParticleVec & cerke = parent[0] -> getDaughters(); // Traže se potomci roditelja Higsa

					for (int l = 0; l < (int) cerke.size(); l++) // Petlja po potomcima Higsa koji su naš signal
					{
						int PDG = cerke[l] -> getPDG();

						TLorentzVector cerkaTemp (TVector3 (cerke[l] -> getMomentum()), cerke[l] -> getEnergy());

						// Ako je potomak elektron 6
						if (PDG == 11) // promenljive za elektron 6
						{
							pt_KG_elektrona_6 = cerkaTemp.Pt();
							Energija_KG_elektrona_6 = cerkaTemp.E();
							Theta_KG_elektrona_6 = cerkaTemp.Theta() * 180/M_PI;
 							N_KG_elektrona_6++;
 							Ee6 = cerkaTemp.E();
 							esps.push_back(cerkaTemp);


						}

						// Ako je potomak pozitron 7
						if (PDG == -11) // promenljive za pozitron 7
						{
							pt_KG_pozitrona_7 = cerkaTemp.Pt();
							Energija_KG_pozitrona_7 = cerkaTemp.E();
							Theta_KG_pozitrona_7 = cerkaTemp.Theta() * 180/M_PI;
 							N_KG_pozitrona_7++;
 							Ep7 = cerkaTemp.E();
 							esps.push_back(cerkaTemp);


						}

					}


				} // Kraj IF uslova za Higsov bozon koji se raspada na ma koja dva kvarka
			} // Kraj FOR petlje u kojoj se trazi signal, mc loop

				// FOR petlja po ParticleFlowObjects (PFOs)
				for (Int_t j = 0; j < recParticles -> getNumberOfElements(); j++)
				{
					EVENT::ReconstructedParticle* pfoj = (EVENT::ReconstructedParticle*) recParticles -> getElementAt(j);
					bool provera = false;
					for (Int_t p = 0; p < (int)reclep_vector.size(); p++)
					{
					//	if (pfoj == reclep_vector[0] || pfoj == reclep_vector[1]) continue;
						if (pfoj == reclep_vector[p] ) provera =true;
					} // for za reclep_vector

					if (provera) continue;
						TLorentzVector pfoTemp (TVector3 (pfoj -> getMomentum()), pfoj -> getEnergy());

						Energija_KG_pfo = pfoj -> getEnergy();
						Theta_KG_pfo = pfoTemp.Theta() * 180/M_PI;

						const EVENT::TrackVec & trkvec = (const EVENT::TrackVec) pfoj -> getTracks();

						if (trkvec.size() > 0)
						{
							d0_KG_pfo = fabs (trkvec[0] -> getD0());
							z0_KG_pfo = fabs (trkvec[0] -> getZ0());
							r0_KG_pfo = sqrt (d0_KG_pfo * d0_KG_pfo + z0_KG_pfo * z0_KG_pfo);
						}

						std::vector<EVENT::Cluster*> clusters = (std::vector<EVENT::Cluster*>) pfoj -> getClusters();

						Float_t Ecal_leptona = 0;
						Float_t Hcal_leptona = 0;

						// energija PFOs deponovana u kalorimetrima
						for (std::vector<EVENT::Cluster*>::const_iterator iCluster = clusters.begin(); iCluster != clusters.end(); ++iCluster)
						{
							Ecal_leptona += (*iCluster) -> getSubdetectorEnergies()[0];
							Hcal_leptona += (*iCluster) -> getSubdetectorEnergies()[1];
						}

						odnos_Ecal_TOTALcal_KG_pfo = Ecal_leptona / (Ecal_leptona + Hcal_leptona);

						pfoTree.Fill();


				} // end of pfo loop


				mee_gen = (esps[0]+esps[1]).M();

			ptTree.Fill();
			elektronsko6i7Tree.Fill();

			if (pozitronac && elektronac) elektronskoRekonstrTree.Fill();

			if(N_leptona_energija_ugao_po_dogadjaju == 2) N_za_dva++;

			if (leptoni.size()==2/*elektron2 == true && pozitron2 == true*/) {
										brDogadjajazaKrivu++;
										histoTrackEnergija_konusa.Fill(temp_e.E(), eke);
										histoTrackEnergija_konusa.Fill(temp_p.E(), ekp);
										    }

							if (elektron2cutted == true && pozitron2cutted == true) {

										histoTrackEnergija_konusa_cutted.Fill(temp_e.E(), eke);
										histoTrackEnergija_konusa_cutted.Fill(temp_p.E(), ekp);
														}




				//test


		if (es.E() > 100 && ps.E()>100){
			//		cout <<"energija estarac: "<< estarac.E()<<endl;
				//	cout <<"energija pstarac: "<< pstarac.E()<<endl;
				//	cout << " masa para ep : "<<(es+ps).M()<<endl;
			mHg = (estarac+pstarac -es-ps).M();
			mZ1g = (estarac - es).M();
			mZ2g = (pstarac - ps).M();
			mee = (es+ps).M();
			Eee = (es+ps).E();
			thetaZ1 = (estarac - es).Theta();
			thetaZ2 = (pstarac-ps).Theta();
			pxZ1 = (estarac - es).Px();
			pyZ1 = (estarac - es).Py();
			pzZ1 = (estarac - es).Pz();
			pxZ2 = (pstarac - ps).Px();
			pyZ2 = (pstarac - ps).Py();
			pzZ2 = (pstarac - ps).Pz();
			thetaH = (estarac+pstarac -es-ps).Theta();
			pxH = (estarac+pstarac -es-ps).Px();
			pyH = (estarac+pstarac -es-ps).Py();
			pzH = (estarac+pstarac -es-ps).Pz();
			Epstarac = pstarac.E();
			Eestarac = estarac.E();

			Ees = es.E();
			Eps = ps.E();


			if (mee>200)testTree.Fill();
		}


		} // Kraj WHILE petlje po dogadjajima




/*		cout << "Zatvara se petlja po dogadjajima\n";
		cout << "Broj dogadjaja po fajlu iznosi " << nevent << endl; */

		lcReader -> close();

	} // Kraj petlje po fajlovima

	cout << "Ukupan broj elektrona signala: " << N_KG_elektrona << endl;
	cout << "Ukupan broj dogadjaja: " << N_ukupnih_dogadjaja << endl;
	cout << "Ukupan broj pocetnih elektrona 6: " << N_KG_elektrona_6 << endl;
	cout << "Ukupan broj pocetnih pozitrona 7: " << N_KG_pozitrona_7 << endl;
	cout << "Broj dogadjaja sa dva generisana leptona iznosi " << N_za_dva << endl;
	cout << "broj dogadjaja za krivu: " << brDogadjajazaKrivu << endl;
	//cout << "broj elektrona sa e > 100 GeV: "<< event_number<<endl;


/*	cout << "Ukupan broj kvarkova nastalih od Higsa iznosi " << N_KG_kvarkova << endl;
	cout << "Ukupan broj dogadjaja: " << nevent_uk << endl;
	cout << "Ukupan broj leptona signala na generisanom nivou: " << nsignala_mc*2 << endl;
	cout << "Ukupan broj truth-linkovanih leptona od rekonstruisanih leptona: " << n_recomc_link_lep << endl;
	// cout << "Ukupan broj dogadjaja sa 2 truth-leptona/evt posle cuts: " << n_truthlep_kriva << endl;
	cout << "Ukupan broj pfo cestica bez truth-linkovanih leptona: " << n_ukpfo << endl;
	cout << "Broj dogadjaja nakon odsecanja za Pt: " << n_lep_evt_pt_cutted << endl;
	cout << "Broj dogadjaja nakon odsecanja za Energija_leptonskog_traga: " << n_lep_evt_cutted << endl;
	cout << "Broj dogadjaja nakon odsecanja za d0: " << n_lep_evt_d0_cutted << endl;
	cout << "Broj dogadjaja nakon odsecanja za z0: " << n_lep_evt_z0_cutted << endl;
	cout << "Broj dogadjaja nakon odsecanja za r03d: " << n_lep_evt_r03d_cutted << endl;
	cout << "Broj dogadjaja nakon odsecanja za ratio: " << n_lep_evt_ratio_cutted << endl;
	cout << "Broj prezivelih leptona: " << brojac_prezivelih_leptona << endl;

	cout << "Efikasnost za d0 uslov: " << n_lep_evt_d0_cutted/n_recomc_link_lep * 100 << "%" << endl;
	cout << "Efikasnost za z0 uslov: " << n_lep_evt_z0_cutted/n_recomc_link_lep * 100 << "%" << endl;
	cout << "Efikasnost za r03d uslov: " << n_lep_evt_r03d_cutted/n_recomc_link_lep * 100 << "%" << endl;
	cout << "Efikasnost za Energiju leptonskog traga uslov: " << n_lep_evt_cutted/n_recomc_link_lep * 100 << "%"  << endl;
	cout << "Efikasnost za Ratio (Rcal) uslov: " << n_lep_evt_ratio_cutted/n_recomc_link_lep * 100 << "%" << endl;
	cout << "Efikasnost za sve uslove: " << brojac_prezivelih_leptona/n_recomc_link_lep * 100 << "%" << endl; */

	TString tfName(rfn);

	if(!tfName.EndsWith(".root")) tfName.Append(".root");

	TFile rootFile(tfName.Data(),"RECREATE");

	histo_niso_lep.Write();
	// histo_mcvect_size.Write();
	histoTrackEnergija_konusa.Write();
	histoTrackEnergija_konusa_cutted.Write();
	// histoTrackEnergija_konusa_pfo.Write();

	leptonTree.Write();
	ptTree.Write();
	pfoTree.Write();

	elektronskoTree.Write();
	elektronsko6i7Tree.Write();
	elektronskoRekonstrTree.Write();
	testTree.Write();

//	pozitronskoRekonstrTree.Write();

	rootFile.Write();
	rootFile.Close();

	return 0;

} // Kraj slcio2appTree funkcije

// -----------------------------------------------------------------------------------------------------------------------------

Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;

	UInt_t nPrviFajl = 1;
	if(argc > iarg)
	{
		nPrviFajl = atoi(argv[iarg]);
		iarg++;
	}

	UInt_t nZadnjiFajl = 50;    //  < ------------------------------------
	if(argc > iarg)
	{
		nZadnjiFajl  = atoi(argv[iarg]);
		iarg++;
	}

	TString fName = "test_";    // < -----------------------------------------
	if(argc > iarg)
	{
		 fName = argv[iarg]; 	// iarg++;
	}

	TString rfName = "KG_Signal_izolacija_leptona_uslovi.root";  // < ------------------------------------

	if(argc > iarg) rfName = argv[0];

	return slcio2appTree(nPrviFajl, nZadnjiFajl, fName.Data(), rfName.Data());
}
