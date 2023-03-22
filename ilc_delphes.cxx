/*
 * signal_eeH_mvaVar.cxx
 *
 *  Created on: Apr 16, 2021
 *      Author: tanja
 */


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
#include "math.h"


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
#endif

#include "stdlib.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <array>
#include <string>
using namespace std;

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	///TTree leptonTree ("leptonTree", "Generator particle tree");
	TTree Tree ("events", "Generator particle tree");

	////leptonTree.AddFriend("ptTree", "");

	Float_t m_onshell, m_offshell, m_ll, m_qq, m_H, m_h, Mass_higgsl3, Mass_higgsl3a, Mee, Mee_ZZ_fuzija, Mee_ZH;
	Float_t m_Z1, m_Z2, Mass_Z1, Mass_Z2, p_e1, p_e2, E_e1, E_e2, e_gen, p_gen, pt_q1, pt_q2, E_q1, E_q2, E_Z1, E_Z2, P_Z1, P_Z2;
	Float_t n_signal, n_signal_uk, n_evt_mva;
	Float_t E_vis, E_vis_E_H,pt_H, E_H;
	Float_t pt_e_sistema, pt_p_sistema, pt_e1, pt_e2;
	Float_t n_pfo, pt_miss, theta_H;
	Float_t logy_12, logy_23, logy_34;
	Float_t Btag1, Btag2;
	Float_t Phi;

	// Promenljive za POCETNI elektron 6
	//	Float_t brojac_preostalih_KG_elektrona_6 = 0;
	Float_t pt_KG_elektrona_6 = 0, p_KG_elektrona_6 = 0;
	Float_t Energija_KG_elektrona_6 = 0;
	Float_t Theta_KG_elektrona_6 = 0;
	//	Float_t odnos_Ecal_TOTALcal_KG_elektrona_6 = 0;
	Float_t d0_KG_elektrona_6 = 0, z0_KG_elektrona_6 = 0, r0_KG_elektrona_6 = 0;
	//	Float_t Energija_konusa_KG_elektrona_6 = 0;

	// Promenljive za POCETNI pozitron 7
	//	Float_t brojac_preostalih_KG_pozitrona_7 = 0;
	Float_t pt_KG_pozitrona_7 = 0, p_KG_pozitrona_7 = 0;
	Float_t Energija_KG_pozitrona_7 = 0;
	Float_t Theta_KG_pozitrona_7 = 0;
	//	Float_t odnos_Ecal_TOTALcal_KG_pozitrona_7 = 0;
	Float_t d0_KG_pozitrona_7 = 0, z0_KG_pozitrona_7 = 0, r0_KG_pozitrona_7 = 0;
	//	Float_t Energija_konusa_KG_pozitrona_7 = 0;


	//TanjaTree.Branch("m_H", &m_H, "m_H");
	//TanjaTree.Branch("m_qq", &m_qq, "m_qq"); // invarijatne mase b-dzeta???
	//TanjaTree.Branch("m_ll", &m_ll, "m_ll"); // invarijatne mase b-dzeta???
	Tree.Branch("n_signal", &n_signal, "n_signal");
	Tree.Branch("n_signal_uk", &n_signal_uk, "n_signal_uk");
	Tree.Branch("n_evt_mva", &n_evt_mva, "n_evt_mva");
	Tree.Branch("m_qq", &m_qq, "m_qq"); // invarijatne mase b-dzeta???
	Tree.Branch("m_H", &m_H, "m_H"); // masa Higzovog bozona
	Tree.Branch("m_Z1", &m_Z1, "m_Z1"); // masa 1. onshell Z-bozona
	Tree.Branch("m_Z2", &m_Z2, "m_Z2"); // masa 2. onshell Z-bozona
	Tree.Branch("Mass_Z1", &Mass_Z1, "Mass_Z1"); // masa 1. onshell Z-bozonka
	Tree.Branch("Mass_Z2", &Mass_Z2, "Mass_Z2"); // masa 2. onshell Z-bozonka
	Tree.Branch("E_Z1", &E_Z1, "E_Z1"); // masa 1. onshell Z-bozonka
	Tree.Branch("E_Z2", &E_Z2, "E_Z2"); // masa 2. onshell Z-bozonka
	Tree.Branch("P_Z1", &P_Z1, "P_Z1"); // masa 1. onshell Z-bozonka
	Tree.Branch("P_Z2", &P_Z2, "P_Z2"); // masa 2. onshell Z-bozonka
	Tree.Branch("p_e1", &p_e1, "p_e1"); // impuls elektrona
	Tree.Branch("p_e2", &p_e2, "p_e2"); // impuls pozitrona
	Tree.Branch("E_e1", &E_e1, "E_e1"); // energija elektrona
	Tree.Branch("E_e2", &E_e2, "E_e2"); // energija pozitrona
    Tree.Branch("pt_e1", &pt_e1, "pt_e1"); // energija elektrona
	Tree.Branch("pt_e2", &pt_e2, "pt_e2"); // energija pozitrona
	Tree.Branch("e_gen", &e_gen, "e_gen"); // energija pozitrona
	Tree.Branch("p_gen", &p_gen, "p_gen"); // energija pozitrona
	Tree.Branch("pt_e_sistema", &pt_e_sistema, "pt_e_sistema");
	Tree.Branch("pt_p_sistema", &pt_p_sistema, "pt_p_sistema");
	Tree.Branch("pt_q1", &pt_q1, "pt_q1"); // transverzalni impuls 1. kvarka
	Tree.Branch("pt_q2", &pt_q2, "pt_q2"); // transverzalni impuls 2. kvarka
	Tree.Branch("E_q1", &E_q1, "E_q1"); // energija 1. kvarka
	Tree.Branch("E_q2", &E_q2, "E_q2"); // energija2. kvarka
	Tree.Branch("theta_H", &theta_H, "theta_H"); // Polarni ugao Higzovog bozona
	Tree.Branch("E_vis", &E_vis, "E_vis"); // Vidljiva energija svih cestica dzumle
	Tree.Branch("E_vis_E_H", &E_vis_E_H, "E_vis_E_H"); // Razlika Vidljive energije svih cestica dzumle i Energije Higzovog bozona
	Tree.Branch("E_H", &E_H, "E_H"); // Energija Higzovog bozona
	Tree.Branch("pt_H", &pt_H, "pt_H"); // transverzalni impuls Higzovog bozona
	Tree.Branch("pt_miss", &pt_miss, "pt_miss"); // nedostajucí transverzalni impuls (otisao neutrinu?)
	//TanjaTree.Branch("logy_12", &logy_12, "logy_12");
	//TanjaTree.Branch("logy_23", &logy_23, "logy_23");
	Tree.Branch("Phi", &Phi, "Phi");
	Tree.Branch("M_H", &Mass_higgsl3, "M_H");
	Tree.Branch("M_H1", &Mass_higgsl3a, "M_H1");
	Tree.Branch("M_{ee}", &Mee, "M_{ee}");

	// Grane drveta za POCETNI generisani elektron 6
	Tree.Branch("pt_KG_elektrona_6", &pt_KG_elektrona_6, "pt_KG_elektrona_6");
	Tree.Branch("p_KG_elektrona_6", &p_KG_elektrona_6, "p_KG_elektrona_6");
	Tree.Branch("Theta_KG_elektrona_6", &Theta_KG_elektrona_6, "Theta_KG_elektrona_6");
//	TanjaTree.Branch("odnos_Ecal_TOTALcal_KG_elektrona_6", &odnos_Ecal_TOTALcal_KG_elektrona_6, "odnos_Ecal_TOTALcal_KG_elektrona_6");
	Tree.Branch("d0_KG_elektrona_6", &d0_KG_elektrona_6, "d0_KG_elektrona_6");
	Tree.Branch("z0_KG_elektrona_6", &z0_KG_elektrona_6, "z0_KG_elektrona_6");
	Tree.Branch("r0_KG_elektrona_6", &r0_KG_elektrona_6, "r0_KG_elektrona_6");
	Tree.Branch("Energija_KG_elektrona_6", &Energija_KG_elektrona_6, "Energija_KG_elektrona_6");
//	TanjaTree.Branch("Energija_konusa_KG_elektrona_6", &Energija_konusa_KG_elektrona_6, "Energija_konusa_KG_elektrona_6");

	// Grane drveta za POCETNI generisani pozitron 7
	Tree.Branch("pt_KG_pozitrona_7", &pt_KG_pozitrona_7, "pt_KG_pozitrona_7");
	Tree.Branch("p_KG_pozitrona_7", &p_KG_pozitrona_7, "p_KG_pozitrona_7");
	Tree.Branch("Theta_KG_pozitrona_7", &Theta_KG_pozitrona_7, "Theta_KG_pozitrona_7");
//	TanjaTree.Branch("odnos_Ecal_TOTALcal_KG_pozitrona_7", &odnos_Ecal_TOTALcal_KG_pozitrona_7, "odnos_Ecal_TOTALcal_KG_pozitrona_7");
	Tree.Branch("d0_KG_pozitrona_7", &d0_KG_pozitrona_7, "d0_KG_pozitrona_7");
	Tree.Branch("z0_KG_pozitrona_7", &z0_KG_pozitrona_7, "z0_KG_pozitrona_7");
	Tree.Branch("r0_KG_pozitrona_7", &r0_KG_pozitrona_7, "r0_KG_pozitrona_7");
	Tree.Branch("Energija_KG_pozitrona_7", &Energija_KG_pozitrona_7, "Energija_KG_pozitrona_7");
//	TanjaTree.Branch("Energija_konusa_KG_pozitrona_7", &Energija_konusa_KG_pozitrona_7, "Energija_konusa_KG_pozitrona_7");

	TH1F histo_y12 ("Histogram_y12", " ; -log y_{12}", 100, 0, 10);
	TH1F histo_y23 ("Histogram_y23", " ; -log y_{23}", 100, 0, 10);
	TH1F histo_y34 ("Histogram_y34", " ; -log y_{34}", 100, 0, 10);

	TH1F histo_nlep ("Histogram_nlep", " ; nlep", 100, 0, 2000);
	TH1F histo_npfo ("Histogram_npfo", " ; n_pfo", 100, 0, 10000);

	TH1F histo_Mee_ZZ_fusion ("h_Mee_ZZ_fusion", " ; M_{ee}", 100, 0, 1000);
	TH1F histo_Mee_ZH ("h_Mee_ZH", " ; M_{ee}", 100, 0, 1000);
    TH1F histo_Mee_all ("h_Mee_all", " ; M_{ee}", 100, 0, 1000);

	TH1F histo_Phi_ZZ_fusion ("h_Phi_ZZ_fusion", " ; #phi", 100, -3.7, 3.7);

	Int_t brojac = 0;
	Int_t Nov_dogadjaj = 0;
	Int_t nlep1 = 0;
	Int_t nlep2 = 0;
	Int_t brojac1 = 0;
	Int_t brojac2 = 0;
	Int_t brojac3 = 0;
	Int_t brojac4 = 0;
	Int_t brojac5 = 0;
	Int_t brojac6 = 0;
	Int_t brojac7 = 0;
	Int_t brojac8 = 0;
	Int_t brojac9 = 0;
	Int_t brojac10 = 0;
	Int_t brojac11 = 0;



	Int_t nevent_uk = 0;
	//Int_t n_signal_uk = 0;
	//Int_t n_evt_mva = 0;
	Int_t n_evt_2lep = 0;
	/////// bool right_evt = false;

	// Petlja koja iscitava .slcio fajlove
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{
		cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);
		try
		{
			lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));
		}
		catch(lcio::IOException &ex)
		{
			cout << ". Could not open.\n";
			continue;
		}
			cout << ". Reading.\n";

		Int_t nevent = 0;

	// WHILE petlja po dogadjajima

		EVENT::LCEvent* evt = 0;

		while( (evt = lcReader -> readNextEvent()) != 0 && nevent< 30000)
		{
			// KG: ovde su Kragujevcani definisali unutrasnje brojace

			/*
			Nov_dogadjaj++;
			Int_t Broj_H_bozona = 0;
			Int_t Broj_b_kvarkova = 0;
			Int_t nlep = 0;

			//vector<TLorentzVector> Z_niz;                              // prikupljaju se Z bozoni
			vector<TLorentzVector> bkvark_niz;                           // prikupljaju se b-kvarkovi

			//vector<TLorentzVector> Niz_1, Niz_2, Niz_H;                // deklarisu se nizovi 4vektora Niz_1, Niz_2 i Niz_H
			vector<TLorentzVector> Niz_1b, Niz_2b, Niz_H;                // deklarisu se nizovi 4vektora Niz_1b, Niz_2b i Niz_H

			//TLorentzVector Niz_1p, Niz_1n, Niz_2p, Niz_2n, Z1, Z2, OnShell4_Z, OffShell4_Z, Higs4_H;   // deklarisu se 4vektori Z1, Z2, H
			TLorentzVector Niz_el_in, Niz_el_in2, Niz_poz_in, Niz_el_out, Niz_poz_out, b1, b2, Higs4_H, Z1_boson, Z2_boson, Niz_elpoz_in;   // deklarisu se 4vektori b1, b2, H

			//EVENT::LCCollection* col = evt->getCollection("MCParticle");
			//std::vector<std::string> colNames = *evt -> getCollectionNames();
			//IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");
			// Int_t nz = 0;

			TLorentzVector higgs;   // 4vektor u koji se smestaju informacije o Higsovom bozonu i koji se koristi za boost
			TLorentzVector firstHiggs;
			TLorentzVector bkvark1;
			TLorentzVector bkvark2;
			TLorentzVector higgs1;

			TLorentzVector cerke_higsa;
*/



			nevent++;
			nevent_uk ++;

			Int_t n_isolep_evt = 0;
			Int_t n_signal = 0;

			std::vector<std::string> colNames = *evt -> getCollectionNames();
			vector <TLorentzVector> nizq_mc;

	//		IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticles");
	//		IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*)evt->getCollection("PFOs");

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticles");//MCParticlesSkimmed, MCParticles
			IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*) evt -> getCollection("PFOs");//PFOs, IsolatedElectrons
			//IMPL::LCCollectionVec* colJet = (IMPL::LCCollectionVec*)evt->getCollection("RefinedJets");
			//IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*)evt->getCollection("FJ_Jets_2");
			//IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*)evt->getCollection("Isolep_Selected_new");

			Int_t nqvarkova = 0;
			bool signal = false;
			TLorentzVector e_starac, p_starac; //dodao 4-vektore elektrona i pozitrona koji se sudaraju
			TLorentzVector e_cerka, p_cerka; //četvorovektor u koji sakupljamo informacije o svakoj čestici
			Float_t visibleEnergy = 0;

// JAJAJAJAJAAJ
			// Petlja preko koje prolazimo kroz sve cestice po svakom dogadjaju
		   	for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		   	{

				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				if (mcParticle -> getPDG() == 25 && mcParticle -> getParents()[0] -> getPDG() == 11)
				{
					const EVENT::MCParticleVec & parent = mcParticle -> getParents();
					for (Int_t h = 0; h < parent.size(); h++)
					{
						//TLorentzVector e_starac1, p_starac1; //četvorovektor u koji sakupljamo informacije o svakoj čestici

						Int_t PDG_starac = parent[h] -> getPDG();

						const double *pgen = parent[h] -> getMomentum(); // impuls čestice
						double egen = parent[h] -> getEnergy();	//energija čestice

						// cout << "Energija: " << parent[h] -> getEnergy() << endl;

						if (PDG_starac == 11 ) e_starac.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen); // KG definisana imena u liniji 157
						if (PDG_starac == -11 ) p_starac.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen);

					}//end za roditelje 163


					// KG traže se potomci Higsa, koristi se za traženje signala u liniji 1617
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					//const EVENT::MCParticleVec & parent = mcParticle -> getParents(); // Traži se roditelj Higsa
					const EVENT::MCParticleVec & cerke = parent[0] -> getDaughters(); // Traže se potomci roditelja Higsa

					for (int l = 0; l < (int) cerke.size(); l++) // Petlja po potomcima Higsa koji su naš signal
					{

						int PDG_cerka = cerke[l] -> getPDG();

						TLorentzVector cerkaTemp (TVector3 (cerke[l] -> getMomentum()), cerke[l] -> getEnergy());

						const double *pgen = cerke[l] -> getMomentum(); // impuls čestice
						double egen = cerke[l] -> getEnergy();	//energija čestice

						// cout << "Energija: " << cerke[h] -> getEnergy() << endl;

						if (PDG_cerka == 11 ) e_cerka.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen); // KG definisana imena u liniji 157
						if (PDG_cerka == -11 ) p_cerka.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen);


						// Ako je potomak elektron 6
						if (PDG_cerka == 11) // promenljive za elektron 6
						{
							pt_KG_elektrona_6 = cerkaTemp.Pt();
							Energija_KG_elektrona_6 = cerkaTemp.E();
							Theta_KG_elektrona_6 = cerkaTemp.Theta() * 180/M_PI;
 							//N_KG_elektrona_6++;


						}

						// Ako je potomak pozitron 7
						if (PDG_cerka == -11) // promenljive za pozitron 7
						{
							pt_KG_pozitrona_7 = cerkaTemp.Pt();
							Energija_KG_pozitrona_7 = cerkaTemp.E();
							Theta_KG_pozitrona_7 = cerkaTemp.Theta() * 180/M_PI;
 							//N_KG_pozitrona_7++;


						}

					}


				}///end if za čestice koje se sudaraju

				if (!(mcParticle -> getPDG() == 25 && mcParticle -> getParents()[0] -> getPDG() == 25 && mcParticle -> getParents()[0] -> getParents()[0]-> getPDG() == 11 && mcParticle -> getParents()[0] -> getParents()[0]-> getParents()[0]-> getParents().size() == 0))
				{
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t PDG_cestica = daughter[k] -> getPDG();
						TLorentzVector temp (TVector3(daughter[k] -> getMomentum()), daughter[k] -> getEnergy());
/*
						TLorentzVector bkvark1, bkvark2;
						Int_t PDG_bkvark = daughter[k] -> getPDG();
						const double *pgen_b = daughter[k] -> getMomentum(); // impuls čestice
						double egen_b = daughter[k] -> getEnergy();	//energija čestice
*/
						if (abs(PDG_cestica) ==5 )
						{
							nqvarkova++;
							nizq_mc.push_back(temp);
						}
/*
						if (PDG_bkvark == 5 ) bkvark1.SetPxPyPzE(pgen_b[0], pgen_b[1], pgen_b[2], egen_b); // KG definisana imena u liniji 157
						if (PDG_bkvark == -5 ) bkvark2.SetPxPyPzE(pgen_b[0], pgen_b[1], pgen_b[2], egen_b);
*/
					}
				}
			} // enf for za MC particles


/*
				//IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) col -> getElementAt(i);
				EVENT::MCParticle* mcParticle = (EVENT::MCParticle*) mcParticles -> getElementAt(i);
				const EVENT::MCParticleVec & parent = mcParticle -> getParents();
				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

				//TLorentzVector temp;                             // 4vektor u koji se smestaju informacije o svakoj cestici
				//TVector3 temp1;
				vector<Int_t> PDG_cestica;			 // vektor u koji se smestaju PDGovi svake cestice

				//const double *p = mcParticle -> getMomentum();   // impuls cestice
				//double e = mcParticle -> getEnergy();            // energija cestice
				//int status = mcParticle -> getGeneratorStatus();
				//int pdg = mcParticle -> getPDG();
				//temp.SetPxPyPzE(p[0], p[1], p[2], e);  	         // zapisuju se vrednosti energije i impulsa u 4vektor temp
				//temp1.SetPxPyPz(p[0], p[1], p[2]);  	         // zapisuju se vrednosti energije i impulsa u 4vektor temp


				if (mcParticle -> getPDG() == 11)  //upadni elektron
				{

					nlep1++;
					//Niz_1b.push_back(temp);
					TLorentzVector temp_roditelj1;
					const double *p = mcParticle -> getMomentum();   // impuls cestice
					double e = mcParticle -> getEnergy();            // energija cestice
					int status = mcParticle -> getGeneratorStatus();
					int pdg = mcParticle -> getPDG();
					temp_roditelj1.SetPxPyPzE(p[0], p[1], p[2], e);  	         // zapisuju se vrednosti energije i impulsa u 4vektor temp
					Niz_el_in = temp_roditelj1;

					cout << "PDG roditelja 11 je: " << pdg << endl;
					cout << "Status roditelja 11 je: " << status << endl;

					for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t pdgdaughters = daughter[k] -> getPDG();
						Int_t status = daughter[k] -> getGeneratorStatus();
						cout << "Status cerke je: " << status << endl;
						cout << "PDG cerke je CHECK: " << pdgdaughters << endl;
						brojac2++;
						if (pdgdaughters == 11)
						{
							brojac3++;
							TLorentzVector temp_cerke1;                     // 4vektor u koji se skupljaju informacije o svakoj cestici
							const double *p = daughter[k] -> getMomentum();   // impuls cestice
							double e = daughter[k] -> getEnergy();	          // energija cestice
							temp_cerke1.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor temp_cerke
							Niz_el_out = temp_cerke1;
							//cerke_niz.push_back(temp_cerke);

						//cout << "PDG cerke 11 je: " << pdgdaughters << endl;
						//cout << "Status cerke 11 je: " << status << endl;

						}
						if (pdgdaughters == -11)
						{
							brojac9++;

						}

					}

				}

				if (mcParticle -> getPDG() == -11) //upadni pozitron
				{
					nlep2++;
					//Niz_1b.push_back(temp);
					TLorentzVector temp_roditelj2;
					const double *p = mcParticle -> getMomentum();   // impuls cestice
					double e = mcParticle -> getEnergy();            // energija cestice
					int status = mcParticle -> getGeneratorStatus();
					int pdg = mcParticle -> getPDG();
					temp_roditelj2.SetPxPyPzE(p[0], p[1], p[2], e);  	         // zapisuju se vrednosti energije i impulsa u 4vektor temp
					Niz_poz_in = temp_roditelj2;

					cout << "PDG roditelja -11 je: " << pdg << endl;
					cout << "Status roditelja -11 je: " << status << endl;

					for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t pdgdaughters = daughter[k] -> getPDG();
						Int_t status = daughter[k] -> getGeneratorStatus();
						if (pdgdaughters == -11)
						{
							TLorentzVector temp_cerke2;                     // 4vektor u koji se skupljaju informacije o svakoj cestici
							const double *p = daughter[k] -> getMomentum();   // impuls cestice
							double e = daughter[k] -> getEnergy();	          // energija cestice
							temp_cerke2.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor temp_cerke
							Niz_poz_out = temp_cerke2;
							//cerke_niz.push_back(temp_cerke);

						cout << "PDG cerke -11 je: " << pdgdaughters << endl;
						cout << "Status cerke -11 je: " << status << endl;

						}

					}

				}


				if (mcParticle -> getPDG() == 11 || mcParticle -> getPDG() == -11)
				{

				for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t pdgdaughters = daughter[k] -> getPDG();
						Int_t status = daughter[k] -> getGeneratorStatus();
						if (pdgdaughters == 25)
						{
							TLorentzVector temp_cerke2;                     // 4vektor u koji se skupljaju informacije o svakoj cestici
							const double *p = daughter[k] -> getMomentum();   // impuls cestice
							double e = daughter[k] -> getEnergy();	          // energija cestice
							temp_cerke2.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor temp_cerke
							higgs = temp_cerke2;
							//cerke_niz.push_back(temp_cerke);

						cout << "PDG cerke 25 je: " << pdgdaughters << endl;
						cout << "Status cerke 25 je: " << status << endl;

						}
					}
				}

				if (mcParticle -> getPDG() == 25)
				{
					for (Int_t m = 0; m < daughter.size(); m++)

								{
									Int_t pdgdaughters1 = daughter[m] -> getPDG();
									Int_t status1 = daughter[m] -> getGeneratorStatus();
									if (pdgdaughters1 == 5 || pdgdaughters1 == -5)
									{
										TLorentzVector temp_cerke3;                     // 4vektor u koji se skupljaju informacije o svakoj cestici
										const double *p = daughter[m] -> getMomentum();   // impuls cestice
										double e = daughter[m] -> getEnergy();	          // energija cestice
										temp_cerke3.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor temp_cerke
										cerke_higsa = temp_cerke3;
										//cerke_niz.push_back(temp_cerke);

										cout << "PDG cerke 5 je: " << pdgdaughters1 << endl;
										cout << "Status cerke 5 je: " << status1 << endl;

										nqvarkova++;
										//nizq_mc.push_back(temp);
										cout << "Ukupan br. kvarkova u kolekciji"<< nqvarkova << endl;

									}

										//cout << "Ukupan br. kvarkova u kolekciji"<< nqvarkova << endl;
								}

				}

			}

*/

			if (nqvarkova == 2 )
		  	{
				signal = true;
		  	}

			if (signal)
		  	{

/*
			TLorentzVector Z1;
			TLorentzVector Z2;

	//		   Double_t Mass_OnShell_Z = 0;
	//		   Double_t Mass_OffShell_Z = 0;
			m_Z1 = 0;
			m_Z2 = 0;

			Z1 = Niz_el_in + Niz_el_out;
			Z2 = Niz_poz_in + Niz_poz_out;


			m_Z1 = Z1.M();
			m_Z2 = Z2.M();

			TLorentzVector Niz_el_in4, Niz_poz_in4, Niz_el_out4, Niz_poz_out4, Higs4_H4;
			//TLorentzVector Z1;
			//TLorentzVector Z2;
			//Niz_el_in.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());
			//Niz_poz_in.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());
			//Niz_el_out.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());
			//Niz_poz_out.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());

			TLorentzVector onshellZ;
			//onshellZ.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());
			TVector3 BoostToHiggs = -(higgs.BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
			//TVector3 BoostToHiggs = -((bkvark1 + bkvark2).BoostVector());
			TVector3 BoostToLab = -((Niz_el_out+higgs+Niz_poz_out).BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
			TVector3 BoostToOnShellZ = -(Z1.BoostVector());     // prelazi se u koordinatni sistem OnShell Z bozona
			TVector3 BoostToOffShellZ = -(Z2.BoostVector());   // prelazi se u koordinatni sistem OffShell Z bozona
			TVector3 BoostToFinalee = -((Niz_el_out+Niz_poz_out).BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
			//TVector3 BoostToFinalee = ((Niz_el_in+Niz_poz_in).BoostVector());

			//TVector3 BoostToLab((Niz_el_in+Niz_poz_in).Px(), (Niz_el_in+Niz_poz_in).Py(), (Niz_el_in+Niz_poz_in).Pz());
			//TVector3 BoostToLab (0, 0, (Niz_el_in+Niz_el_out).Pz()/(Niz_el_in+Niz_el_out).E());
			//TVector3 boost (0, 0, p->P()/p->E());

			TVector3 Niz_el_in3;   // lepton 1 od Z
			TVector3 Niz_poz_in3;   // lepton 2 od Z

			TVector3 Niz_el_out3;   // lepton 1 od Z*
			TVector3 Niz_poz_out3;   // lepton 2 od Z*

			TVector3 OnShell3_l1;   // lepton 1 od Z
			TVector3 OnShell3_l2;   // lepton 2 od Z

			TVector3 OffShell3_l1;   // lepton 1 od Z*
			TVector3 OffShell3_l2;   // lepton 2 od Z*

			TVector3 OnShell3_Z;   // ovo je Z
			TVector3 OffShell3_Z;  // ovo je Z*
			TVector3 OnShell3_Z1;

			TVector3 higgsl3;

			TVector3 nz;

			//lokalac

			TLorentzVector lokalac_Niz_el_in4, lokalac_Niz_poz_in4, lokalac_Niz_el_out4, lokalac_Niz_poz_out4, lokalacZ1, lokalacZ2, lokalacZ1_1, lokalac_higgs;
		//	TLorentzVector lokalacZ1, lokalacZ2;
			lokalac_Niz_el_in4 = Niz_el_in;
			lokalac_Niz_el_out4 = Niz_el_out;
		    	lokalac_Niz_poz_in4 = Niz_poz_in;
			lokalac_Niz_poz_out4 = Niz_poz_out;
			lokalacZ1 = Niz_el_in-Niz_el_out; //4-vektor Z1-bozona Ovo je OK!
			lokalacZ2 = Niz_poz_in-Niz_poz_out; //4-vektor Z2-bozona OK je OK!
			//lokalacZ1_1 = lokalacZ2-higgs;
			//higgs = 2 * lokalacZ1;
			//lokalac_higgs = bkvark1 + bkvark2;

/*
			//*****BOOST TO LAB*****
			lokalac_Niz_el_in4.Boost(BoostToLab);
			lokalac_Niz_poz_in4.Boost(BoostToLab);
			lokalac_Niz_el_out4.Boost(BoostToLab);
			lokalac_Niz_poz_out4.Boost(BoostToLab);
			lokalacZ1.Boost(BoostToLab);
			lokalacZ2.Boost(BoostToLab);
			lokalacZ1_1.Boost(BoostToLab);
			higgs.Boost(BoostToLab);
			////Niz_el_in.Boost(BoostToHiggs);
			////Niz_el_out.Boost(BoostToHiggs);
*/
/*
			//*****BOOST TO Z*****
			lokalac_Niz_el_in4.Boost(BoostToOnShellZ);
			lokalac_Niz_poz_in4.Boost(BoostToOnShellZ);
			lokalac_Niz_el_out4.Boost(BoostToOnShellZ);
			lokalac_Niz_poz_out4.Boost(BoostToOnShellZ);
			lokalacZ1.Boost(BoostToOnShellZ);
			lokalacZ2.Boost(BoostToOnShellZ);
			lokalacZ1_1.Boost(BoostToOnShellZ);
			lokalac_higgs.Boost(BoostToOnShellZ);
			////Niz_el_in.Boost(BoostToHiggs);
			////Niz_el_out.Boost(BoostToHiggs);
*/
/*
			//*****BOOST TO HIGGS****
			lokalac_Niz_el_in4.Boost(BoostToHiggs);
			lokalac_Niz_poz_in4.Boost(BoostToHiggs);
			lokalac_Niz_el_out4.Boost(BoostToHiggs);
			lokalac_Niz_poz_out4.Boost(BoostToHiggs);
			lokalacZ1.Boost(BoostToHiggs);
			lokalacZ2.Boost(BoostToHiggs);
			lokalacZ1_1.Boost(BoostToHiggs);
			higgs.Boost(BoostToHiggs);
*/
/*
			//*****BOOST TO FINAL e-e+****
			lokalac_Niz_el_in4.Boost(BoostToFinalee);
			lokalac_Niz_poz_in4.Boost(BoostToFinalee);
			lokalac_Niz_el_out4.Boost(BoostToFinalee);
			lokalac_Niz_poz_out4.Boost(BoostToFinalee);
			lokalacZ1.Boost(BoostToFinalee);
			lokalacZ2.Boost(BoostToFinalee);
			lokalacZ1_1.Boost(BoostToFinalee);
			higgs.Boost(BoostToFinalee);
*/
/*
			Niz_el_in3.SetXYZ(lokalac_Niz_el_in4.X(), lokalac_Niz_el_in4.Y(), lokalac_Niz_el_in4.Z());
			Niz_el_out3.SetXYZ(lokalac_Niz_el_out4.X(), lokalac_Niz_el_out4.Y(), lokalac_Niz_el_out4.Z());
		    	Niz_poz_in3.SetXYZ(lokalac_Niz_poz_in4.X(), lokalac_Niz_poz_in4.Y(), lokalac_Niz_poz_in4.Z());
		    	Niz_poz_out3.SetXYZ(lokalac_Niz_poz_out4.X(), lokalac_Niz_poz_out4.Y(), lokalac_Niz_poz_out4.Z());

			//OnShell3_l1.SetXYZ(lokalac_OnShell4_l1.X(), lokalac_OnShell4_l1.Y(), lokalac_OnShell4_l1.Z());
			//OnShell3_l2.SetXYZ(lokalac_OnShell4_l2.X(), lokalac_OnShell4_l2.Y(), lokalac_OnShell4_l2.Z());

			OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z()); // impuls Z1
			OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z()); // impuls Z2
			///OnShell3_Z1.SetXYZ(lokalacZ1_1.X(),lokalacZ1_1.Y(),lokalacZ1_1.Z()); // impuls OnShell Z

  	       	        higgsl3.SetXYZ(higgs.X(), higgs.Y(), higgs.Z());

			nz.SetXYZ(0,0,1);   // z-osa

		    	//angle_onsh = (OnShell3_Z.Angle(higgsl3))* 180/M_PI;
			//angle_offsh = (OffShell3_Z.Angle(higgsl3))* 180/M_PI;
			//angle_q1q2 = (OnShell3_Z.Angle(OffShell3_Z)) * 180/M_PI;

			Double_t thetaH = higgsl3.Theta();
			thetaH1 = TMath::Cos(thetaH);

			//------ NASA DEFINICIJA UGLA PHI - u ovom slucaju su normale na ravni usmerene isto ------//

			TVector3 brojilac1 = Niz_el_in3.Cross(Niz_el_out3);   // brojilac vektora n1
			Double_t imenilac1 = brojilac1.Mag();   // sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));   // imenilac vektora n1

			TVector3 n1 = brojilac1 * pow(imenilac1,-1);           // vektor n1
			//sin_n1 = sin(Niz_el_in3.Angle(Niz_el_out3));
		        //sin_n11= brojilac1.Mag() * (pow (Niz_el_in3.Mag(), -1)) * pow(Niz_el_out3.Mag(), -1);

			TVector3 brojilac2 = Niz_poz_in3.Cross(Niz_poz_out3); // brojilac vektora n2
			Double_t imenilac2 = brojilac2.Mag(); // sqrt(pow(brojilac2.X(),2) + pow(brojilac2.Y(),2) + pow(brojilac2.Z(),2));   // imenilac vektora n1

		        TVector3 n2 = brojilac2 * pow(imenilac2,-1);           // vektor n2
		        //sin_n2 = sin(Niz_poz_in3.Angle(Niz_poz_out3));
		        //sin_n21 = brojilac2.Mag() * (pow(Niz_poz_in3.Mag(), -1)) * pow(Niz_poz_out3.Mag(), -1);
		        //cos_n12 = n1.Dot(n2) * pow(n1.Mag(), -1) * pow(n2.Mag(), -1);

			fi1 = acos(n1.Dot(n2));      // * 180 / M_PI prvi ugao
			fi2 = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));

						//------ OGAWINA DEFINICIJA UGLA NA OSNOVU FORMULE I GRAFICKE INTERPRETACIJE SLIKE - u ovom slucaju su normale na ravni usmerene razlicito ------//
			TVector3 brojilac3 = OnShell3_Z.Cross(Niz_el_in3);   // brojilac vektora n1
			Double_t imenilac3 = brojilac3.Mag();   // sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));   // imenilac vektora n1

			TVector3 n3 = brojilac3 * pow(imenilac3,-1);          // vektor n1
			//sin_n1 = sin(OnShell3_Z.Angle(OffShell3_Z));
		        //sin_n11= brojilac3.Mag() * (pow (OnShell3_Z.Mag(), -1)) * pow(OffShell3_Z.Mag(), -1);

			TVector3 brojilac4 = OnShell3_Z.Cross(Niz_poz_out3); // brojilac vektora n2
			Double_t imenilac4 = brojilac4.Mag(); // sqrt(pow(brojilac4.X(),2) + pow(brojilac4.Y(),2) + pow(brojilac4.Z(),2));   // imenilac vektora n1

		        TVector3 n4 = brojilac4 * pow(imenilac4,-1);          // vektor n2
		        //sin_n2 = sin(Niz_el_out3.Angle(Niz_poz_out3));
		        //sin_n21 = brojilac4.Mag() * (pow(Niz_el_out3.Mag(), -1)) * pow(Niz_poz_out3.Mag(), -1);
		        //cos_n12 = n1.Dot(n2) * pow(n1.Mag(), -1) * pow(n2.Mag(), -1);

			fi3 = acos(-n3.Dot(n4));
			fi4 = OnShell3_Z.Dot(n3.Cross(n4)) * fabs (pow(OnShell3_Z.Dot(n3.Cross(n4)),-1)) * acos(-n3.Dot(n4));

			//------ SA UZIMANJEM Z2 BOZONA, A NE Z1 BOZONA U DRUGOJ RAVNI - u ovom slucaju su normale na ravni usmerene isto ------//
			TVector3 brojilac5 = OnShell3_Z.Cross(Niz_el_in3);   // brojilac vektora n1
			Double_t imenilac5 = brojilac5.Mag();
			TVector3 n5 = brojilac5 * pow(imenilac5,-1);
			TVector3 brojilac6 = OffShell3_Z.Cross(Niz_poz_out3); // brojilac vektora n2
			Double_t imenilac6 = brojilac6.Mag();
			TVector3 n6 = brojilac6 * pow(imenilac6,-1);

			fi5 = acos(n5.Dot(n6));
			fi6 = OnShell3_Z.Dot(n5.Cross(n6)) * fabs (pow(OnShell3_Z.Dot(n5.Cross(n6)),-1)) * acos(n5.Dot(n6)); 	*/

			// FOR petlja u kojoj se traze signali kroz sve redove jedne kolekcije MCParticlesSkimmmed
/*
			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		   	{

				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

				if (mcParticle -> getPDG() == 25 && mcParticle -> getParents()[0] -> getPDG() == 25 && mcParticle -> getParents()[0] -> getParents()[0]-> getPDG() == 11 && mcParticle -> getParents()[0] -> getParents()[0]-> getParents()[0]-> getParents().size() == 0)
				{
					const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();
					for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t PDG_cestica = abs (daughter[k] -> getPDG());
						TLorentzVector temp (TVector3(daughter[k]->getMomentum()), daughter[k]->getEnergy());

						if (  PDG_cestica == 5 )
						{
							nqvarkova++;
							nizq_mc.push_back(temp);
						}

					}

				}

			}

			if (nqvarkova == 2 )
		   	{
				signal = true;
		  	 }

		 	 if (signal)
		 	 {
*/
				n_signal++;
				n_signal_uk++;
		//		if (n_signal_uk > 6192 && n_signal_uk <= 18576)
		//		{
				n_evt_mva++;
				Int_t nlep = 0;
				Int_t nqu = 0;
				vector <TLorentzVector> v1l, v2j;
				TLorentzVector es, ps;

				E_vis = 0;
				n_pfo = 0;
				pt_miss = 0;
				TLorentzVector sum_vec;

				for (Int_t k = 0; k < pfos->getNumberOfElements(); k++)
			 	{
					IMPL::ReconstructedParticleImpl* pfo = (IMPL::ReconstructedParticleImpl*) pfos->getElementAt(k);
					TLorentzVector pfoLV (TVector3(pfo->getMomentum()), pfo->getEnergy());
					E_vis += pfoLV.Energy();
					Double_t E_track = pfo->getEnergy();
					n_pfo++;

					sum_vec += pfoLV;
			 	}

					//histo_npfo.Fill(n_pfo);
		//			cout << "Ukupan br. pfos u kolekciji" << n_pfo << endl;

					pt_miss = sum_vec.Pt();

					vector<TLorentzVector> pfo_jets;
					vector<TLorentzVector> iso_lep;

				Float_t m_y12 = 0;
				Float_t m_y23 = 0;

		/*		for (Int_t i = 0; i < jets2 -> getNumberOfElements(); i++)
			 	{
				 IMPL::ReconstructedParticleImpl* recJet = (IMPL::ReconstructedParticleImpl*) jets2 -> getElementAt(i);
				 TLorentzVector tempjet;
				 nqu++;
				 const double *p = recJet -> getMomentum();
				 Double_t e = recJet -> getEnergy();
				 tempjet.SetPxPyPzE(p[0], p[1], p[2], e);
				 Int_t pidjet = recJet -> getType();
				 v2j.push_back(tempjet);

	     			 float yMinus = static_cast<double>(jets2 -> parameters().getFloatVal("y_{n-1,n}"));
				 float yPlus  = static_cast<double>(jets2 -> parameters().getFloatVal("y_{n,n+1}"));
				 // float dMinus = static_cast<double>(jets2 -> parameters().getFloatVal("d_{n-1,n}"));
	     			 // float dPlus  = static_cast<double>(jets2 -> parameters().getFloatVal("d_{n,n+1}"));

				 m_y12 = -log10(yMinus);
				 m_y23 = -log10(yPlus);

				//mqq = (v2j[0] + v2j[1]).M();

			//	 cout << "pidjet" << pidjet << endl;
				 }

			//	cout << "Ukupan br. jets u kolekciji"<< nqu << endl;

				histo_y12.Fill(m_y12);
				histo_y23.Fill(m_y23);
*/
	/*			for (Int_t j = 0; j < isolep->getNumberOfElements(); j++)
				{
				 	IMPL::ReconstructedParticleImpl* ilf = (IMPL::ReconstructedParticleImpl*) isolep->getElementAt(j);
				 	if ( abs (ilf -> getType()) != 11) continue;

					TLorentzVector templep;
					Int_t pid = ilf->getType();
					nlep++;

				 	const double *p = ilf->getMomentum();
					Double_t e = ilf->getEnergy();
					templep.SetPxPyPzE(p[0], p[1], p[2], e);
					v1l.push_back(templep);
					if (pid == 11) es = templep;
			 		if (pid == -11)ps = templep;

			//		cout << "pid" << pid << endl;

				//cout << "Ukupan br. leptona u kolekciji"<< nlep << endl;

				}*/

		//		cout << "Ukupan br. leptona1 u kolekciji"<< nlep << endl;

				//histo_nlep.Fill(nlep);

			if (nlep == 2)
			{

			 	m_ll = 0; m_qq = 0;

				m_ll = (v1l[0] + v1l[1]).M();
			 	m_qq = (v2j[0] + v2j[1]).M();

				p_e1 = es.P(); // KG dodao e1 je elektron, e2 pozitron
				E_e1 = es.E();// KG dodao e1 je elektron, e2 pozitron
				pt_e1 = es.Pt();

				p_e2 = ps.P();// KG dodao e1 je elektron, e2 pozitron
				E_e2 = ps.E();// KG dodao e1 je elektron, e2 pozitron
				pt_e2 = ps.Pt();

				pt_q1 = v2j[0].Pt();
				E_q1 = v2j[0].E();

				pt_q2 = v2j[1].Pt();
				E_q2 = v2j[1].E();

				pt_e_sistema = (e_starac + es).Pt();
				pt_p_sistema = (p_starac + ps).Pt();

				m_Z1 = (e_starac + es).M();
				m_Z2 = (p_starac + ps).M();

				//E_Z1 = (e_starac + es).E();
				//E_Z2 = (p_starac + ps).E();

				e_gen = e_starac.E();
				p_gen = p_starac.E();

				TLorentzVector Z1;
				TLorentzVector Z2;

				Z1 =  e_starac - e_cerka;
				Z2 =  p_starac - p_cerka;

/*
				 if (mll <= mqq)
				 {
					 m_onshell = mqq;
					 m_offshell = mll;

				 }
				 else
				 {
					 m_onshell = mll;
					 m_offshell = mqq;
				 }
*/
				 TLorentzVector on, off, higgs;

				 higgs = v2j[0] + v2j[1];
				 m_h = higgs.M();
				//if (m_h < 140 && m_h > 110)
				//{
					m_H = higgs.M();
				//}
				 E_H = higgs.E();
				 E_vis_E_H = abs (E_vis - E_H);
				 theta_H = higgs.Theta() * 180/M_PI;
				 pt_H = higgs.Pt(); // KG ovo dodati u drvo
				 //logy_12 = m_y12;
				 //logy_23 = m_y23;

				//DEFINICIJA UGLA PHI -->> Tanja

				TVector3 BoostToHiggs = -(higgs.BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
				//TVector3 BoostToHiggs = -((bkvark1 + bkvark2).BoostVector());
				TVector3 BoostToLab = -((es+higgs+ps).BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
				//TVector3 BoostToZ1Z2 = -((Z1+Z2).BoostVector());
				TVector3 BoostToFinalee = -((es+ps).BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona

				TVector3 Niz_el_in3;   // lepton 1 od Z
				TVector3 Niz_poz_in3;   // lepton 2 od Z

				TVector3 Niz_el_out3;   // lepton 1 od Z*
				TVector3 Niz_poz_out3;   // lepton 2 od Z*

				TVector3 OnShell3_l1;   // lepton 1 od Z
				TVector3 OnShell3_l2;   // lepton 2 od Z

				TVector3 OffShell3_l1;   // lepton 1 od Z*
				TVector3 OffShell3_l2;   // lepton 2 od Z*

				TVector3 OnShell3_Z;   // ovo je Z
				TVector3 OffShell3_Z;  // ovo je Z*
				TVector3 OnShell3_Z1;

				TVector3 higgsl3;

				//lokalac

				TLorentzVector lokalac_Niz_el_in4, lokalac_Niz_poz_in4, lokalac_Niz_el_out4, lokalac_Niz_poz_out4, lokalacZ1, lokalacZ2, lokalacZ1_1, lokalac_higgs;
				//TLorentzVector lokalacZ1, lokalacZ2;
				lokalac_Niz_el_in4 = e_starac;
				lokalac_Niz_el_out4 = es;
		    		lokalac_Niz_poz_in4 = p_starac;
				lokalac_Niz_poz_out4 = ps;
				lokalacZ1 = Z1; //Niz_el_in-Niz_el_out; //4-vektor Z1-bozona Ovo je OK!
				lokalacZ2 = Z2; //Niz_poz_in-Niz_poz_out; //4-vektor Z2-bozona OK je OK!
				//lokalacZ1_1 = lokalacZ2-higgs;
				//higgs = 2 * lokalacZ1;
				////lokalac_higgs = bkvark1 + bkvark2;

				//*****BOOST TO HIGGS****
				lokalac_Niz_el_in4.Boost(BoostToHiggs);
				lokalac_Niz_poz_in4.Boost(BoostToHiggs);
				lokalac_Niz_el_out4.Boost(BoostToHiggs);
				lokalac_Niz_poz_out4.Boost(BoostToHiggs);
				lokalacZ1.Boost(BoostToHiggs);
				lokalacZ2.Boost(BoostToHiggs);
				lokalacZ1_1.Boost(BoostToHiggs);
				//////lokalac_higgs.Boost(BoostToHiggs);

/*
				//*****BOOST TO LAB*****
				lokalac_Niz_el_in4.Boost(BoostToLab);
				lokalac_Niz_poz_in4.Boost(BoostToLab);
				lokalac_Niz_el_out4.Boost(BoostToLab);
				lokalac_Niz_poz_out4.Boost(BoostToLab);
				lokalacZ1.Boost(BoostToLab);
				lokalacZ2.Boost(BoostToLab);
				lokalacZ1_1.Boost(BoostToLab);
				higgs.Boost(BoostToLab);
				////Niz_el_in.Boost(BoostToHiggs);
				////Niz_el_out.Boost(BoostToHiggs);
*/
/*
				//*****BOOST TO FINAL e-e+****
				lokalac_Niz_el_in4.Boost(BoostToFinalee);
				lokalac_Niz_poz_in4.Boost(BoostToFinalee);
				lokalac_Niz_el_out4.Boost(BoostToFinalee);
				lokalac_Niz_poz_out4.Boost(BoostToFinalee);
				lokalacZ1.Boost(BoostToFinalee);
				lokalacZ2.Boost(BoostToFinalee);
				lokalacZ1_1.Boost(BoostToFinalee);
				higgs.Boost(BoostToFinalee);
*/
/*
				//*****BOOST TO Z1Z2*****
				lokalac_Niz_el_in4.Boost(BoostToZ1Z2);
				lokalac_Niz_poz_in4.Boost(BoostToZ1Z2);
				lokalac_Niz_el_out4.Boost(BoostToZ1Z2);
				lokalac_Niz_poz_out4.Boost(BoostToZ1Z2);
				lokalacZ1.Boost(BoostToZ1Z2);
				lokalacZ2.Boost(BoostToZ1Z2);
				lokalacZ1_1.Boost(BoostToZ1Z2);
				lokalac_higgs.Boost(BoostToZ1Z2);
				////Niz_el_in.Boost(BoostToHiggs);
				////Niz_el_out.Boost(BoostToHiggs);
*/
				Niz_el_in3.SetXYZ(lokalac_Niz_el_in4.X(), lokalac_Niz_el_in4.Y(), lokalac_Niz_el_in4.Z());
				Niz_el_out3.SetXYZ(lokalac_Niz_el_out4.X(), lokalac_Niz_el_out4.Y(), lokalac_Niz_el_out4.Z());
		   		Niz_poz_in3.SetXYZ(lokalac_Niz_poz_in4.X(), lokalac_Niz_poz_in4.Y(), lokalac_Niz_poz_in4.Z());
		    		Niz_poz_out3.SetXYZ(lokalac_Niz_poz_out4.X(), lokalac_Niz_poz_out4.Y(), lokalac_Niz_poz_out4.Z());

				//OnShell3_l1.SetXYZ(lokalac_OnShell4_l1.X(), lokalac_OnShell4_l1.Y(), lokalac_OnShell4_l1.Z());
				//OnShell3_l2.SetXYZ(lokalac_OnShell4_l2.X(), lokalac_OnShell4_l2.Y(), lokalac_OnShell4_l2.Z());

				OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z()); // impuls Z1
				OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z()); // impuls Z2
				///OnShell3_Z1.SetXYZ(lokalacZ1_1.X(),lokalacZ1_1.Y(),lokalacZ1_1.Z()); // impuls OnShell Z

 	 	       	        higgsl3.SetXYZ(higgs.X(), higgs.Y(), higgs.Z());

				//------ NASA DEFINICIJA UGLA PHI - u ovom slucaju su normale na ravni usmerene isto ------//

				TVector3 brojilac1 = Niz_el_in3.Cross(Niz_el_out3); // brojilac vektora n1
				Double_t imenilac1 = brojilac1.Mag();   // sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));   // imenilac vektora n1

				TVector3 n1 = brojilac1 * pow(imenilac1,-1);           // vektor n1
				//sin_n1 = sin(Niz_el_in3.Angle(Niz_el_out3));
		        	//sin_n11= brojilac1.Mag() * (pow (Niz_el_in3.Mag(), -1)) * pow(Niz_el_out3.Mag(), -1);

				TVector3 brojilac2 = Niz_poz_in3.Cross(Niz_poz_out3); // brojilac vektora n2
				Double_t imenilac2 = brojilac2.Mag(); // sqrt(pow(brojilac2.X(),2) + pow(brojilac2.Y(),2) + pow(brojilac2.Z(),2));   // imenilac vektora n1

		        	TVector3 n2 = brojilac2 * pow(imenilac2,-1);           // vektor n2
		        	//sin_n2 = sin(Niz_poz_in3.Angle(Niz_poz_out3));
		        	//sin_n21 = brojilac2.Mag() * (pow(Niz_poz_in3.Mag(), -1)) * pow(Niz_poz_out3.Mag(), -1);
		        	//cos_n12 = n1.Dot(n2) * pow(n1.Mag(), -1) * pow(n2.Mag(), -1);

				//fi1 = acos(n1.Dot(n2));      // * 180 / M_PI prvi ugao
				/////OK////Phi = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));


				// MASA HIGGS BOZONA

				double cme = 1000;
				double P_dl_in = (e_starac + p_starac).P();// ukupan impuls para leptona u inicijalnom stanju
				double P_dl_out = (es + ps).P();// ukupan impuls para leptona u finalnom stanju
				double Pt_dl_out = (es + ps).Pt();
				//double P_dl_in = (Niz_el_in + Niz_poz_in).P();// ukupan impuls para leptona u inicijalnom stanju
				//double P_dl_out = (Niz_el_out + Niz_poz_out).P();// ukupan impuls para leptona u finalnom stanju
				double E_dl_in = (e_starac + p_starac).E(); // ukupna energija para leptona u finalnom stanju
				double E_dl_out = (es + ps).E(); // ukupna energija para leptona u finalnom stanju
				double M_dl_in = (e_starac + p_starac).M(); // invarijantna masa para leptona u finalnom stanju
				double M_dl_out = (es + ps).M(); // invarijantna masa para leptona u finalnom stanju
				double diff_E_H = E_dl_in - E_dl_out;
				double diff_P_H = P_dl_in - P_dl_out;

				double Mrec21 = sqrt((E_dl_out * E_dl_out) - (P_dl_out * P_dl_out)); // Masa Z-bozona OK

				//double Mrec2 = ((sqrt(2) - lokalacZ1.E())*(sqrt(2) - lokalacZ1.E())) - (abs(PZ)*abs(PZ));
				double Mrec2 = ((cme * cme) - (2 * cme * E_dl_out) + (M_dl_out * M_dl_out));
				/////double Mrec2 = ((E_dl_out * E_dl_out) - (M_dl_out * M_dl_out));  // Masa H-bozona OK
				//double Mrec2 = diff_E_H * diff_E_H - diff_P_H * diff_P_H;  // Masa H-bozona OK

				//double Mrec2 = 1000*1000 + M_dl_out * M_dl_out - 2 * E_dl_out * 1000;


				////////Mass_higgsl3 = (lokalacZ1 + lokalacZ2).M();

				Mass_higgsl3 = (sqrt(Mrec2));

				if (Mrec21 > 300 && Mrec21 < 870)
				{

				Mass_higgsl3a = (sqrt(Mrec2));

				Mass_Z1 = lokalacZ1.M();
				Mass_Z2 = lokalacZ2.M();

				E_Z1 = lokalacZ1.E();
				E_Z2 = lokalacZ2.E();
				P_Z1 = lokalacZ1.P();
				P_Z2 = lokalacZ2.P();

				//Mass_higgsl3a = - diff_P_H;

				}

				if (Mrec21 > 200)
				{
				Mee_ZZ_fuzija = Mrec21;
				histo_Mee_ZZ_fusion.Fill(Mee_ZZ_fuzija);
				}

				if (Mrec21 < 200)
				{
				Mee_ZH = Mrec21;
				histo_Mee_ZH.Fill(Mee_ZH);
				}

				if (Mrec21 > 200)
				{
				Phi = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));
				histo_Phi_ZZ_fusion.Fill(Phi);
				Tree.Fill();
				}

				histo_Mee_all.Fill(Mrec21);
				// leptonTree.Fill();
				//TanjaTree.Fill();

				//cout << "Ukupan br. leptona1 u kolekciji"<< nlep << endl;

			//TanjaTree.Fill();

			} // if za nlep = 2

			histo_nlep.Fill(nlep);
			histo_npfo.Fill(n_pfo);
			//TanjaTree.Fill();

			/////////n_evt_2lep++;

		 } // if za signal

		//TanjaTree.Fill();
		//leptonTree.Fill();

		} // End of event loop

			cout << "Zatvaram petlju po dogadjajima\n";
	//		cout << "Ukupan br. dogadjaja sa 2 leptona " << nevt2l << endl;
			cout << "Broj dogadjaja po fajlu " << nevent << endl;

			lcReader->close();
 	//       if (right_evt) break;

		} // Kraj petlje po fajlovima

		 ///// ptTree.Fill();

	cout << "Ukupan broj dogadjaja " << nevent_uk << endl;
	//cout << "Ukupan br. dogadjaja sa 2 leptona "<< n_evt_2lep << endl;
	cout << "Ukupan br. dogadjaja signala "<< n_signal_uk << endl;

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	Tree.Write();
	histo_y12.Write();
	histo_y23.Write();
	histo_nlep.Write();
	histo_npfo.Write();
	histo_Mee_ZZ_fusion.Write();
	histo_Mee_ZH.Write();
	histo_Phi_ZZ_fusion.Write();
        histo_Mee_all.Write();
	//leptonTree.Write();
	//ptTree.Write();
	rootFile.Write();
	rootFile.Close();

return 0;

} // Kraj funkcije

Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;
	UInt_t nFirstJob = 1;
	if(argc>iarg) nFirstJob = atoi(argv[iarg]); iarg++;
	UInt_t nLastJob = 1;
	if(argc>iarg) nLastJob  = atoi(argv[iarg]); iarg++;

	TString fName = "signal_1tev_";
	if(argc>iarg) fName = argv[iarg]; //iarg++;

	TString rfName = "signal_1tev.root";
	if(argc>iarg) rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
