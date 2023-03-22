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

	Float_t m_onshell, m_offshell, m_ll, m_qq, m_H;
	Float_t m_Z1, m_Z2, p_e1, p_e2, E_e1, E_e2, pt_q1, pt_q2, E_q1, E_q2;
	Float_t n_signal, n_signal_uk, n_evt_mva;
	Float_t E_vis, E_vis_E_H,pt_H, E_H;
	Float_t pt_e_sistema, pt_p_sistema, pt_e1, pt_e2;
	Float_t n_pfo, pt_miss, theta_H, theta_el, theta_poz, theta_el_reko, theta_poz_reko;
	Float_t logy_12, logy_23, logy_34;
	Float_t Btag1, Btag2;
	Float_t Phi;

	//Tree.Branch("m_H", &m_H, "m_H");
	//Tree.Branch("m_qq", &m_qq, "m_qq"); // invarijatne mase b-dzeta???
	//Tree.Branch("m_ll", &m_ll, "m_ll"); // invarijatne mase b-dzeta???
	Tree.Branch("n_signal", &n_signal, "n_signal");
	Tree.Branch("n_signal_uk", &n_signal_uk, "n_signal_uk");
	Tree.Branch("n_evt_mva", &n_evt_mva, "n_evt_mva");
	//Tree.Branch("nlep", &nlep, "nlep");

	Tree.Branch("m_qq", &m_qq, "m_qq"); // invarijatne mase b-dzeta???
	Tree.Branch("m_H", &m_H, "m_H"); // masa Higzovog bozona
	Tree.Branch("m_Z1", &m_Z1, "m_Z1"); // masa 1. onshell Z-bozonka
	Tree.Branch("m_Z2", &m_Z2, "m_Z2"); // masa 2. onshell Z-bozonka
	Tree.Branch("p_e1", &p_e1, "p_e1"); // impuls elektrona
	Tree.Branch("p_e2", &p_e2, "p_e2"); // impuls pozitrona
	Tree.Branch("E_e1", &E_e1, "E_e1"); // energija elektrona
	Tree.Branch("E_e2", &E_e2, "E_e2"); // energija pozitrona
    Tree.Branch("pt_e1", &pt_e1, "pt_e1"); // energija elektrona
	Tree.Branch("pt_e2", &pt_e2, "pt_e2"); // energija pozitrona
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
	//Tree.Branch("logy_12", &logy_12, "logy_12");
	//Tree.Branch("logy_23", &logy_23, "logy_23");
	Tree.Branch("Phi", &Phi, "Phi");

/*
	Float_t sin_n1, sin_n2;
	Float_t sin_n11, sin_n21, cos_n12;
	Float_t cos_ksi, sin_ksi, sin_ksi1, tan_ksi, tan_ksi1;
	Float_t fi1, fi2, fi3, fi4, fi5, fi6, fi7, fiOgawa, fiOgawa1, fiOgawa2, fiOgawa3, fiOgawa4;
	Float_t theta1, theta2, theta3, theta4, thetaH, thetaH1, theta1a, polarni_ugao, theta_tracker_el_in, theta_tracker_el_out, theta_el_all, theta_tracker_poz_in, theta_tracker_poz_out, theta_poz_all, thetael2, thetael3, thetael2st, thetael3st, PI, theta_tracker_all_in, theta_tracker_all_out, theta_tracker_all_in2, theta_tracker_all_out2;
	Float_t m_Z1, m_Z2;
	Float_t Pt, Px_el_in, Py_el_in, Pz_el_in;
	Float_t Px_poz_in, Py_poz_in, Pz_poz_in;
	Float_t Px_el_out, Py_el_out, Pz_el_out, Pt_el_out;
	Float_t Px_poz_out, Py_poz_out, Pz_poz_out, Pt_poz_out;
	Float_t Px_lokalacZ1, Py_lokalacZ1, Pz_lokalacZ1, P_lokalacZ1, M_lokalacZ1, M_lokalacZ2;
	Float_t thetaZ11, thetaZ22;
	Float_t Px_higgsl3, Py_higgsl3, Pz_higgsl3, P_higgsl3, P_higgsl3_1, P_higgsl3_2, P_higgsl3_3;
	Float_t Mass_higgsl3, Mee;
	Float_t Ene_el_out, Ene_el_out_all;

	Tree.Branch("#phi_{1}", &fi1, "fi1");
	Tree.Branch("#phi_{2}", &fi2, "fi2");
	Tree.Branch("#phi_{3}", &fi3, "fi3");
	Tree.Branch("#phi_{4}", &fi4, "fi4");
	Tree.Branch("#phi_{5}", &fi5, "fi5");
	Tree.Branch("#phi_{6}", &fi6, "fi6");
	Tree.Branch("#phi_{7}", &fi7, "fi7");
	Tree.Branch("#phi_{Ogawa}", &fiOgawa, "fi6");
	Tree.Branch("#phi_{Ogawa1}", &fiOgawa1, "fi7");
	Tree.Branch("#phi_{Ogawa2}", &fiOgawa2, "fi8");
	Tree.Branch("#phi_{Ogawa3}", &fiOgawa3, "fi9");
	Tree.Branch("#phi_{Ogawa4}", &fiOgawa4, "fi10");
	Tree.Branch("sin_n1", &sin_n1, "sin_n1");
	Tree.Branch("sin_n2", &sin_n2, "sin_n2");
	Tree.Branch("sin_n11", &sin_n11, "sin_n11");
	Tree.Branch("polarni_ugao", &polarni_ugao, "polarni_ugao");
	Tree.Branch("#theta_{H1}", &thetaH1, "theta_{H1}");
	Tree.Branch("#theta_1", &theta1, "#theta_1");
	Tree.Branch("#theta_1a", &theta1a, "#theta_1a");
	Tree.Branch("theta2", &theta2, "theta2");
	Tree.Branch("theta_el_all", &theta_el_all, "theta_el_all");
	Tree.Branch("theta_tracker_el_in", &theta_tracker_el_in, "theta_tracker_el_in");
	Tree.Branch("theta_tracker_el_out", &theta_tracker_el_out, "theta_tracker_el_out");
	Tree.Branch("theta_poz_all", &theta_poz_all, "theta_poz_all");
	Tree.Branch("theta_tracker_poz_in", &theta_tracker_poz_in, "theta_tracker_poz_in");
	Tree.Branch("theta_tracker_poz_out", &theta_tracker_poz_out, "theta_tracker_poz_out");
	Tree.Branch("theta_tracker_all_in", &theta_tracker_all_in, "theta_tracker_all_in");
	Tree.Branch("theta_tracker_all_out", &theta_tracker_all_out, "theta_tracker_all_out");
	Tree.Branch("theta_tracker_all_in2", &theta_tracker_all_in2, "theta_tracker_all_in2");
	Tree.Branch("theta_tracker_all_out2", &theta_tracker_all_out2, "theta_tracker_all_out2");
	Tree.Branch("#theta_3", &theta3, "#theta_3");
	Tree.Branch("#theta_{poz}", &theta4, "#theta_{poz}");
	Tree.Branch("#theta_Z_1", &thetaZ11, "theta_Z_1");
	Tree.Branch("#theta_Z_2", &thetaZ22, "theta_Z_2");
	Tree.Branch("Pt", &Pt, "Pt");
	Tree.Branch("Px_el_in", &Px_el_in, "Px_el_in");
	Tree.Branch("Py_el_in", &Py_el_in, "Py_el_in");
	Tree.Branch("Pz_el_in", &Pz_el_in, "Pz_el_in");
	Tree.Branch("Px_el_out", &Px_el_out, "Px_el_out");
	Tree.Branch("Py_el_out", &Py_el_out, "Py_el_out");
	Tree.Branch("Pz_el_out", &Pz_el_out, "Pz_el_out");
	Tree.Branch("Px_poz_in", &Px_poz_in, "Px_poz_in");
	Tree.Branch("Py_poz_in", &Py_poz_in, "Py_poz_in");
	Tree.Branch("Pz_poz_in", &Pz_poz_in, "Pz_poz_in");
	Tree.Branch("Px_poz_out", &Px_poz_out, "Px_poz_out");
	Tree.Branch("Py_poz_out", &Py_poz_out, "Py_poz_out");
	Tree.Branch("Pz_poz_out", &Pz_poz_out, "Pz_poz_out");
	Tree.Branch("Px_lokalacZ1", &Px_lokalacZ1, "Px_lokalacZ1");
	Tree.Branch("Py_lokalacZ1", &Py_lokalacZ1, "Py_lokalacZ1");
	Tree.Branch("Pz_lokalacZ1", &Pz_lokalacZ1, "Pz_lokalacZ1");
	Tree.Branch("P_lokalacZ1", &P_lokalacZ1, "P_lokalacZ1");
	Tree.Branch("M_lokalacZ1", &M_lokalacZ1, "M_lokalacZ1");
	Tree.Branch("M_lokalacZ2", &M_lokalacZ2, "M_lokalacZ2");
	Tree.Branch("Px_higgsl3", &Px_higgsl3, "Px_higgsl3");
	Tree.Branch("Py_higgsl3", &Py_higgsl3, "Py_higgsl3");
	Tree.Branch("Pz_higgsl3", &Pz_higgsl3, "Pz_higgsl3");
	Tree.Branch("P_higgsl3", &P_higgsl3, "P_higgsl3");
	Tree.Branch("P_higgsl3_1", &P_higgsl3_1, "P_higgsl3_1");
	Tree.Branch("P_higgsl3_2", &P_higgsl3_2, "P_higgsl3_2");
	Tree.Branch("P_higgsl3_3", &P_higgsl3_3, "P_higgsl3_3");
	Tree.Branch("M_{H}", &Mass_higgsl3, "M_{H}");
	Tree.Branch("Mee", &Mee, "Mee");
	Tree.Branch("E_{tracker}", &Ene_el_out, "E_{tracker}");
	Tree.Branch("E_{all}", &Ene_el_out_all, "E_{all}");
	Tree.Branch("thetael2", &thetael2, "thetael2");
	Tree.Branch("thetael3", &thetael3, "thetael3");
	Tree.Branch("thetael2st", &thetael2st, "thetael2st");
	Tree.Branch("thetael3st", &thetael3st, "thetael3st");
	Tree.Branch("PI", &PI, "PI");
	//Tree.Branch("PDG", &PDG, "PDG");
*/
	//Tree.Branch("n_signal_uk", &n_signal_uk, "n_signal_uk");
	//Tree.Branch("n_evt_mva", &n_evt_mva, "n_evt_mva");

	TH1F histo_y12 ("Histogram_y12", " ; -log y_{12}", 100, 0, 10);
	TH1F histo_y23 ("Histogram_y23", " ; -log y_{23}", 100, 0, 10);
	TH1F histo_y34 ("Histogram_y34", " ; -log y_{34}", 100, 0, 10);

	TH1F histo_nlep ("Histogram_nlep", " ; nlep", 100, 0, 2000);
	TH1F histo_npfo ("Histogram_npfo", " ; n_pfo", 100, 0, 10000);

	TH1F histo_theta_el ("theta_el", " ; theta_el", 100, 0, 180); // Polarni ugao prvog izracenog elektrona posle sudara
	TH1F histo_theta_poz ("theta_poz", " ; theta_poz", 100, 0, 180); // Polarni ugao prvog izracenog pozitrona posle sudara
	TH1F histo_theta_el_reko ("theta_el_reko", " ; theta_el_reko", 100, 0, 180); // Polarni ugao rekonstruisanog, krajnjeg elektrona
	TH1F histo_theta_poz_reko ("theta_poz_reko", " ; theta_poz_reko", 100, 0, 180); // Polarni ugao rekonstruisanog, krajnjeg pozitrona

	TH1F histo_Phi_ZZ_fusion ("h_Phi_ZZ_fusion", " ; #phi", 100, -3.7, 3.7);

	Int_t nevent_uk = 0;


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

		while( (evt = lcReader -> readNextEvent()) != 0)
		{

			nevent++;
			nevent_uk ++;

			Int_t n_isolep_evt = 0;
			Int_t n_signal = 0;

			std::vector<std::string> colNames = *evt -> getCollectionNames();
			vector <TLorentzVector> nizq_mc;

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");
			IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*)evt->getCollection("FJ_Jets_2");
			IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*)evt->getCollection("Isolep_Selected");
		//	IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*)evt->getCollection("Isolep_Selected_new");
			IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*)evt->getCollection("PandoraPFOs");
			//IMPL::LCCollectionVec* colJet = (IMPL::LCCollectionVec*)evt->getCollection("RefinedJets");

			Int_t nqvarkova = 0;
			bool signal = false;
			TLorentzVector e_starac, p_starac; //dodao 4-vektore elektrona i pozitrona koji se sudaraju
			Float_t visibleEnergy = 0;

// JAJAJAJAJAAJ
			// Petlja preko koje prolazimo kroz sve cestice po svakom dogadjaju
		   	for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		   	{

				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				if (mcParticle -> getPDG() == 11 && mcParticle -> getParents().size() == 0)
				{
					const EVENT::MCParticleVec & cerke = mcParticle -> getDaughters();
					for (Int_t h = 0; h < cerke.size(); h++)
					{
						Int_t PDG_starac = cerke[h] -> getPDG();

						const double *pgen = cerke[h] -> getMomentum(); // impuls čestice
						double egen = cerke[h] -> getEnergy();	//energija čestice

						// cout << "Energija: " << parent[h] -> getEnergy() << endl;

						if (PDG_starac == 11 ) e_starac.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen); // KG definisana imena u liniji 157
						if (PDG_starac == -11 ) p_starac.SetPxPyPzE(pgen[0], pgen[1], pgen[2], egen);

					}//end za roditelje 163

				}///end if za čestice koje se sudaraju

			}

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
			//		cout << "Ukupan br. pfos u kolekciji" << n_pfo << endl;

					pt_miss = sum_vec.Pt();

					vector<TLorentzVector> pfo_jets;
					vector<TLorentzVector> iso_lep;

				Float_t m_y12 = 0;
				Float_t m_y23 = 0;

				for (Int_t i = 0; i < jets2 -> getNumberOfElements(); i++)
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

		//		cout << "Ukupan br. jets u kolekciji"<< nqu << endl;

				histo_y12.Fill(m_y12);
				histo_y23.Fill(m_y23);

				for (Int_t j = 0; j < isolep->getNumberOfElements(); j++)
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

				//	cout << "pid" << pid << endl;

				//cout << "Ukupan br. leptona u kolekciji"<< nlep << endl;

				}

		//		cout << "Ukupan br. leptona1 u kolekciji"<< nlep << endl;

				//histo_nlep.Fill(nlep);

			theta_el = e_starac.Theta() * 180/M_PI;
			theta_poz = p_starac.Theta() * 180/M_PI;

			histo_theta_el.Fill(theta_el);
			histo_theta_poz.Fill(theta_poz);

			if (nlep == 2) //&& theta_el > 8 && theta_poz < 172)
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

				TLorentzVector Z1;
				TLorentzVector Z2;

				Z1 = e_starac + es;
				Z2 = p_starac + ps;

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
				 m_H = higgs.M();
				 E_H = higgs.E();
				 E_vis_E_H = abs (E_vis - E_H);
				 theta_H = higgs.Theta() * 180/M_PI;
				 pt_H = higgs.Pt(); // KG ovo dodati u drvo
				 theta_el_reko = es.Theta() * 180/M_PI;
				 theta_poz_reko = ps.Theta() * 180/M_PI;
				 //logy_12 = m_y12;
				 //logy_23 = m_y23;

				 histo_theta_el_reko.Fill(theta_el_reko);
				 histo_theta_poz_reko.Fill(theta_poz_reko);

				//DEFINICIJA UGLA PHI -->> Tanja

				TVector3 BoostToHiggs = -(higgs.BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
				//TVector3 BoostToHiggs = -((bkvark1 + bkvark2).BoostVector());
				TVector3 BoostToLab = -((es+higgs+ps).BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona

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
				Phi = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));

				histo_Phi_ZZ_fusion.Fill(Phi);

				// leptonTree.Fill();
				Tree.Fill();

			} // if za nlep = 2

			histo_nlep.Fill(nlep);
			histo_npfo.Fill(n_pfo);
			//Tree.Fill();

			/////////n_evt_2lep++;

			// } // if za signal

		//Tree.Fill();
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
	histo_theta_el.Write();
	histo_theta_poz.Write();
	histo_theta_el_reko.Write();
	histo_theta_poz_reko.Write();
	histo_Phi_ZZ_fusion.Write();
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

	TString fName = "/mnt/myBook/ilc_cpv/aa_4f/qqqq/backg_aa_4f_eB_pB_1tev_";
	//TString fName = "/media/tjovin/VERBATIM/ILC_ZZ_fusion_ana/data/background_files/qqee_1tev/backg_qqee_1tev_";
	//TString fName = "/media/tjovin/VERBATIM/ILC_ZZ_fusion_ana/data/background_files/qqenu_1tev/qqenu_1tev_jets_1.1_";
	if(argc>iarg) fName = argv[iarg]; //iarg++;

	TString rfName = "aa_qqqq_eB_pB.root"; //"backg_aa_4f_qqqq_1tev.root"; // backg_qqenu_1tev.root // backg_qqee_1tev.root
	if(argc>iarg) rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
