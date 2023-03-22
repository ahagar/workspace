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
	#include <UTIL/PIDHandler.h>

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

	TTree ptTree ("ptTree", "Generator particle tree");
	TTree pfoTree ("pfoTree", "Generator particle tree");
	TTree mcTree ("mcTree", "Generator particle tree"); // definicija drveta gde upisujemo histograme
	TTree elektronskoTree("elektronskoTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za sve elektrone i pozitrone
	TTree elektronsko6i7Tree("elektronsko6i7Tree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami samo za elektron 6 i pozitron 7
	TTree elektronskoRekonstrTree("elektronskoRekonstrTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za rekonstruisane elektrone i pozitrone
//	TTree pozitronskoRekonstrTree("pozitronskoRekonstrTree", "Generator particle tree"); // ovde su Kragujevcani definisali drvce gde se upisuju histogrami za rekonstruisane pozitrone
//	elektronskoRekonstrTree.AddFriend("pozitronskoRekonstrTree", "");
	Float_t m_ll, m_qq, m_H, m_Z1, m_Z2, p_e1, p_e2, E_e1, E_e2, pt_q1, pt_q2, E_q1, E_q2, pt_e1, pt_e2;
	Float_t E_Z1, E_Z2, theta_Z1, theta_Z2, m_e1e2, theta_e, theta_p;
	Float_t E_vis, E_vis_E_H,pt_H, E_H, E_miss;
	Float_t pt_e_sistema, pt_p_sistema;
	Float_t n_pfo, pt_miss, theta_H;
	Float_t logy_12, logy_23, logy_34;
	Float_t Btag1, Btag2;
	// Float_t Ctag1, Ctag2;
///////////////////////////////////////////////
	Float_t fi, theta1, theta2;
	Float_t n_signal, n_signal_uk, n_evt_mva, nevent_uk;
	TTree leptonTree ("leptonTree", "Generator particle tree");
	leptonTree.Branch("m_qq", &m_qq, "m_qq"); // invarijatne mase b-dzeta???
	leptonTree.Branch("m_H", &m_H, "m_H"); // masa Higzovog bozona
	leptonTree.Branch("m_Z1", &m_Z1, "m_Z1"); // masa 1. onshell Z-bozonka
	leptonTree.Branch("m_Z2", &m_Z2, "m_Z2"); // masa 2. onshell Z-bozonka
	leptonTree.Branch("E_Z1", &E_Z1, "E_Z1");
	leptonTree.Branch("E_Z2", &E_Z2, "E_Z2");
	leptonTree.Branch("theta_Z1", &theta_Z1, "theta_Z1");
	leptonTree.Branch("theta_Z2", &theta_Z2, "theta_Z2");
	leptonTree.Branch("m_e1e2", &m_e1e2, "m_e1e2");
	leptonTree.Branch("p_e1", &p_e1, "p_e1"); // impuls elektrona
	leptonTree.Branch("p_e2", &p_e2, "p_e2"); // impuls pozitrona
	leptonTree.Branch("theta_e", &theta_e, "theta_e"); // impuls elektrona
	leptonTree.Branch("theta_p", &theta_p, "theta_p"); // impuls pozitrona
	leptonTree.Branch("pt_e1", &pt_e1, "pt_e1"); // impuls elektrona
	leptonTree.Branch("pt_e2", &pt_e2, "pt_e2"); // impuls pozitrona
	leptonTree.Branch("E_e1", &E_e1, "E_e1"); // energija elektrona
	leptonTree.Branch("E_e2", &E_e2, "E_e2"); // energija pozitrona
	leptonTree.Branch("pt_e_sistema", &pt_e_sistema, "pt_e_sistema");
	leptonTree.Branch("pt_p_sistema", &pt_p_sistema, "pt_p_sistema");
	leptonTree.Branch("pt_q1", &pt_q1, "pt_q1"); // transverzalni impuls 1. kvarka
	leptonTree.Branch("pt_q2", &pt_q2, "pt_q2"); // transverzalni impuls 2. kvarka
	leptonTree.Branch("E_q1", &E_q1, "E_q1"); // energija 1. kvarka
	leptonTree.Branch("E_q2", &E_q2, "E_q2"); // energija2. kvarka
	leptonTree.Branch("theta_H", &theta_H, "theta_H"); // Polarni ugao Higzovog bozona
	leptonTree.Branch("E_vis", &E_vis, "E_vis"); // Vidljiva energija svih cestica dzumle
	leptonTree.Branch("E_vis_E_H", &E_vis_E_H, "E_vis_E_H"); // Razlika Vidljive energije svih cestica dzumle i Energije Higzovog bozona
	leptonTree.Branch("E_H", &E_H, "E_H"); // Energija Higzovog bozona
	leptonTree.Branch("pt_H", &pt_H, "pt_H"); // transverzalni impuls Higzovog bozona
	leptonTree.Branch("pt_miss", &pt_miss, "pt_miss"); // nedostajucí transverzalni impuls (otisao neutrinu?)
	leptonTree.Branch("E_miss", &E_miss, "E_miss"); // nedostajucí transverzalni impuls (otisao neutrinu?)
	leptonTree.Branch("logy_12", &logy_12, "logy_12");
	leptonTree.Branch("logy_23", &logy_23, "logy_23");
	leptonTree.Branch("n_pfo", &n_pfo, "n_pfo");
	leptonTree.Branch("Btag1", &Btag1, "Btag1"); // Verovatnoca da je cestica b-kvark 1
	leptonTree.Branch("Btag2", &Btag2, "Btag2"); // Verovatnoca da je cestica b-kvark 2
	leptonTree.Branch("fi", &fi, "fi");




//	Float_t brojac_preostalih_KG_elektrona_6 = 0;
	Float_t pt_KG_elektrona_6 = 0, p_KG_elektrona_6 = 0;
	Float_t Energija_KG_elektrona_6 = 0;



	// Promenljive za PFO rekonstruisanih elektrona/pozitrona
//	Float_t brojac_preostalih_KG_elektrona = 0;


	// Brojacke promenljive za KRAJNJI elektron/pozitron
	Int_t N_ukupnih_dogadjaja = 0;
	Int_t N_KG_elektrona = 0;

	Int_t nsignala_mc = 0;

	Int_t N_za_dva = 0;




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
		while( (evt = lcReader -> readNextEvent()) != 0 && N_dogadjaja <= 31000)//31e+03
		{
			// KG: ovde su Kragujevcani definisali unutrasnje brojace
			int c_n_pfo = 0;
			N_dogadjaja++;
			N_ukupnih_dogadjaja++;
			TLorentzVector temp_e; // definisan privremeni četvorovektor u koji se zapisuju impuls i energija svakog pozitrona
			TLorentzVector temp_p;
			EVENT::ReconstructedParticle* reco_link_lepton = 0;

			//cout <<"Dogadjaj "<< N_dogadjaja << endl;
			Int_t N_leptona_energija_ugao_po_dogadjaju = 0;
			vector <TLorentzVector> leptoni;
			vector <EVENT::ReconstructedParticle*> niz_leptons, niz_pfos_GMD;

	/*	    	bool elektronac = false;
			bool pozitronac = false;

		    	bool elektron2 = false;
			bool pozitron2 = false;

		    	bool elektron2cutted = false;
			bool pozitron2cutted = false;

		    	bool mcElektron = false;
			bool mcPozitron = false;

			double ekp = 0 ;
			double eke = 0;*/

			vector <TLorentzVector> vec_2isol;

			vector <EVENT::ReconstructedParticle*> niz_linklep, niz_pfos;

			// KG: ovde su Kragujevcani definisali potrebne nizove kvadri-vektora
			vector <TLorentzVector> niz_KG_elektrona, niz_KG_pozitrona, niz_KG_kvarkova;

//			vector <TLorentzVector> nizl, nizq;

			std::vector<std::string> colNames = *evt -> getCollectionNames();

		//	EVENT::LCCollection* links = evt -> getCollection("RecoMCTruthLink");

			/*for(int i = 0; i < colNames.size(); i++)
			{
				if (colNames[i] != "RecoMCTruthLink") break;

			}
			cout <<"kolekcija "<< colNames[0] << endl;*/
			//if (colNames != colNames) continue;
			//if (!links) continue;

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticles");//MCParticlesSkimmed, MCParticles
			IMPL::LCCollectionVec* rec_Particles = (IMPL::LCCollectionVec*) evt -> getCollection("PFOs");//PFOs, IsolatedElectrons
			IMPL::LCCollectionVec* recParticles = (IMPL::LCCollectionVec*) evt -> getCollection("IsolatedElectrons");//PFOs,
			IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*)evt -> getCollection("Durham2Jets");

			/*UTIL::PIDHandler pidHandler(evt -> getCollection("PFOs"));
			int trackParamsAlgId = pidHandler.getAlgorithmID("TrackParameters");
			int d0ParamIndex = pidHandler.getParameterIndex(trackParamsAlgId, "D0");*/


			bool signal = false;
			vector <EVENT::ReconstructedParticle*> reclep_vector; // Definisan niz u koji se zapisuju rekonstrusani elektroni/pozitroni

			// Definisani broj elektrona, broj pozitrona i broj kvarkova
//			Int_t N_KG_kvarkova = 0;

/*			Int_t Broj_leptona = 0;
			Int_t Broj_kvarkova = 0; */

//			Float_t Brojac_za_dvicu = 0;
			EVENT::MCParticle* e_starac ;
			EVENT::MCParticle* p_starac ;
			TLorentzVector e_mc_final, p_mc_final;
			int c_e_mc_final = 0, c_p_mc_final = 0;
			int cElec = 0;
			int cPosi = 0;

			bool hbb = false;
			bool e_tracker = false;
			bool p_tracker = false;
			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
				const double *p = mcParticle->getMomentum(); // impuls čestice
				double e = mcParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
				Int_t particlePDG = mcParticle->getPDG();


				if (mcParticle->getGeneratorStatus() == 1) {
					if(particlePDG == 11 && e > 60) {
						e_mc_final = temp;
						c_e_mc_final++;
					}
					if(particlePDG == -11 && e > 60) {
						p_mc_final = temp;
						c_p_mc_final++;
					}

				}



				if (mcParticle -> getPDG() == 25 && mcParticle ->getDaughters().size() == 2 && abs (mcParticle ->getDaughters()[0] -> getPDG()) == 5   )
				{// if za hbb
					hbb = true;
					if (mcParticle->getParents()[0]->getPDG()==25){
					//	cout << "ćale mi je Higgs, a stric mi je:  "<< mcParticle->getParents()[0]->getParents()[0]->getDaughters()[0]->getPDG()<<endl;
					TLorentzVector e_6, p_7;
					EVENT::MCParticle* elecTemp =  mcParticle->getParents()[0]->getParents()[0]->getDaughters()[0];
					EVENT::MCParticle* posiTemp =  mcParticle->getParents()[0]->getParents()[0]->getDaughters()[1];
					e_starac = mcParticle->getParents()[0]->getParents()[0];
					p_starac = mcParticle->getParents()[0]->getParents()[1];


				//	LCCollection* trackCollection = event->getCollection(trackCollectionName);
				//	Track* track = dynamic_cast<Track*>(trackCollection->getElementAt(i));

			//		cout << "ele: "<< e_starac->getPDG()<<endl; //proverio da li je cestica elektron
				//	cout << "posi: "<< p_starac->getPDG()<<endl; // proverio da li je cestica pozitron


					const double *ep = elecTemp->getMomentum();
					double ee = elecTemp->getEnergy();
					const double *pp = posiTemp->getMomentum();
					double pe = posiTemp->getEnergy();

					e_6.SetPxPyPzE(ep[0], ep[1], ep[2], ee);
					p_7.SetPxPyPzE(pp[0], pp[1], pp[2], pe);

					if (e_6.Theta()*180/M_PI > 8 && e_6.Theta()*180/M_PI < 172) e_tracker = true;
					if (p_7.Theta()*180/M_PI > 8 && p_7.Theta()*180/M_PI < 172) p_tracker = true;

			//		cout<<"redni broj dogadjaja: "<<N_dogadjaja <<endl;

			//		cout << "polarni ugao elektrona je : " <<e_6.Theta()*180/M_PI<<endl;
			//		cout << "polarni ugao pozitrona je : " <<p_7.Theta()*180/M_PI<<endl;
				/*	cout << "energija  elektrona je : " <<e_6.E()<<endl;
					cout << "energija   pozitrona je : " <<p_7.E()<<endl;*/
				//	cout <<"----------------------------------------"<<endl;
					} // end of mcParticle->getParents()[0]->getPDG()==25
				}// if za hbb

			}// end of MCParticles

			if (hbb == false || e_tracker == false || p_tracker == false) continue; // izdvajanje signala
			nsignala_mc++;
		//
		//	cout << "redni broj dogadjaja: "<<N_dogadjaja<<endl;

			E_vis = 0;
			pt_miss = 0;
			TLorentzVector sum_vec;
			for (Int_t i = 0; i < rec_Particles -> getNumberOfElements(); i++)
			{
				IMPL::ReconstructedParticleImpl* rec_Particle = (IMPL::ReconstructedParticleImpl*) rec_Particles -> getElementAt(i);
				IMPL::ReconstructedParticleImpl* iso_lep;
				c_n_pfo++;
			//	cout << c_n_pfo<<endl;
			    E_vis += rec_Particle->getEnergy();
			    TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
			    const double *p = rec_Particle->getMomentum(); // impuls čestice
				double e = rec_Particle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
			    sum_vec += temp;
			//    int no_of_elements = rec_Particles->getNumberOfElements();

			   if (abs(rec_Particle->getType()) ==11)  iso_lep = rec_Particle;

		//		UTIL::PIDHandler pidHandler(evt -> getCollection("PFOs"));


			 /*   int algo = pidHandler.getAlgorithmID( "TrackParameters" );// output 1
			    float d0tag = pidHandler.getParameterIndex(algo, "D0");// output 2
			    float d0_par = pidHandler.getParticleID(rec_Particle, algo).getParameters()[d0tag];
				//auto particleId = 	pidHandler.getParticleID(rec_Particle, algo).getParameters()[d0tag];
			//    cout <<N_dogadjaja<<". , "<<"size: "<<pidHandler.getParticleID(rec_Particle, algo).getParameters().size()<<endl;
			    cout <<N_dogadjaja<<". , "<<"D0: "<<d0_par<<endl;*/
			/*	if (particleId){

				}*/

			/*	UTIL::PIDHandler pidHandler(evt -> getCollection("PFOs"));
				int trackParamsAlgId = pidHandler.getAlgorithmID("TrackParameters");
				int d0ParamIndex = pidHandler.getParameterIndex(trackParamsAlgId, "D0");*/


		//	    cout << "D0 ???" << rec_Particle->getTracks().size() << endl;

	     	   //  float D0 = static_cast<double>(rec_Particles -> parameters().getFloatVal("D0"));
			   //cout <<"d0: "<< static_cast<float>(jets2 -> parameters().getFloatVal("y_{n-1,n}"))<<endl;
			//	cout <<" D0: "<<endl;

			    //auto particleIds = rec_Particle->getParticleIDs();
			    //cout << "particleIds.size(): " << particleIds.size() << endl;
			    //float flvtag1 = pidHandler.getParticleID(particleIds[0], trackParamsAlgId).getParameters()[d0ParamIndex];

			/*    auto particleIds = pidHandler.getParticleIDs(rec_Particle->getParticleIDUsed(), trackParamsAlgId);
			    if (particleIds.size() > 0){
			    	float test = particleIds[0]-getParameters()[d0ParamIndex];
			    	cout << "test --> " << test << endl;
			    }


			    EVENT::StringVec recPartParams;
			    rec_Particles->parameters().getStringVals("ParameterNames_TrackParameters", recPartParams);
			    for(int p = 0;  p < (int)recPartParams.size(); p++) {
			   // 	cout <<" D0[" << p << "]: " << recPartParams[p] << endl;
			    }*/

			}
		    pt_miss += sum_vec.Pt();


			Int_t N_KG_elektrona_6 = 0;
			Int_t N_KG_pozitrona_7 = 0;

			TLorentzVector e_rec, p_rec;
			for (Int_t i = 0; i < recParticles -> getNumberOfElements(); i++)
			{
				IMPL::ReconstructedParticleImpl* recParticle = (IMPL::ReconstructedParticleImpl*) recParticles -> getElementAt(i);
				TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
				const double *p = recParticle->getMomentum(); // impuls čestice
				double e = recParticle->getEnergy();	//energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor



				//cout <<"Nasao elektron \n";

				//if (recParticle -> getType() == 11) {Energija_KG_elektrona_6 = recParticle -> getEnergy(); elektronsko6i7Tree.Fill(); }
				if (recParticle -> getEnergy() > 60)
				{
					if (recParticle -> getType() == 11)
					{N_KG_elektrona_6 ++; Energija_KG_elektrona_6 = recParticle -> getEnergy();
					e_rec = temp;
					cElec++;
					}

					if (recParticle -> getType() == -11)
					{N_KG_pozitrona_7 ++;
					p_rec = temp;
					cPosi++;
					/*cout <<"nasao pozitron\n";*/}
				}

			} // for po PFOs
			/*cout << "Ukupan broj pocetnih elektrona 6: " << N_KG_elektrona_6 << endl;
			cout << "Ukupan broj pocetnih pozitrona 7: " << N_KG_pozitrona_7 << endl;
			cout <<"____________________________________________________________"<<endl;*/

		/*	cout <<N_dogadjaja<< ". po dogadjaju broj MC elektrona: "<<c_e_mc_final<< ", broj reco IsoLep elektona: "<< cElec<<endl;
			cout << "po dogadjaju broj MC pozitrona: "<<c_e_mc_final<< ", broj reco IsoLep pozitrona: "<< cPosi<<endl;

			bool e_t = false, p_t = false;
			if (abs (e_mc_final.Theta()*180/M_PI -e_rec.Theta()*180/M_PI ) < 0.1) e_t = true;
			if (abs (p_mc_final.Theta()*180/M_PI -p_rec.Theta()*180/M_PI ) < 0.1) p_t = true;

			if (e_t) cout << "theta MC elektrona: "<< e_mc_final.Theta()*180/M_PI<< ", isolovanog elektrona: " << e_rec.Theta()*180/M_PI<<endl;

			if (p_t) cout << "theta MC pozitrona: "<< p_mc_final.Theta()*180/M_PI<< ", isolovanog pozitrona: " << p_rec.Theta()*180/M_PI<<endl;
			cout <<"----------------------------------------"<<endl;*/



			vector <TLorentzVector> v1l, v2j;
			double m_y12;
			double m_y23;
			 for (Int_t i = 0; i < jets2 -> getNumberOfElements(); i++)
			 {
				 IMPL::ReconstructedParticleImpl* recJet = (IMPL::ReconstructedParticleImpl*) jets2 -> getElementAt(i);
				 TLorentzVector tempjet;
				// nqu++;
				 const double *p = recJet -> getMomentum();
				 Double_t e = recJet -> getEnergy();
				 tempjet.SetPxPyPzE(p[0], p[1], p[2], e);
				 Int_t pidjet = recJet -> getType();
				 v2j.push_back(tempjet);

				 float yMinus = static_cast<double>(jets2 -> parameters().getFloatVal("y_{n-1,n}"));
				 float yPlus  = static_cast<double>(jets2 -> parameters().getFloatVal("y_{n,n+1}"));
				 float dMinus = static_cast<double>(jets2 -> parameters().getFloatVal("d_{n-1,n}"));
				 float dPlus  = static_cast<double>(jets2 -> parameters().getFloatVal("d_{n,n+1}"));
				// jets2->parameters().


				 m_y12 = -log10(yMinus);
				 m_y23 = -log10(yPlus);
			//	 cout << "y12 " << m_y12<<endl;
			//	 cout << "y23 " << m_y23<<endl;

			 }


			if(N_KG_elektrona_6 == 1 && N_KG_pozitrona_7 == 1) {

				TLorentzVector e_s, p_s;
				const double *e_es = e_starac->getMomentum();
				double e_e = e_starac->getEnergy();
				const double *p_ps = p_starac->getMomentum();
				double p_e = p_starac->getEnergy();

				e_s.SetPxPyPzE(e_es[0], e_es[1], e_es[2], e_e);// 4Vektor starca elektrona
				p_s.SetPxPyPzE(p_ps[0], p_ps[1], p_ps[2], p_e); // 4V starca pozitrona
				//e_rec , p_rec su rekontruisani elektron i pozitron

	//			cout << "id: "<< e_starac->id() <<endl;



				TLorentzVector Z1, Z2, higgs;
				Z1 = e_s  - e_rec;
				Z2 = p_s - p_rec;
				//higgs = e_s - e_rec + p_s - p_rec;
				higgs = v2j[0] + v2j[1];
				if (higgs.M() < 0 ) continue;
				if (Z1.E() < 0 || Z2.E()<0) continue;


				m_Z1 = abs(Z1.M());
				m_Z2 = abs(Z2.M());
				E_Z1 = Z1.E();
				E_Z2 = Z2.E();
				m_e1e2 = (e_rec+p_rec).M();
				theta_Z1 = Z1.Theta()*180/M_PI;
				theta_Z2 = Z2.Theta()*180/M_PI;
				theta_e = e_rec.Theta()*180/M_PI;
				theta_p = p_rec.Theta()*180/M_PI;
				p_e1 = e_rec.P();
				p_e2 = p_rec.P();
				pt_e1 = e_rec.Pt();
				pt_e2 = p_rec.Pt();
				E_e1 = e_rec.E();
				E_e2 = p_rec.E();
				pt_e_sistema = (e_s + e_rec).Pt();
				pt_p_sistema = (p_s + p_rec).Pt();
				pt_q1 = v2j[0].Pt();
				pt_q2 = v2j[1].Pt();
				E_q1 = v2j[0].E();
				E_q2 = v2j[1].E();
				theta_H = higgs.Theta();
				E_H = higgs.E();
				pt_H = higgs.Pt();
				m_H = (e_s - e_rec + p_s - p_rec).M();
				m_qq  = (v2j[0] + v2j[1]).M();
				n_pfo = c_n_pfo;
				E_vis_E_H = E_vis - E_H;
				E_miss = 1000 - E_vis;
				logy_12 = m_y12;
				logy_23 = m_y23;

		//		if (e_rec.E() > e_s.E()) cout <<"finalni elektron ima višu energiju!!!\n";
		//		if (p_rec.E() > p_s.E()) cout <<"finalni pozitron ima višu energiju!!!\n";

				TLorentzVector lokalac_Niz_el_in4, lokalac_Niz_poz_in4, lokalac_Niz_el_out4, lokalac_Niz_poz_out4, lokalacZ1, lokalacZ2;

						TVector3 BoostToHiggs = -(higgs.BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona

						lokalac_Niz_el_in4 = e_s;
						lokalac_Niz_el_out4 = e_rec;
						lokalac_Niz_poz_in4 = p_s;
						lokalac_Niz_poz_out4 = p_rec;
						lokalacZ1 = Z1; // 4-vektor Z1-bozona ++++//////--------
						lokalacZ2 = Z2; // 4-vektor Z2-bozona

						//*****BOOST TO HIGGS****
						lokalac_Niz_el_in4.Boost(BoostToHiggs);
						lokalac_Niz_poz_in4.Boost(BoostToHiggs);
						lokalac_Niz_el_out4.Boost(BoostToHiggs);
						lokalac_Niz_poz_out4.Boost(BoostToHiggs);
						lokalacZ1.Boost(BoostToHiggs);
						lokalacZ2.Boost(BoostToHiggs);


						TVector3 Niz_el_in3;
						TVector3 Niz_poz_in3;

						TVector3 Niz_el_out3;
						TVector3 Niz_poz_out3;

						TVector3 OnShell3_Z;
						TVector3 OffShell3_Z;


						Niz_el_in3.SetXYZ(lokalac_Niz_el_in4.X(), lokalac_Niz_el_in4.Y(), lokalac_Niz_el_in4.Z());
						Niz_el_out3.SetXYZ(lokalac_Niz_el_out4.X(), lokalac_Niz_el_out4.Y(), lokalac_Niz_el_out4.Z());
						Niz_poz_in3.SetXYZ(lokalac_Niz_poz_in4.X(), lokalac_Niz_poz_in4.Y(), lokalac_Niz_poz_in4.Z());
						Niz_poz_out3.SetXYZ(lokalac_Niz_poz_out4.X(), lokalac_Niz_poz_out4.Y(), lokalac_Niz_poz_out4.Z());

						OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z()); // impuls Z1
						OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z()); // impuls Z2


						//------ NASA DEFINICIJA UGLA PHI - u ovom slucaju su normale na ravni usmerene isto ------//

						TVector3 brojilac1 = Niz_el_in3.Cross(Niz_el_out3);   // brojilac vektora n1
						Double_t imenilac1 = brojilac1.Mag();   // sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));   // imenilac vektora n1
		//				cout <<"brojilac1.px " <<brojilac1.Px()<< endl;
		//				cout <<"imenilac1 " <<imenilac1<< endl;
						TVector3 n1 = brojilac1 * pow(imenilac1,-1);           // vektor n1

						TVector3 brojilac2 = Niz_poz_in3.Cross(Niz_poz_out3); // brojilac vektora n2
						Double_t imenilac2 = brojilac2.Mag(); // sqrt(pow(brojilac2.X(),2) + pow(brojilac2.Y(),2) + pow(brojilac2.Z(),2));   // imenilac vektora n2

						TVector3 n2 = brojilac2 * pow(imenilac2,-1);           // vektor n2

		//				fi1 = acos(n1.Dot(n2));      // * 180 / M_PI prvi ugao ---> znak plus u argumentu funkcije kada su vektori n1 i n2 u istom smeru i obrnuto
						fi = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));

				leptonTree.Fill();

			//	cout <<"inv M  z1 rec  je : " <<Z1.M()<<endl;
			//	cout <<"inv M  z2  je : " << Z2.M()<<endl;

		//		cout <<"Energija e rec  je : " <<e_rec.E()<<endl;
		//		cout <<"Energija p rec je : " <<p_rec.E()<<endl;

			//	cout <<"inv M Higsovog pozona je  : " <<  (Z1+Z2).M()  <<endl;


				N_za_dva++;
			}//cout <<"Nasao sam event \n";

			//if ( N_KG_elektrona_6 >= 2) cout <<"Broj el. "<< N_KG_elektrona_6 << endl;
			//if ( N_KG_pozitrona_7 >= 2) cout <<"Broj poz. "<< N_KG_pozitrona_7 << endl;
			//cout <<"\\\\\\\\\\\\\\\\\\\\\\\\\\\\n"<<endl;


		} // Kraj WHILE petlje po dogadjajima

/*		cout << "Zatvara se petlja po dogadjajima\n";
		cout << "Broj dogadjaja po fajlu iznosi " << nevent << endl; */

		lcReader -> close();

	} // Kraj petlje po fajlovima

	cout << "Ukupan broj dogadjaja: " << N_ukupnih_dogadjaja << endl;
	cout << "Ukupan broj dogadjaja signala: " << nsignala_mc << endl;
	cout << "Ukupan broj elektrona signala, nelinkovanih: " << N_KG_elektrona << endl;
//	cout << "Ukupan broj pocetnih elektrona 6: " << N_KG_elektrona_6 << endl;
//	cout << "Ukupan broj pocetnih pozitrona 7: " << N_KG_pozitrona_7 << endl;
	//cout<<"broj reko elektrona je : "<< cElec<<endl;
	//cout<<"broj reko pozitrona je : "<< cPosi<<endl;

	cout << "Broj dogadjaja sa dva rekonstruisana leptona sa E > 60 GeV " << N_za_dva << endl;
//	cout << "broj dogadjaja za krivu: " << brDogadjajazaKrivu << endl;
/*cout << "broj dogadjaja sa dva leptona new: " << broj_dogadjaja_2lep << endl;
cout << "broj dogadjaja sa dva leptona new posle Etrack cut: " << n_etrack << endl;
cout << "broj dogadjaja sa dva leptona new posle ratio cut: " << n_ratio << endl;
cout << "broj dogadjaja sa dva leptona new posle IP cut: " << n_IP_brojac << endl;
cout << "broj dogadjaja sa dva leptona new posle kriva cut: " << n_brojac_kriva << endl;
cout << "broj dogadjaja sa dva leptona finalno: " << nlep_final_new << endl;*/







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


	leptonTree.Write();
//	ptTree.Write();
	//pfoTree.Write();

	//elektronskoTree.Write();
	//elektronsko6i7Tree.Write();
	//elektronskoRekonstrTree.Write();

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
