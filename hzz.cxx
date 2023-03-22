//////////////////////////////////// 6.06.2018.
// Ovaj program obradjuje uzorke hzqq_rec/dst_9593_1,2, ... .slcio signala
// radna verzija programa
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
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <UTIL/LCTOOLS.h>
#include <Exceptions.h>
#endif

#include "stdlib.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include "varList.h"
using namespace std;

// ----------------------------------------------------------------------------------------------------------------------

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn/*, const char * rfn*/)
{
#ifdef __CINT__
	gSystem -> Load("${LCIO}/lib/liblcio.so");
	gSystem -> Load("${LCIO}/lib/liblcioDict.so");
#endif

	TFile rootFile("uglovi.root","RECREATE", "uglovi", 1);
	TTree hzz ("hzz", "Generator particle tree");
	varListPairMonitor vl; /* SL specific */
	hzz.Branch("fi", &(vl.fi), "fi");
	hzz.Branch("fi1", &(vl.fi1), "fi1");
	hzz.Branch("psi", &(vl.psi), "psi");
	hzz.Branch("theta1", &(vl.theta1), "theta1");
	hzz.Branch("theta2", &(vl.theta2), "theta2");
	hzz.Branch("mz1", &(vl.mz1), "mz1");
	hzz.Branch("mz2", &(vl.mz2), "mz2");
	hzz.Branch("n1", &(vl.n1), "n1");
	hzz.Branch("n2", &(vl.n2), "n2");
	hzz.Branch("nsc", &(vl.nsc), "nsc");
	hzz.Branch("n1_n2", &(vl.n1_n2), "n1_n2");
	hzz.Branch("n1_nsc", &(vl.n1_nsc), "n1_nsc");
	hzz.Branch("acos_n1_nsc", &(vl.acos_n1_nsc), "acos_n1_nsc");
	hzz.Branch("acos_n1_n2", &(vl.acos_n1_n2), "acos_n1_n2");






	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

	TH1F histofi("histofi", "fi; #Phi", 100, -3.14, 3.14);
	TH1F histofi1 ("histofi1", "fi1; #Phi_{1}", 100, -3.14, 3.14);
	TH1F histocostheta1 ("histocostheta1", "costheta1; cos#theta_{1}", 100, -1, 1);
	TH1F histocostheta2 ("histocostheta2", "costheta2; cos#theta_{2}", 100, -1, 1);
	TH1F histoPsi ("histoPsi", "Psi; Psi", 100, -3.14, 3.14);
	TH1F histoMass_OnShell_Z ("histoMass_OnShell_Z", "Mass_OnShell_Z; m_{Z} (GeV)", 210, 40, 110);
	TH1F histoMass_OffShell_Z ("histoMass_OffShell_Z", "Mass_OffShell_Z; m_{Z*} (GeV)", 210, 0, 70);
	TH1F cosZ ("cosZ", "cosZ; cos#theta*", 100, -1, 1);

	int brojac = 0;


	// Petlja koja iscitava .slcio fajlove
		for(UInt_t iJob = nFirstJob; iJob <= nLastJob; iJob++)
		{
		   cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);
		   try
		   {
			  lcReader -> open(Form("%s%i.slcio", fName.Data(), iJob));
		   }
		   catch(lcio :: IOException &ex)
		   {
			  cout << ". Could not open.\n";
			  continue;
		   }

		cout << ". Reading.\n";

		int broj_Dogadjaja = lcReader -> getNumberOfEvents();   // Ukupan broj događaja

		cout << "Broj dogadjaja po fajlu je: " << broj_Dogadjaja << endl;

		Int_t Nov_dogadjaj = 0;

		// Petlja po dogadjajima
		EVENT::LCEvent* evt = 0;

		while((evt = lcReader -> readNextEvent()) != 0 )    /* && ndogadjaja < 10 */
		{
			Nov_dogadjaj++;
			Int_t Broj_Z_bozona = 0;


			vector<TLorentzVector> Z_niz;                              // prikupljaju se Z bozoni

			vector<TLorentzVector> Niz_2,Niz_1, Niz_H;                // deklarisu se nizovi 4vektora Niz_1, Niz_2 i Niz_H

			TLorentzVector Niz_1p, Niz_1n, Niz_2p, Niz_2n,Z1, Z2, OnShell4_Z, OffShell4_Z, Higs4_H;   // deklarisu se 4vektori Z1, Z2, H

			std::vector<std::string> colNames = *evt -> getCollectionNames();
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");
			//Int_t nz = 0;

			TLorentzVector higgs; // 4vektor u koji se smestaju informacije o Higsovom bozonu i koji se koristi za boost
			TLorentzVector firstHiggs;
    	//	vector <mcParticles>testing;

			bool kvarkovi = false;
			bool leptoni = false;
		    int b_lepton ;
		    int b_kvark ;


		// Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
		for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i); // uzimamo MCCestice
				const EVENT::MCParticleVec & parent = mcParticle -> getParents();//definisemo roditelje mccestice
				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();//definisemo cerke mc cestice

				TLorentzVector temp;                             // 4vektor u koji se smestaju informacije o svakoj cestici
				TLorentzVector majka;							 //4 vektor majke koji je sluzio kao provera.
				const double *p = mcParticle -> getMomentum();   // impuls cestice
				double e = mcParticle -> getEnergy();            // energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	         // 4vektor mc cestice

				if (mcParticle -> getPDG() == 25 && daughter.size()==2) // Higgs koji je predak Z bozonima
				{
					Niz_H.push_back(temp); // niz koji su napravili mirko i jasna
			//s		higgs = temp;//cetvorovektor higgsa
				    for (Int_t k = 0; (int) k < daughter.size(); k++)//prolazimo kroz cerke higgsa
				    {
					TLorentzVector temp_cerke;                        // 4vektor u koji se skupljaju informacije o svakoj cestici
					const double *p = daughter[k] -> getMomentum();   // impuls cestice
					double e = daughter[k] -> getEnergy();	          // energija cestic
					temp_cerke.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor

					Z_niz.push_back(temp_cerke);//niz u kom cuvako cerke od H. sluzi kao test.
				    }


				}// end Higgs


				if (mcParticle -> getPDG() == 25 && daughter.size()==1)//intermediarni Higsov bozon.
				{
					firstHiggs = temp;
				}


				if (mcParticle -> getPDG() == 23)//petlja po z bozonima
				{
				    Broj_Z_bozona++;		 // ovaj brojač može da se iskoristi za Z bozone
				   // Z_niz.push_back(temp);   // prikupljaju se Z bozoni
				    vector<TLorentzVector> cerke_niz;//prikupljaju se potomci od Z bozona
				    TLorentzVector pozitivan;//cetvorovektor fermiona
				    TLorentzVector negativan;//cetvorovektor antifermiona
				    b_lepton = 0 ;
				    b_kvark = 0;

				    for (Int_t k = 0; (int) k < daughter.size(); k++)
				    {
				    	int pdgdaughters = daughter[k] -> getPDG();
				    	if (pdgdaughters < 17 && pdgdaughters > - 17 )//&& pdgdaughters != 12 && pdgdaughters !=14 && pdgdaughters !=16 && pdgdaughters != -12 && pdgdaughters !=-14 && pdgdaughters !=-16)
				    	{
							TLorentzVector temp_cerke;                        // 4vektor u koji se skupljaju informacije o svakoj cestici
							const double *p = daughter[k] -> getMomentum();   // impuls cestice
							double e = daughter[k] -> getEnergy();	          // energija cestice
							temp_cerke.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor

				    		cerke_niz.push_back(temp_cerke);
				    		if (pdgdaughters > 0) pozitivan=temp_cerke;

				    		if (pdgdaughters < 0) negativan=temp_cerke;

				    		if(fabs(pdgdaughters) == 11 ||  fabs(pdgdaughters) ==  13) b_lepton++;

				    		if(fabs(pdgdaughters) < 7 ) b_kvark++;
							// cout << "naelektrisanje "<<daughter[k]->getCharge()<<endl;
			//	    		testing.push_back(daughter[k]);
							const double *p1 = parent[0] -> getMomentum();   // impuls cestice
							double e1 = parent[0] -> getEnergy();	          // energija cestice
							majka.SetPxPyPzE(p1[0], p1[1], p1[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor
				    	}//skupljamo fermione i antifermione.
				    }

					if(Broj_Z_bozona == 1)   // ako je brojac Z bozona = 1 radi se o podacima za Z1 bozon
					{
						 Z1 = temp;
						 Niz_1p = pozitivan;//cetvorovektor fermiona od prvog Z bozona
						 Niz_1n = negativan;//cetvorovektor antifermion od prvog z bozona
						 Niz_1 = cerke_niz;// vektor koji sadrži fermion i antifermion
					//	 cout << "Z1 energija: "<<Z1.E()<< "pdg roditelja: "<<parent[0]->getPDG()<< ", invarijantna masa roditelja: "<< majka.M()<<", invarijantna masa higgsa: "<<higgs.M() <<endl;
					}

					if(Broj_Z_bozona == 2)   // ako je brojac Z bozona = 2 radi se o podacima za Z2 bozon
					{
						Z2 = temp;
						 Niz_2p = pozitivan;
						 Niz_2n = negativan;
					 	 Niz_2 = cerke_niz;
					//	 cout << "Z2 energija: "<<Z2.E()<< " ,pdg roditelja: "<<parent[0]->getPDG()<< ", invarijantna masa roditelja: "<< majka.M()<<", invarijantna masa higgsa: "<<higgs.M()<<"Z1 + Z2 : "<< (Z1+Z2).M()<<endl;

					}
				}//kraj petlje po Z bozonima

		}   // kraj petlje po cesticama

		//cout <<"firstHiggs M: "<<firstHiggs.Pt()<<", higgs.M : "<<higgs.Pt()<<endl;
		//cout<<testing.size();
	//	cout<<"test "<<(Z1+Z2).M()<<endl;
//	if ((Niz_1p+ Niz_1n+ Niz_2p+Niz_2n).M() > 100)	cout << "test2 "<< (Niz_1p+ Niz_1n+ Niz_2p+Niz_2n).M()<<endl;





	if (Niz_1.size() == 2 && Niz_2.size() == 2)   // proverava se da li je broj cestica u oba niza = 2
	{
		//cout<<"u nizu sam"<<endl;
		brojac++;
		TLorentzVector OnShell4_l1;
		TLorentzVector OnShell4_l2;

		TLorentzVector OffShell4_l1;
		TLorentzVector OffShell4_l2;

		TLorentzVector on;
		TLorentzVector off;

		Double_t Mass_OnShell_Z = 0;
		Double_t Mass_OffShell_Z = 0;


		if (b_lepton ==2 && b_kvark==2 )cout <<"broj leptona: "<<b_lepton<<", broj kvarkova: "<<b_kvark<<endl;



		if((Niz_1p + Niz_1n).M() > (Niz_2p + Niz_2n).M())    // if Z1.M() > Z2.M() onda su čestice iz Niz_1 OnShell, a cestice iz Niz_2 su OffShell
		{
		    OnShell4_l1 =  Niz_1p;
		    OnShell4_l2 =  Niz_1n;

		    OffShell4_l1 = Niz_2p;
		    OffShell4_l2 = Niz_2n;

		}

		if((Niz_1p + Niz_1n).M() <= (Niz_2p + Niz_2n).M())   // if Z1.M() <= Z2.M() onda su čestice iz Niz_1 OffShell, a cestice iz Niz_2 su OnShell
		{
		    OnShell4_l1 = Niz_2p;
		    OnShell4_l2 = Niz_2n;

		    OffShell4_l1 = Niz_1p;
		    OffShell4_l2 = Niz_1n;
		}

		on = OnShell4_l1 + OnShell4_l2;
		off= OffShell4_l1 + OffShell4_l2;
		higgs = on + off ;
		Mass_OnShell_Z = on.M();
		Mass_OffShell_Z = off.M();



//		TLorentzVector onshellZ;
	//	onshellZ.SetPxPyPzE(OnShell4_l1.Px()+OnShell4_l2.Px(), OnShell4_l1.Py()+ OnShell4_l2.Py(), OnShell4_l1.Pz()+OnShell4_l2.Pz(),OnShell4_l1.E()+ OnShell4_l2.E());
		vl.mz1 = (OnShell4_l1 + OnShell4_l2).M();
		vl.mz2 = (OffShell4_l1+ OffShell4_l2).M();
		//cout << "test2 "<< (Niz_1p+ Niz_1n+ Niz_2p+Niz_2n).M()<<", higgs mass: "<<higgs.M()<<endl;

	/*	cout << "pozitivne + negativne "<< (Niz_1p+ Niz_1n+ Niz_2p+Niz_2n).M()<<endl;
		cout<<"onshell+ offshell: "<<(OnShell4_l1 + OnShell4_l2+ OffShell4_l1 + OffShell4_l2).M() <<endl;

		cout << "onshell preko 4v3ktora "<< (OnShell4_l1 + OnShell4_l2).M()<<"onshell preko niz1p: " <<", onshell: "<<on.M()<<endl;
		cout << "offshell preko 4v3ktora "<< (OffShell4_l1 + OffShell4_l2).M()<<", offshell: "<<off.M()<<endl;
		cout <<"onshellZ =" << onshellZ.M()<<endl;

		cout <<"zi z2 inv masa: "<<(on+off).M()<<endl;
		cout<<"_________________________________________________________"<<endl;*/


		TVector3 BoostToHiggs = -(higgs.BoostVector());     // prelazi se u koordinatni sistem Higsovog bozona
		TVector3 BoostToOnShellZ = -(on.BoostVector());     // prelazi se u koordinatni sistem OnShell Z bozona
	    TVector3 BoostToOffShellZ = -(off.BoostVector());   // prelazi se u koordinatni sistem OffShell Z bozona

	    TVector3 OnShell3_l1;   // lepton 1 od Z
	 	TVector3 OnShell3_l2;   // lepton 2 od Z

	    TVector3 OffShell3_l1;   // lepton 1 od Z*
	    TVector3 OffShell3_l2;   // lepton 2 od Z*

	    TVector3 nz;

	    TVector3 OnShell3_Z;   // ovo je Z
	  	TVector3 OffShell3_Z;  // ovo je Z*

//**********************************************************************************************************************
	  	//Boost u Higgs
				// 4vektori leptona su bustovani u sistem Higsovog bozona, potom su zadati 3vektori leptona i na kraju su izracunati uglovi fi, fi1 i psi

				TLorentzVector lokalac_OnShell4_l1, lokalac_OnShell4_l2, lokalac_OffShell4_l1, lokalac_OffShell4_l2, lokalacZ1, lokalacZ2;

				lokalac_OnShell4_l1 = OnShell4_l1;
				lokalac_OnShell4_l2 = OnShell4_l2;
				lokalac_OffShell4_l1 = OffShell4_l1;
				lokalac_OffShell4_l2 = OffShell4_l2;
				lokalacZ1 = on;
				lokalacZ2 = off;
			//	cout <<"brojilac1: "<<Niz_2n.X()<<endl;




			//	cout <<"z1 + z2 impuls: "<< (lokalacZ1 + lokalacZ2).M()<<endl;

				lokalac_OnShell4_l1.Boost(BoostToHiggs);
				lokalac_OnShell4_l2.Boost(BoostToHiggs);
				lokalac_OffShell4_l1.Boost(BoostToHiggs);
				lokalac_OffShell4_l2.Boost(BoostToHiggs);
				lokalacZ1.Boost(BoostToHiggs);
				lokalacZ2.Boost(BoostToHiggs);
				higgs.Boost(BoostToHiggs);
		//		cout<<"lokalac_OnShell4_l1 posle boosta"<<lokalac_OnShell4_l1.Pt()<<endl;
		//		cout <<"z1 + z2 impuls: "<< (higgs.P())<<endl;
			//	cout <<"z1 + z2 impuls: "<< (lokalacZ1 + lokalacZ2).P()<<endl;

				double thetaZ = lokalacZ1.Theta();
				cosZ.Fill(TMath::Cos(thetaZ));//cos*theta Radi lepo

				//trovektor impulsa
				OnShell3_l1.SetXYZ(lokalac_OnShell4_l1.X(), lokalac_OnShell4_l1.Y(), lokalac_OnShell4_l1.Z());
				OnShell3_l2.SetXYZ(lokalac_OnShell4_l2.X(), lokalac_OnShell4_l2.Y(), lokalac_OnShell4_l2.Z());

				OffShell3_l1.SetXYZ(lokalac_OffShell4_l1.X(), lokalac_OffShell4_l1.Y(), lokalac_OffShell4_l1.Z());
				OffShell3_l2.SetXYZ(lokalac_OffShell4_l2.X(), lokalac_OffShell4_l2.Y(), lokalac_OffShell4_l2.Z());

				nz.SetXYZ(0,0,1);   // z-osa

				OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z()); // impuls OnShell Z
				OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z()); // impuls OnShell Z

	    	//    cout <<"OnShell3_Z: "<<OnShell3_Z.Mag() - OffShell3_Z.Mag()<<endl;

				TVector3 brojilac1 = OnShell3_l1.Cross(OnShell3_l2);   // brojilac vektora n1
				Double_t imenilac1 = brojilac1.Mag();//sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));                  // imenilac vektora n1
				TVector3 unit_n1 = brojilac1.Unit();

	     		TVector3 n1 = brojilac1 * pow(imenilac1,-1);           // vektor n1

				TVector3 brojilac2 = OffShell3_l1.Cross(OffShell3_l2); // brojilac vektora n2
		    	Double_t imenilac2 = brojilac2.Mag();//sqrt(pow(brojilac2.X(),2) + pow(brojilac2.Y(),2) + pow(brojilac2.Z(),2)); // brojilac2.Mag();                  // imenilac vektora n1
				TVector3 unit_n2 = brojilac2.Unit();

				TVector3 n2 = brojilac2 * pow(imenilac2,-1);           // vektor n2

				TVector3 brojilac3 = nz.Cross(OnShell3_Z);             // brojilac vektora n_sc
				Double_t imenilac3 = brojilac3.Mag();//sqrt(pow(brojilac3.X(),2) + pow(brojilac3.Y(),2) + pow(brojilac3.Z(),2));  // brojilac3.Mag();                  // brojilac vektora n_sc

				TVector3 nsc = brojilac3 * pow(imenilac3,-1);          // vektor n_sc

				double fiBrojilac = OnShell3_Z.Dot(unit_n1.Cross(unit_n2));
				double arcuscosinusN1 = acos(-unit_n1.Dot(unit_n2));
				double vrednost = -unit_n1.Dot(unit_n2);

			//	double dotovi = -n1.Dot(n2);
			//	cout <<"dotovi: "<<dotovi<<endl;

				//double x1 =  TMath::ACos(dotovi);
				//double x2 = 1/TMath::Cos(dotovi);

			//		cout<< "x1= "<<x1 <<", x2 = "<<x2<<endl;
			//	double ver1 = 0;
				//double ver2 = 0;

				 //ver1 = OnShell3_Z.Dot(n1.Cross(n2)) * pow(abs(OnShell3_Z.Dot(n1.Cross(n2))),-1);
				 //ver2 = fiBrojilac / abs(fiBrojilac);

		//		cout<<"vrednost: "<< (fiBrojilac/abs(fiBrojilac))*vrednost<<endl;
	//			cout<<"TMATH::ACOS: "<<TMath::ACos(vrednost)<<endl;

	//	if (ver1 != ver2)		cout << " verzija 1: "<<ver1<<", verzija2: "<<ver2<<endl;

			    Double_t fi = (fiBrojilac/ abs(fiBrojilac)) *TMath::ACos(vrednost); //ver2 * x1;//( OnShell3_Z.Dot(n1.Cross(n2)) * pow(abs(OnShell3_Z.Dot(n1.Cross(n2))),-1) )* (acos(-n1.Dot(n2)));      // * 180 / M_PI prvi ugao, Fi

				Double_t fi1 = (OnShell3_Z.Dot(n1.Cross(nsc)) * pow(abs(OnShell3_Z.Dot(n1.Cross(nsc))),-1)) * acos(n1.Dot(nsc));   // * 180 / M_PI drugi ugao, Fi1

				Double_t Psi = fi1 + fi/2;   // peti ugao, Psi

			//	cout <<"vrednost : "<<vrednost<<endl;

	//			cout <<"fi : "<<fi<<endl;

			//	cout << "redni broj dogadjaja je : "<<Nov_dogadjaj<<endl;
//
			//	cout<<"_______________________________________________"<<endl;

				histofi.Fill(fi);
				histofi1.Fill(fi1);
				histoPsi.Fill(Psi);
				vl.fi = fi;
				vl.fi1 = fi1;


				//cout<<"acos vrednosi: " <<TMath::ACos(-unit_n1.Dot(unit_n2))<<endl;
			//	cout<<"fi : "<<(-n1.Dot(n2))<<", arccos : "<< acos(-n1.Dot(n2))<<endl;

			//	cout<<"fi1 : "<<(n1.Dot(nsc))<<", arccos : "<< acos(n1.Dot(nsc))<<endl;
//
				histoMass_OnShell_Z.Fill(Mass_OnShell_Z);
				histoMass_OffShell_Z.Fill(Mass_OffShell_Z);




				// 4vektori leptona su bustovani u sistem Higsa + OnShell Z bozona, potom su zadati 3vektori leptona i na kraju je izracunat kosinus ugla teta1

				TLorentzVector lokal_OnShell4_l1, lokal_OnShell4_l2, lokal_OffShell4_l1, lokal_OffShell4_l2;

				lokal_OnShell4_l1 = OnShell4_l1;
				lokal_OnShell4_l2 = OnShell4_l2;
				lokal_OffShell4_l1 = OffShell4_l1;
				lokal_OffShell4_l2 = OffShell4_l2;



				lokal_OnShell4_l1.Boost(BoostToOnShellZ);
				lokal_OnShell4_l2.Boost(BoostToOnShellZ);
				lokal_OffShell4_l1.Boost(BoostToOnShellZ);
				lokal_OffShell4_l2.Boost(BoostToOnShellZ);

				OnShell3_l1.SetXYZ(lokal_OnShell4_l1.X(), lokal_OnShell4_l1.Y(), lokal_OnShell4_l1.Z());
				OnShell3_l2.SetXYZ(lokal_OnShell4_l2.X(), lokal_OnShell4_l2.Y(), lokal_OnShell4_l2.Z());

				OffShell3_l1.SetXYZ(lokal_OffShell4_l1.X(), lokal_OffShell4_l1.Y(), lokal_OffShell4_l1.Z());
				OffShell3_l2.SetXYZ(lokal_OffShell4_l2.X(), lokal_OffShell4_l2.Y(), lokal_OffShell4_l2.Z());

				OnShell3_Z = OnShell3_l1 + OnShell3_l2;
				OffShell3_Z = OffShell3_l1 + OffShell3_l2;

			//	cout << "OffShell3_l1 P :"<< OffShell3_l1.Mag()<<endl;

				Double_t theta1 = -OffShell3_Z.Dot(OnShell3_l1) / (sqrt(pow(OffShell3_Z.X(),2) + pow(OffShell3_Z.Y(),2) + pow(OffShell3_Z.Z(),2)) * sqrt(pow(OnShell3_l1.X(),2) + pow(OnShell3_l1.Y(),2) + pow(OnShell3_l1.Z(),2)));   // kosinus treceg ugla Teta1

				histocostheta1.Fill(theta1);
				vl.theta1 = theta1;



				// 4vektori leptona su bustovani u sistem Higsa + u sistem OffShell Z bozona, potom su zadati 3vektori leptona i na kraju je izracunat kosinus ugla teta2

				TLorentzVector lok_OnShell4_l1, lok_OnShell4_l2, lok_OffShell4_l1, lok_OffShell4_l2, lokZ1, lokZ2;

				lok_OnShell4_l1 = OnShell4_l1;
				lok_OnShell4_l2 = OnShell4_l2;
				lok_OffShell4_l1 = OffShell4_l1;
				lok_OffShell4_l2 = OffShell4_l2;
				lokZ1 = on;
				lokZ2= off;

				lok_OnShell4_l1.Boost(BoostToOffShellZ);
				lok_OnShell4_l2.Boost(BoostToOffShellZ);
				lok_OffShell4_l1.Boost(BoostToOffShellZ);
				lok_OffShell4_l2.Boost(BoostToOffShellZ);
				lokZ1.Boost(BoostToOffShellZ);
				lokZ2.Boost(BoostToOffShellZ);

				OnShell3_l1.SetXYZ(lok_OnShell4_l1.X(), lok_OnShell4_l1.Y(), lok_OnShell4_l1.Z());
				OnShell3_l2.SetXYZ(lok_OnShell4_l2.X(), lok_OnShell4_l2.Y(), lok_OnShell4_l2.Z());

				OffShell3_l1.SetXYZ(lok_OffShell4_l1.X(), lok_OffShell4_l1.Y(), lok_OffShell4_l1.Z());
				OffShell3_l2.SetXYZ(lok_OffShell4_l2.X(), lok_OffShell4_l2.Y(), lok_OffShell4_l2.Z());

				OnShell3_Z = OnShell3_l1 + OnShell3_l2;
				OffShell3_Z = OffShell3_l1 + OffShell3_l2;

			//	cout<<"OnShell3_l2.Px = "<< OnShell3_l2.Px()<<",lok_OnShell4_l2.Px:  "<<lok_OnShell4_l2.Px()<<endl;//slaže se za px
			//	cout <<"p Z2= "<<lokZ2.P()<<endl;

		//		cout <<"p lokZ1: " <<lokZ1.P()<<", p OnshellZ: "<<OnShell3_Z.Mag()<<endl;

				Double_t theta2 = (-OnShell3_Z.Dot(OffShell3_l1) / (OnShell3_Z.Mag()*OffShell3_l1.Mag()));//(sqrt(pow(OnShell3_Z.X(),2) + pow(OnShell3_Z.Y(),2) + pow(OnShell3_Z.Z(),2)) * sqrt(pow(OffShell3_l1.X(),2) + pow(OffShell3_l1.Y(),2) + pow(OffShell3_l1.Z(),2)));  // kosinus cetvrtog ugla Teta2

				histocostheta2.Fill(theta2);
				vl.theta2 = theta2;
				vl.n1 = unit_n1.Mag();
				vl.n2 = unit_n2.Mag();
				vl.nsc = nsc.Mag();
				vl.n1_n2 = -unit_n1.Dot(unit_n2);
				vl.n1_nsc = unit_n1.Dot(nsc);
				vl.acos_n1_n2 = TMath::ACos(-unit_n1.Dot(unit_n2))*180/M_PI;
				vl.acos_n1_nsc = TMath::ACos(unit_n1.Dot(nsc))*180/M_PI;



				hzz.Fill();

		}
	} // kraj petlje po dogadjajima

		lcReader -> close();

 } // Kraj petlje po fajlovima

		cout <<"brojac: "<<brojac<<endl;

		histofi.Write();
		histofi1.Write();
		histocostheta1.Write();
		histocostheta2.Write();
		histoPsi.Write();
		histoMass_OnShell_Z.Write();
		histoMass_OffShell_Z.Write();
		cosZ.Write();
		hzz.Write();

		rootFile.Close();

return 0;

} // Kraj slcio2appTree funkcije

// ----------------------------------------------------------------------------------------------------------------------

Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;

	UInt_t nFirstJob = 1;
	if(argc > iarg)
	{
	  nFirstJob = atoi(argv[iarg]);
	  iarg++;
	}

	UInt_t nLastJob = 107;
	if(argc > iarg)
	{
	  nLastJob = atoi(argv[iarg]);
	  iarg++;
	}

	TString fName = "hzqq_dst_9593_";
	if(argc > iarg)
	{
	  fName = argv[iarg];   // iarg++;
	}

	return slcio2appTree(nFirstJob, nLastJob, fName.Data()/*, rfName.Data()*/);
}

