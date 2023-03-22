//////////////////////////////////// 11.06.2018.
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

	TFile rootFile("MirkoUglovi3.root","RECREATE", "MirkoUglovi3", 1);
/*	TTree hzz("hzz", "Generator of the particle tree");
	varListPairMonitor vl;               SL specific
	hzz.Branch("fi", &(vl.fi), "fi");
	hzz.Branch("fi1", &(vl.fi1), "fi1");
	hzz.Branch("psi", &(vl.psi), "psi");
	hzz.Branch("theta1", &(vl.theta1), "theta1");
	hzz.Branch("theta2", &(vl.theta2), "theta2");
	hzz.Branch("mz1", &(vl.mz1), "mz1");
	hzz.Branch("mz2", &(vl.mz2), "mz2"); */

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader();
	TString fName = fn;
	stringstream fNameStream;

	TH1F histoFi("Histogram #it{#phi}", " ; #it{#phi}", 100, -3.14, 3.14);
	TH1F histoFi1 ("Histogram #it{#phi}_{1}", " ; #it{#phi}_{1}", 100, -3.14, 3.14);
	TH1F histoCosTheta1 ("Histogram cos #it{#theta}_{1}", " ; cos #it{#theta}_{1}", 100, -1, 1);
	TH1F histoCosTheta2 ("Histogram cos #it{#theta}_{2}", " ; cos #it{#theta}_{2}", 100, -1, 1);
	TH1F histoPsi ("Histogram #it{#psi}", " ; #it{#psi}", 100, -3.14, 3.14);
	TH1F histoMass_OnShell_Z ("Histogram #it{m}_{Z}", " ; #it{m}_{Z}", 210, 40, 110);
	TH1F histoMass_OffShell_Z ("Histogram #it{m}_{Z*}", " ; #it{m}_{Z*}", 210, 0, 70);
	TH1F histoCosTheta ("Histogram cos #it{#theta}*", " ; cos #it{#theta}*", 100, -1, 1);

	Int_t brojac = 0;

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

		while((evt = lcReader -> readNextEvent()) != 0)   /* && ndogadjaja < 10 */
		{
			Nov_dogadjaj++;
			Int_t Broj_Z_bozona = 0;

			vector<TLorentzVector> Z_niz;                              // prikupljaju se Z bozoni

			vector<TLorentzVector> Niz_1, Niz_2, Niz_H;                // deklarisu se nizovi 4vektora Niz_1, Niz_2 i Niz_H

			TLorentzVector Niz_1p, Niz_1n, Niz_2p, Niz_2n, Z1, Z2, OnShell4_Z, OffShell4_Z, Higs4_H;   // deklarisu se 4vektori Z1, Z2, H

			vector<Int_t> Identifikacioni_broj_cestice;                // deklarise se niz koji cuva celobrojne PDG-ove cestice

			bool kvarkovi = false;                					   // deklarise se logicka promenljiva za par kvarkova
			bool leptoni = false;                					   // deklarise se logicka promenljiva za par leptona

			std::vector<std::string> colNames = *evt -> getCollectionNames();
			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");
			// Int_t nz = 0;

			TLorentzVector higgs;   // 4vektor u koji se smestaju informacije o Higsovom bozonu i koji se koristi za boost
			TLorentzVector firstHiggs;

		   // Petlja preko koje prolazimo kroz sve čestice po svakom dogadjaju
		   for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
		   {
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				const EVENT::MCParticleVec & parent = mcParticle -> getParents();
				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

				TLorentzVector temp;                             // 4vektor u koji se smestaju informacije o svakoj cestici
				TLorentzVector majka;

				//Identifikacioni_broj_cestice[i] = mcParticle -> getPDG();

				const double *p = mcParticle -> getMomentum();   // impuls cestice
				double e = mcParticle -> getEnergy();            // energija čestice
				temp.SetPxPyPzE(p[0], p[1], p[2], e);  	         // zapisuju se vrednosti energije i impulsa u 4vektor

				if (mcParticle -> getPDG() == 25 && daughter.size() == 2)
				{
					Niz_H.push_back(temp);
					higgs = temp;

					for (Int_t k = 0; (int) k < daughter.size(); k++)
					{
						TLorentzVector temp_cerke;                        // 4vektor u koji se skupljaju informacije o svakoj cestici
						const double *p = daughter[k] -> getMomentum();   // impuls cestice
						double e = daughter[k] -> getEnergy();	          // energija cestice
						temp_cerke.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor

						Z_niz.push_back(temp_cerke);
					}
				}   // kraj petlje za Higsov bozon

				if (mcParticle -> getPDG() == 25 && daughter.size() == 1)
				{
					firstHiggs = temp;
				}

				if (mcParticle -> getPDG() == 23)
				{
				    Broj_Z_bozona++;		 // ovaj brojac moze da se iskoristi za Z bozone

				    vector<TLorentzVector> cerke_niz;
				    TLorentzVector pozitivan;
				    TLorentzVector negativan;

				    for (Int_t k = 0; (int) k < daughter.size(); k++)
				    {
				    	Int_t pdgdaughters = daughter[k] -> getPDG();
				    	if (pdgdaughters < 17 && pdgdaughters > -17 && pdgdaughters != 12 && pdgdaughters != 14 && pdgdaughters != 16 && pdgdaughters != -12 && pdgdaughters != -14 && pdgdaughters != -16)
				    	{
							TLorentzVector temp_cerke;                        // 4vektor u koji se skupljaju informacije o svakoj cestici
							const double *p = daughter[k] -> getMomentum();   // impuls cestice
							double e = daughter[k] -> getEnergy();	          // energija cestice
							temp_cerke.SetPxPyPzE(p[0], p[1], p[2], e);  	  // zapisuju se vrednosti energije i impulsa u 4vektor

				    		cerke_niz.push_back(temp_cerke);

				    		if (pdgdaughters > 0) pozitivan = temp_cerke;
				    		if (pdgdaughters < 0) negativan = temp_cerke;

				    		const double *p1 = parent[0] -> getMomentum();    // impuls cestice
				    		double e1 = parent[0] -> getEnergy();	          // energija cestice
				    		majka.SetPxPyPzE(p1[0], p1[1], p1[2], e1);        // zapisuju se vrednosti energije i impulsa u 4vektor
				    	}

				    	if (pdgdaughters == -6 || pdgdaughters == 6 || pdgdaughters == -5 || pdgdaughters == 5 || pdgdaughters == -4 || pdgdaughters == 4 || pdgdaughters == -3 || pdgdaughters == 3 || pdgdaughters == -2 || pdgdaughters == 2 || pdgdaughters == -1 || pdgdaughters == 1)
				    	{
				    		kvarkovi = true;
				    	}

				    	if (pdgdaughters == -11 || pdgdaughters == 11 || pdgdaughters == -13 || pdgdaughters == 13)
				    	{
				    		leptoni = true;
				    	}

				    	if (kvarkovi == true && leptoni == true)
				    	{
				    		if (Broj_Z_bozona == 1)   // ako je brojac Z bozona = 1 radi se o podacima za Z1 bozon
				    		{
				    			Z1 = temp;
				    			Niz_1p = pozitivan;
				    			Niz_1n = negativan;
				    			Niz_1 = cerke_niz;
				    		}

				    		if (Broj_Z_bozona == 2)   // ako je brojac Z bozona = 2 radi se o podacima za Z2 bozon
				    		{
				    			Z2 = temp;
				    			Niz_2p = pozitivan;
				    			Niz_2n = negativan;
				    			Niz_2 = cerke_niz;
				    		}
				    	}
				    }   // end FOR petlje
				}

		   }   // kraj petlje po cesticama

		 //  cout<<"izlazim iz petlje po cesticama"<<endl;
		   if (Niz_1.size() == 2 && Niz_2.size() == 2)   // proverava se da li je broj cestica u oba niza = 2
		   {

			   brojac++;
			   TLorentzVector OnShell4_l1;
			   TLorentzVector OnShell4_l2;

			   TLorentzVector OffShell4_l1;
			   TLorentzVector OffShell4_l2;

			   TLorentzVector on;
			   TLorentzVector off;

			   Double_t Mass_OnShell_Z = 0;
			   Double_t Mass_OffShell_Z = 0;

			   if(Z1.M() > Z2.M())    // if Z1.M() > Z2.M() onda su čestice iz Niz_1 OnShell, a cestice iz Niz_2 su OffShell
			   {
				   OnShell4_l1 =  Niz_1p;
				   OnShell4_l2 =  Niz_1n;

				   OffShell4_l1 = Niz_2p;
				   OffShell4_l2 = Niz_2n;
			   }

			   if(Z1.M() <= Z2.M())   // if Z1.M() <= Z2.M() onda su čestice iz Niz_1 OffShell, a cestice iz Niz_2 su OnShell
			   {
				   OnShell4_l1 = Niz_2p;
				   OnShell4_l2 = Niz_2n;

				   OffShell4_l1 = Niz_1p;
				   OffShell4_l2 = Niz_1n;
				}

				on = OnShell4_l1 + OnShell4_l2;
				off = OffShell4_l1 + OffShell4_l2;

				Mass_OnShell_Z = on.M();
				Mass_OffShell_Z = off.M();

				TLorentzVector onshellZ;
				onshellZ.SetPxPyPzE(OnShell4_l1.Px() + OnShell4_l2.Px(), OnShell4_l1.Py() + OnShell4_l2.Py(), OnShell4_l1.Pz() + OnShell4_l2.Pz(),OnShell4_l1.E() + OnShell4_l2.E());
//				vl.mz1 = (OnShell4_l1 + OnShell4_l2).M();

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

			   // 4vektori leptona su bustovani u sistem Higsovog bozona, potom su zadati 3vektori leptona i na kraju su izracunati uglovi fi, fi1 i psi

			   TLorentzVector lokalac_OnShell4_l1, lokalac_OnShell4_l2, lokalac_OffShell4_l1, lokalac_OffShell4_l2, lokalacZ1, lokalacZ2;

			   lokalac_OnShell4_l1 = OnShell4_l1;
			   lokalac_OnShell4_l2 = OnShell4_l2;
			   lokalac_OffShell4_l1 = OffShell4_l1;
			   lokalac_OffShell4_l2 = OffShell4_l2;
			   lokalacZ1 = on;
			   lokalacZ2 = off;

			   lokalac_OnShell4_l1.Boost(BoostToHiggs);
			   lokalac_OnShell4_l2.Boost(BoostToHiggs);
			   lokalac_OffShell4_l1.Boost(BoostToHiggs);
			   lokalac_OffShell4_l2.Boost(BoostToHiggs);
			   lokalacZ1.Boost(BoostToHiggs);
			   lokalacZ2.Boost(BoostToHiggs);
			   higgs.Boost(BoostToHiggs);

			   double thetaZ = lokalacZ1.Theta();
			   histoCosTheta.Fill(TMath::Cos(thetaZ));

			   OnShell3_l1.SetXYZ(lokalac_OnShell4_l1.X(), lokalac_OnShell4_l1.Y(), lokalac_OnShell4_l1.Z());
			   OnShell3_l2.SetXYZ(lokalac_OnShell4_l2.X(), lokalac_OnShell4_l2.Y(), lokalac_OnShell4_l2.Z());

			   OffShell3_l1.SetXYZ(lokalac_OffShell4_l1.X(), lokalac_OffShell4_l1.Y(), lokalac_OffShell4_l1.Z());
			   OffShell3_l2.SetXYZ(lokalac_OffShell4_l2.X(), lokalac_OffShell4_l2.Y(), lokalac_OffShell4_l2.Z());

			   nz.SetXYZ(0,0,1);   // z-osa

			   OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z()); // impuls OnShell Z
			   OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z()); // impuls OffShell Z

			   TVector3 brojilac1 = OnShell3_l1.Cross(OnShell3_l2);   // brojilac vektora n1
			   Double_t imenilac1 = brojilac1.Mag();   // sqrt(pow(brojilac1.X(),2) + pow(brojilac1.Y(),2) + pow(brojilac1.Z(),2));   // imenilac vektora n1

			   TVector3 n1 = brojilac1 * pow(imenilac1,-1);           // vektor n1

			   TVector3 brojilac2 = OffShell3_l1.Cross(OffShell3_l2); // brojilac vektora n2
			   Double_t imenilac2 = brojilac2.Mag(); // sqrt(pow(brojilac2.X(),2) + pow(brojilac2.Y(),2) + pow(brojilac2.Z(),2));   // imenilac vektora n1

			   TVector3 n2 = brojilac2 * pow(imenilac2,-1);           // vektor n2

			   TVector3 brojilac3 = nz.Cross(OnShell3_Z);             // brojilac vektora n_sc
			   Double_t imenilac3 = brojilac3.Mag();    // sqrt(pow(brojilac3.X(),2) + pow(brojilac3.Y(),2) + pow(brojilac3.Z(),2));   // imenilac vektora n_sc

			   TVector3 nsc = brojilac3 * pow(imenilac3,-1);          // vektor n_sc

			   Double_t fi = OnShell3_Z.Dot(n1.Cross(n2)) * pow(abs(OnShell3_Z.Dot(n1.Cross(n2))),-1) * acos(-n1.Dot(n2));      // * 180 / M_PI prvi ugao, Fi

			   Double_t fi1 = OnShell3_Z.Dot(n1.Cross(nsc)) * pow(abs(OnShell3_Z.Dot(n1.Cross(nsc))),-1) * acos(n1.Dot(nsc));   // * 180 / M_PI drugi ugao, Fi1

			   Double_t Psi = fi1 + fi/2;   // peti ugao, Psi

			   histoFi.Fill(fi);
			   histoFi1.Fill(fi1);
			   histoPsi.Fill(Psi);

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

			   Double_t theta1 = -OffShell3_Z.Dot(OnShell3_l1) / (sqrt(pow(OffShell3_Z.X(),2) + pow(OffShell3_Z.Y(),2) + pow(OffShell3_Z.Z(),2)) * sqrt(pow(OnShell3_l1.X(),2) + pow(OnShell3_l1.Y(),2) + pow(OnShell3_l1.Z(),2)));   // kosinus treceg ugla Teta1

			   histoCosTheta1.Fill(theta1);

			   // 4vektori leptona su bustovani u sistem Higsa + u sistem OffShell Z bozona, potom su zadati 3vektori leptona i na kraju je izracunat kosinus ugla teta2

			   TLorentzVector lok_OnShell4_l1, lok_OnShell4_l2, lok_OffShell4_l1, lok_OffShell4_l2,  lokZ1, lokZ2;

			   lok_OnShell4_l1 = OnShell4_l1;
			   lok_OnShell4_l2 = OnShell4_l2;
			   lok_OffShell4_l1 = OffShell4_l1;
			   lok_OffShell4_l2 = OffShell4_l2;
			   lokZ1 = on;
			   lokZ2 = off;

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

			   Double_t theta2 = -OnShell3_Z.Dot(OffShell3_l1) / (sqrt(pow(OnShell3_Z.X(),2) + pow(OnShell3_Z.Y(),2) + pow(OnShell3_Z.Z(),2)) * sqrt(pow(OffShell3_l1.X(),2) + pow(OffShell3_l1.Y(),2) + pow(OffShell3_l1.Z(),2)));   // kosinus cetvrtog ugla Teta2

			   histoCosTheta2.Fill(theta2);

//			   hzz.Fill();
		   }
		} // kraj petlje po dogadjajima

		lcReader -> close();

	} // Kraj petlje po fajlovima

	cout << "Brojac: "<< brojac << endl;

	histoFi.Write();
	histoFi1.Write();
	histoCosTheta1.Write();
	histoCosTheta2.Write();
	histoPsi.Write();
	histoMass_OnShell_Z.Write();
	histoMass_OffShell_Z.Write();
	histoCosTheta.Write();
//	hzz.Write();

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
