// -------------------------
// Autor: NATASA
// Izmuljali: MIRKO i JASNA
// Datum: 2021.12.3.
// -------------------------
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
	#include "lcio.h"
	#include <IOIMPL/LCFactory.h>
	#include <EVENT/LCCollection.h>
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

//#ifndef ROOT_Math_GenVector_Translation3D
#define ROOT_Math_GenVector_Translation3D  1
#include "Math/GenVector/DisplacementVector3D.h"


#include "Math/GenVector/Plane3D.h"

#include "Math/GenVector/PositionVector3Dfwd.h"

#include "Math/GenVector/LorentzVectorfwd.h"

#include <iostream>
#include <type_traits>

using namespace ROOT;

using namespace Math;

using namespace Impl;


using namespace std;

// ----------------------------------------------------------------------------------------------------------------------
Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
	#ifdef __CINT__
		gSystem -> Load("${LCIO}/lib/liblcio.so");
		gSystem -> Load("${LCIO}/lib/liblcioDict.so");
	#endif

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;
	TTree leptonTree ("leptonTree", "Generator particle tree");

	Float_t fi, fi2;
	Float_t fiOgawa;
	Float_t m_Z1, m_Z2;
	Float_t theta1, theta2;
	Float_t mH;

	leptonTree.Branch("fi", &fi, "fi");
	leptonTree.Branch("fiOgawa", &fiOgawa, "fiOgawa");
	leptonTree.Branch("theta1", &theta1, "theta1");
	leptonTree.Branch("theta2", &theta2, "theta2");
	leptonTree.Branch("mH", &mH, "mH");
	leptonTree.Branch("m_z1", &m_Z1, "m_z1");
	leptonTree.Branch("m_z2", &m_Z2, "m_z2");

	Int_t brojac = 0;
	Int_t Nov_dogadjaj = 0;
	Int_t nlep1 = 0;
	Int_t nlep2 = 0;

	// Petlja koja iscitava .slcio fajlove
	for(UInt_t iJob = nFirstJob; iJob <= nLastJob; iJob++)
	{
	   cout << "Otvara se " << Form("%s%i.slcio", fName.Data(), iJob);
	   try
	   {
		  lcReader -> open(Form("%s%i.slcio", fName.Data(), iJob));
	   }
	   catch(lcio :: IOException &ex)
	   {
		  cout << ". Ne mere se otvoriti.\n";
		  continue;
	   }

		cout << ". Ucitavanje.\n";

		int broj_Dogadjaja = lcReader -> getNumberOfEvents();   // Ukupan broj dogadjaja

		cout << "Broj dogadjaja po fajlu je: " << broj_Dogadjaja << endl;

		// Petlja po dogadjajima
		EVENT::LCEvent* evt = 0;

		while((evt = lcReader -> readNextEvent()) != 0)   /* && ndogadjaja < 10 */
		{
			Nov_dogadjaj++;
			Int_t Broj_H_bozona = 0;

			TLorentzVector Niz_el_in, Niz_el_in2, Niz_poz_in, Niz_el_out, Niz_poz_out, Higs4_H, Z1_boson, Z2_boson, Niz_elpoz_in;   // deklarisu se 4vektori elektrona, pozitrona, Higsa

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");

			TLorentzVector higgs;   // 4vektor u koji se smestaju informacije o Higsovom bozonu i koji se koristi za boost

			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);
				const EVENT::MCParticleVec & parent = mcParticle -> getParents();
				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

				if ( mcParticle -> getPDG() == 11 && mcParticle -> getParents().size() == 0)  //upadni elektron
				{
					for (Int_t k = 0; k < daughter.size(); k++)
					{
						Int_t pdgdaughters = daughter[k] -> getPDG();
						TLorentzVector temp_cerke1;
						const double *p = daughter[k] -> getMomentum();
						double e = daughter[k] -> getEnergy();
						temp_cerke1.SetPxPyPzE(p[0], p[1], p[2], e);

						if (pdgdaughters == 11)
						{
							Niz_el_in = temp_cerke1;

							const EVENT::MCParticleVec & granddaughter = daughter[k]->getDaughters();

							for (Int_t h = 0; h < granddaughter.size(); h++)
							{
								Int_t pdggranddaughters = granddaughter[h] -> getPDG();
								TLorentzVector temp_unuke1;
								const double *p = granddaughter[h] -> getMomentum();
								double e = granddaughter[h] -> getEnergy();
								temp_unuke1.SetPxPyPzE(p[0], p[1], p[2], e);
								if (pdggranddaughters == 25) { higgs = temp_unuke1; }

								if (pdggranddaughters == 11) { Niz_el_out = temp_unuke1; }

								if (pdggranddaughters == -11) { Niz_poz_out = temp_unuke1; }

							} // FOR po unukama

						}
						if (pdgdaughters == -11)
						{ Niz_poz_in = temp_cerke1;	}

					} // FOR za inicijalni e-e+
				}
			} // FOR po MC

			TLorentzVector Z1;
			TLorentzVector Z2;

			m_Z1 = 0;
			m_Z2 = 0;

			Z1 = Niz_el_in - Niz_el_out;
			Z2 = Niz_poz_in - Niz_poz_out;

			m_Z1 = Z1.M();
			m_Z2 = Z2.M();

			TLorentzVector Niz_el_in4, Niz_poz_in4, Niz_el_out4, Niz_poz_out4, Higs4_H4;

			TVector3 BoostToHiggs = -(higgs.BoostVector());    // prelazi se u koordinatni sistem Higsovog bozona
			TVector3 BoostToOnShellZ = -(Z1.BoostVector());    // prelazi se u koordinatni sistem Z1 bozona
			TVector3 BoostToOffShellZ = -(Z2.BoostVector());   // prelazi se u koordinatni sistem Z2 bozona

			TVector3 Niz_el_in3;
			TVector3 Niz_poz_in3;

			TVector3 Niz_el_out3;
			TVector3 Niz_poz_out3;

			TVector3 OnShell3_Z;   // ovo je Z1
			TVector3 OffShell3_Z;  // ovo je Z2

			TVector3 higgsl3;

			TVector3 nz;

			TLorentzVector lokalac_Niz_el_in4, lokalac_Niz_poz_in4, lokalac_Niz_el_out4, lokalac_Niz_poz_out4, lokalacZ1, lokalacZ2, lokalac_higgs;
			lokalac_Niz_el_in4 = Niz_el_in;
			lokalac_Niz_el_out4 = Niz_el_out;
		  lokalac_Niz_poz_in4 = Niz_poz_in;
			lokalac_Niz_poz_out4 = Niz_poz_out;
			lokalacZ1 = Niz_el_in - Niz_el_out;		// 4-vektor Z1-bozona
			lokalacZ2 = Niz_poz_in - Niz_poz_out;	// 4-vektor Z2-bozona
			lokalac_higgs = higgs;

			//*** BOOST TO HIGGS ***
			lokalac_Niz_el_in4.Boost(BoostToHiggs);
			lokalac_Niz_poz_in4.Boost(BoostToHiggs);
			lokalac_Niz_el_out4.Boost(BoostToHiggs);
			lokalac_Niz_poz_out4.Boost(BoostToHiggs);
			lokalacZ1.Boost(BoostToHiggs);
			lokalacZ2.Boost(BoostToHiggs);
			higgs.Boost(BoostToHiggs);

			Niz_el_in3.SetXYZ(lokalac_Niz_el_in4.X(), lokalac_Niz_el_in4.Y(), lokalac_Niz_el_in4.Z());
			Niz_el_out3.SetXYZ(lokalac_Niz_el_out4.X(), lokalac_Niz_el_out4.Y(), lokalac_Niz_el_out4.Z());
			Niz_poz_in3.SetXYZ(lokalac_Niz_poz_in4.X(), lokalac_Niz_poz_in4.Y(), lokalac_Niz_poz_in4.Z());
		  	Niz_poz_out3.SetXYZ(lokalac_Niz_poz_out4.X(), lokalac_Niz_poz_out4.Y(), lokalac_Niz_poz_out4.Z());

			OnShell3_Z.SetXYZ(lokalacZ1.X(),lokalacZ1.Y(),lokalacZ1.Z());		// impuls Z1
			OffShell3_Z.SetXYZ(lokalacZ2.X(),lokalacZ2.Y(),lokalacZ2.Z());	// impuls Z2

			higgsl3.SetXYZ(higgs.X(), higgs.Y(), higgs.Z());

			nz.SetXYZ(0,0,1);   // z-osa

			Double_t thetaH = higgsl3.Theta();
			mH = higgs.M();

			//------ NASA DEFINICIJA UGLA FI - u ovom slucaju su normale na ravni usmerene isto ------//

			TVector3 brojilac1 = Niz_el_in3.Cross(Niz_el_out3);   // brojilac vektora n1
			Double_t imenilac1 = brojilac1.Mag();									// imenilac vektora n1
			TVector3 n1 = brojilac1 * pow(imenilac1,-1);          // vektor n1

			TVector3 brojilac2 = Niz_poz_in3.Cross(Niz_poz_out3); // brojilac vektora n2
			Double_t imenilac2 = brojilac2.Mag();									// imenilac vektora n2
        	TVector3 n2 = brojilac2 * pow(imenilac2,-1);          // vektor n2


			TVector3 brojilacOgawa1 = OnShell3_Z.Cross(Niz_el_in3);   	// brojilac orta n1
			Double_t imenilacOgawa1 = brojilacOgawa1.Mag();							// imenilac orta n1
			TVector3 nOgawa1 = brojilacOgawa1 * pow(imenilacOgawa1,-1); // ort normale n1

			TVector3 brojilacOgawa2 = OffShell3_Z.Cross(Niz_poz_out3);	// brojilac orta n2
			Double_t imenilacOgawa2 = brojilacOgawa2.Mag();							// imenilac orta n2
			TVector3 nOgawa2 = brojilacOgawa2 * pow(imenilacOgawa2,-1);	// ort normale n2

			fi = OnShell3_Z.Dot(n1.Cross(n2)) * fabs (pow(OnShell3_Z.Dot(n1.Cross(n2)),-1)) * (acos(n1.Dot(n2)));

			/*goran kod*/

			   TRotation r ();        // r initialized as identity
			 //  XYZVector u (0,1,0);

			   TVector3 yosa  (0,1,0);




			/*kraj goran kod*/

			if(fi < 0) {	fi = fi + 2 * M_PI; }

			Double_t thetaZ1 = lokalacZ1.Theta();
			theta1 = TMath::Cos(thetaZ1);

			Double_t thetaZ2 = lokalacZ2.Theta();
			theta2 = TMath::Cos(thetaZ2);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			Float_t Ipsilon, Iks;

			Ipsilon = OnShell3_Z.Dot(nOgawa1.Cross(nOgawa2));
			Iks = nOgawa1.Dot(nOgawa2);

//			Float_t TanOgawa = tan(Ipsilon/Iks);
//			fiOgawa = atan(TanOgawa);

			fiOgawa = atan2(Ipsilon,Iks);

			if(fiOgawa < 0) {	fiOgawa = fiOgawa + 2 * M_PI; }
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			leptonTree.Fill();

		} // end of event loop

		lcReader -> close();

	} // Kraj petlje po fajlovima

	cout << "Ukupan broj dogadjaja : " << Nov_dogadjaj << endl;

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	leptonTree.Write();

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

	UInt_t nLastJob = 1;
	if(argc > iarg)
	{
	  nLastJob = atoi(argv[iarg]);
	  iarg++;
	}

	TString fName = "ZZH_fuzija_1.4TeV_1ab_";
	if(argc > iarg)
	{
	  fName = argv[iarg];   // iarg++;
	}
	TString rfName = "probaluk.root";
	if(argc > iarg) rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
