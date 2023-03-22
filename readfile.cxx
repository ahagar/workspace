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
	#include <algorithm>
	#include "TMatrixD.h"

// LCIO includes
	#include "lcio.h"
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

vector <double> E_el;  // дефинише се 1. низ (vector = низ)
vector <double> E_po;  // дефинише се 2. низ (vector = низ)
vector <double> Px_el; // дефинише се 3. низ (vector = низ)
vector <double> Py_el; // дефинише се 4. низ (vector = низ)
vector <double> Pz_el; // дефинише се 5. низ (vector = низ)
vector <double> Px_po; // дефинише се 6. низ (vector = низ)
vector <double> Py_po; // дефинише се 7. низ (vector = низ)
vector <double> Pz_po; // дефинише се 8. низ (vector = низ)

vector <double> Px_e;  // дефинише се 9. низ (vector = низ)
vector <double> Py_e;  // дефинише се 10. низ (vector = низ)
vector <double> Pz_e;  // дефинише се 11. низ (vector = низ)
vector <double> E_e;   // дефинише се 12. низ (vector = низ)
vector <double> Px_p;  // дефинише се 13. низ (vector = низ)
vector <double> Py_p;  // дефинише се 14. низ (vector = низ)
vector <double> Pz_p;  // дефинише се 15. низ (vector = низ)
vector <double> E_p;   // дефинише се 16. низ (vector = низ)

Float_t BhlumiElektronPx, BhlumiElektronPy, BhlumiElektronPz, BhlumiElektronE;
Float_t BhlumiPozitronPx, BhlumiPozitronPy, BhlumiPozitronPz, BhlumiPozitronE;

TTree AjvanTree ("AjvanTree", "Generator particle tree");

AjvanTree.Branch("BhlumiElektronPx", &BhlumiElektronPx, "BhlumiElektronPx");
AjvanTree.Branch("BhlumiElektronPy", &BhlumiElektronPy, "BhlumiElektronPy");
AjvanTree.Branch("BhlumiElektronPz", &BhlumiElektronPz, "BhlumiElektronPz");
AjvanTree.Branch("BhlumiElektronE", &BhlumiElektronE, "BhlumiElektronE");

AjvanTree.Branch("BhlumiPozitronPx", &BhlumiPozitronPx, "BhlumiPozitronPx");
AjvanTree.Branch("BhlumiPozitronPy", &BhlumiPozitronPy, "BhlumiPozitronPy");
AjvanTree.Branch("BhlumiPozitronPz", &BhlumiPozitronPz, "BhlumiPozitronPz");
AjvanTree.Branch("BhlumiPozitronE", &BhlumiPozitronE, "BhlumiPozitronE");

int main() // int main ЈОК void main
{
  ifstream infile1; // дефинише се 1. фајл

  infile1.open("GP.ee_500evt.out"); // отвара се 1. Ајванов фајл

  std::string line1;

	if(infile1.fail()) // проверава се да ли је 1. Ајванов фајл успешно отворен
	{
      cout << "Greška poluteška GuineaPiglet!" << endl;
      return 1; // џаба рмбање ако се фајл не мере отворит'...
	}

  while (std::getline(infile1, line1)) // из фајла се учитава линија по линија
  {
      std::stringstream ss1(line1); // дужина линије

      double a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q;  // променљиве у колонама 1. Ајвановог фајла које су типа DOUBLE

			if (ss1 >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j >> k >> l >> m >> n >> o >> p >> q) // учитавају се све колоне из 1. Ајвановог фајла
      {
				E_el.push_back(a);  // уписују се све вредности из 1. колоне у 1. низ
    	  E_po.push_back(b);  // уписују се све вредности из 2. колоне у 2. низ
    	  Px_el.push_back(g*a); // уписују се све вредности из 7. колоне у 3. низ (Impuls je dobijen mnozenjem brzine iz Ajvanovog fajla sa energijom)
    	  Py_el.push_back(h*a); // уписују се све вредности из 8. колоне у 4. низ (Impuls je dobijen mnozenjem brzine iz Ajvanovog fajla sa energijom)
    	  Pz_el.push_back( TMath::Sqrt(a*a - g*a*g*a - h*a*h*a - 511.e-6) ); // уписују се сви израчунати Pz_el у 5. низ
				// cout << "Vrednost Pz_el = " << TMath::Sqrt(a*a - g*a*g*a - h*a*h*a - 511.e-6) << endl;
    	  Px_po.push_back(i*b); // уписују се све вредности из 9. колоне у 6. низ (Impuls je dobijen mnozenjem brzine iz Ajvanovog fajla sa energijom)
    	  Py_po.push_back(j*b); // уписују се све вредности из 10. колоне у 7. низ (Impuls je dobijen mnozenjem brzine iz Ajvanovog fajla sa energijom)
    	  Pz_po.push_back( TMath::Sqrt(b*b - i*b*i*b - j*b*j*b - 511.e-6) ); // уписују се сви израчунати Pz_po у 8. низ
      }
  }

  ifstream infile2; // дефинише се 2. фајл

	infile2.open("bhabha_Zpole_no_photons_510evt.ini"); // отвара се 2. Ајванов фајл

  std::string line2;

	if(infile2.fail()) // проверава се да ли је 2. Ајванов фајл успешно отворен
	{
      cout << "Greška poluteška Bhabha!" << endl;
      return 1; // џаба рмбање ако се фајл не мере отворит'...
	}

  while (std::getline(infile2, line2)) // из фајла се учитава линија по линија
  {
      std::stringstream ss2(line2); // дужина линије

      double r, s, t, u, v, w, x, y, z, z1;  // променљиве у колонама 2. Ајвановог фајла које су типа DOUBLE

      if (ss2 >> r >> s >> t >> u >> v >> w >> x >> y >> z >> z1) // учитавају се све колоне из 2. Ајвановог фајла
      {
    	  Px_e.push_back(s); // уписују се све вредности из 2. колоне у 9. низ
    	  Py_e.push_back(t); // уписују се све вредности из 3. колоне у 10. низ
    	  Pz_e.push_back(u); // уписују се све вредности из 4. колоне у 11. низ
    	  E_e.push_back(v);  // уписују се све вредности из 8. колоне у 12. низ
    	  Px_p.push_back(w); // уписују се све вредности из 5. колоне у 13. низ
    	  Py_p.push_back(x); // уписују се све вредности из 6. колоне у 14. низ
    	  Pz_p.push_back(y); // уписују се све вредности из 7. колоне у 15. низ
    	  E_p.push_back(z);  // уписују се све вредности из 9. колоне у 16. низ
      }
  }


//sve otvaramo u okviru jednog for-a
/* zadatak sledeci:
saberemo 4vektore elektrona i pozitrona iz GuineaPiglet u 1 4vektor (tacka sudara). Boostujemo 4vektore elektrona i pozitrona iz BHlumi fajla u tacku sudara.
Komponente novih 4vektora iz BHlumija upisemo u ROOT. Drvo se kreira isto kao i u prethodnim programima. BoostToVector takodje ima u polarimetrima.
*/

//	vector <TLorentzVector> GP_elektron, GP_pozitron
	TLorentzVector GP_tacka_sudara;
	TLorentzVector Bhlumi_elektron, Bhlumi_pozitron;

	vector <TVector3> BoostToSudar, BhlumiElektron, BhlumiPozitron;

  for (int i = 0; (int) i < Px_el.size(); i++)
	{
    	GP_elektron.SetPxPyPzE( Px_el[i], Py_el[i], Pz_el[i], E_el[i] );
		GP_pozitron.SetPxPyPzE( Px_po[i], Py_po[i], Pz_po[i], E_po[i] );

		GP_tacka_sudara.SetPxPyPzE( Px_el[i] + Px_po[i], Py_el[i] + Py_po[i], Pz_el[i] + Pz_po[i], E_el[i] + E_po[i] ); //u zagradi su vektori u koje su upisane promenljive.

		Bhlumi_elektron.SetPxPyPzE( Px_e[i], Py_e[i], Pz_e[i], E_e[i] );
		Bhlumi_pozitron.SetPxPyPzE( Px_p[i] Py_p[i], Pz_p[i], E_p[i] );

		BoostToSudar = -( GP_tacka_sudara.BoostVector() );

	//	Bhlumi_elektron.Boost(BoostToSudar);
	//	Bhlumi_pozitron.Boost(BoostToSudar);

/*	BhlumiElektron(Bhlumi_elektron.Px(), Bhlumi_elektron.Py(), Bhlumi_elektron.Pz()) ;
		BhlumiPozitron(Bhlumi_pozitron.Px(),Bhlumi_pozitron.Py(),Bhlumi_pozitron.Pz()) ;    */

		/*BhlumiElektronPx = Bhlumi_elektron.Px();
		BhlumiElektronPy = Bhlumi_elektron.Py();
		BhlumiElektronPz = Bhlumi_elektron.Pz();
		BhlumiElektronE = Bhlumi_elektron.E();

		BhlumiPozitronPx = Bhlumi_pozitron.Px();
		BhlumiPozitronPy = Bhlumi_pozitron.Py();
		BhlumiPozitronPz = Bhlumi_pozitron.Pz();
		BhlumiPozitronE = Bhlumi_pozitron.E();*/

		AjvanTree.Fill();

	  cout << "Ispis impulsa elektrona iz 1. Ajvanovog fajla: Px_el = " << Px_el[i] << ", Py_el = " << Py_el[i] << ", Pz_el = " << Pz_el[i] << endl;
	  cout << "Ispis impulsa pozitrona iz 1. Ajvanovog fajla: Px_po = " << Px_po[i] << ",Py_po = " << Py_po[i] << ", Pz_po = " << Pz_po[i] << endl;
  }

  cout << "   " << endl;
  cout << "----------------------" << endl;
  cout << "   " << endl;

/*  for (int i = 0; (int) i < E_p.size(); i++)
	{
	  cout << "Ispis impulsa elektrona iz 2. Ajvanovog fajla: Px_e = " << Px_e[i] << ", Py_e = " << Py_e[i] << ", Pz_e = " << Pz_e[i] << endl;
	  cout << "Ispis impulsa pozitrona iz 2. Ajvanovog fajla: Px_p = " << Px_p[i] << ", Py_p = " << Py_p[i] << ", Pz_p = " << Pz_p[i] << endl;
  } */

	/* int num = 0; // број мора почети од нуле
  infile.open("lumi.ee_Zpole_500evt.ini"); // у фајлу се налазе 3 колоне бројева

	if(infile.fail()) // проверава се да ли је фајл успешно отворен
	{
      cout << "Greška poluteška!" << endl;
      return 1; // џаба рмбање ако се фајл не мере отворит'...
	}

	while(!infile.eof()) // учитава се фајл до самог краја фајла, не до краја текуће линије
	{
		infile >> exam1[num]; // учитава се број у 1. колони
		infile >> exam2[num]; // учитава се број у 2. колони
		infile >> exam3[num]; // учитава се број у 3. колони

		++num; // иде се на следећи број
		cout << "Ispis "<< exam1[num] << endl;

    // Све ово може се урадити и у једној јединој линији кода, овако:
    // infile >> exam1[num] >> exam2[num] >> exam3[num]; ++num;
	} */


  infile1.close();
  infile2.close();

	TFile rootFile("Ucitavanje","RECREATE");

	rootFile.Write();
	rootFile.Close();

  return 0; // Изгледа да је све оки-доки.
}
