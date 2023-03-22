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


vector <double> pxe;// upisujemo sve iz kolone u niz (vector = array)


int main() // int main NOT void main
{
			TVector3 x(1, 0, 0);
			TVector3 y(0, 1, 0);
			TVector3 z(0, 0, 1);
			TVector3 x1(5,3,10);

		//	y.Print();
			cout<<"___________________________________"<<endl;

			TRotation r;
		//	double angle = M_PI/2;
			double angle = x1.Theta();
			TVector3 ort = x1.Cross(z);
			r.Rotate(angle,ort);

			x.Transform(r);
			y.Transform(r);
			z.Transform(r);

			x1.Transform(r);
			x1.Print();




	/*		TVector3 directionMinus = z;
			TVector3 testy = y.Cross(directionMinus);

			testy.Print();
			cout<<"___________________________________"<<endl;


			rTest.AngleAxis(angle, x1);
			z.Transform(rTest);

			z.Print();*/

/*			TRotation rMinus;
			rMinus.SetZAxis(x, x);
			rMinus.Invert();
			x.Transform(rMinus);
			y.Transform(rMinus);
			z.Transform(rMinus);
*/


	/*			x.Print();
				cout<<"X___________________________________"<<endl;
				y.Print();
				cout<<"Y___________________________________"<<endl;
				z.Print();
				cout<<"Z___________________________________"<<endl;*/

  return 0; // everything went right.
}
