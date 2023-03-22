/***********************************************************
 * 
 * 	Read reconstructed muons from a  .slcio file
 *  and plot some distributions
 *
 *
 *  Author G.Milutinovic-Dumbelovic
 *
 ***********************************************************/

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TVectorT.h"
#include "TCanvas.h" 
#include "TStyle.h"
#include "TPostScript.h"
#include "TLatex.h"
#include "TLorentzVector.h" 
#include "TTree.h"

// LCIO includes
#include "lcio.h"
#include <IOIMPL/LCFactory.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <Exceptions.h>
#endif

#include "stdlib.h"
#include <sstream>
#include "varList.h"
#include "hmumu_pdf.h"
//#include "hmumu_pdf.C"
#include <sstream>

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooArgSet.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"

using namespace std;





Int_t main(int argc, char* argv[])
{
	cout << "testing"<<endl;


}


