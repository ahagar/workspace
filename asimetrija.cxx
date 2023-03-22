/***********************************************************
 * meri ugao iznedju rekonstruisanog i generisanog fotona, i ako je ugao manji ode 5, daje odnos energija
 * 	Read reconstructed photons from a  .slcio file, and draws histograms for different variables
 *
 *  Author:Goran Kačarević
 *
 *
 *
 ***********************************************************/



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


// LCIO includes
#include "lcio.h"
#include <IOIMPL/LCFactory.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include "EVENT/ReconstructedParticle.h"
#include <UTIL/LCTOOLS.h>
#include <Exceptions.h>
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/LCRelation.h"
#include "IO/LCReader.h"
#endif

#include "stdlib.h"
#include "CandidateData.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cassert>


#include "varList.h"
#include "utils.h"

using namespace std;



Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
#ifdef __CINT__
	gSystem->Load("${LCIO}/lib/liblcio.so");
	gSystem->Load("${LCIO}/lib/liblcioDict.so");
#endif







	TTree eventList("treeMuons", "ILD event list");
	varListIvan vl; /* SL specific */
	eventList.Branch("thetaMuPlus", &(vl.thetaMuPlus), "#theta_{#mu^{+}}");
	eventList.Branch("thetaMuMinus", &(vl.thetaMuMinus), "thetaMinus");
	eventList.Branch("energyMuPlus", &(vl.energyMuPlus), "energyMuPlus");
	eventList.Branch("energyMuMinus", &(vl.energyMuMinus), "energyMuMinus");
	eventList.Branch("ptMuPlus", &(vl.ptMuPlus), "ptPlus");
	eventList.Branch("ptMuMinus", &(vl.ptMuMinus), "ptMinus");
	eventList.Branch("angleMuons", &(vl.angleMuons), "angleMuons");
	eventList.Branch("angleMuonsRad", &(vl.angleMuonsRad), "angleMuonsRad");
	eventList.Branch("deltaTheta", &(vl.deltaTheta), "deltaTheta");
	eventList.Branch("px2muons", &(vl.px2muons), "px2muons");

	Float_t ex =0, ey=0, ez=0, px=0, py=0, pz=0, eth=0, pth=0, een=0, pen=0, ept=0, ppt=0;
	Float_t  exaxis=0, pxaxis=0, eyaxis=0, pyaxis=0, ezaxis=0, pzaxis=0;
	Float_t NL1mm, NR1mm, NL4mm, NR4mm, NLfid, NRfid, NLlumi, NRlumi, NLfid1mm, NLfid4mm, NRfid1mm,NRfid4mm ;



	TTree bhabha("treeBaba", "babe");
	bhabha.Branch("xosa", &(vl.xosa), "xosa");
	bhabha.Branch("yosa", &(vl.yosa), "yosa");
	bhabha.Branch("zosa", &(vl.zosa), "zosa");
	bhabha.Branch("thetaBhabha", &(vl.thetaBhabha), "thetaBhabha");
	bhabha.Branch("pdgbaba", &(vl.pdgbaba), "pdgbaba");

	TTree bhabha_01("treeFid", "babe");
	bhabha_01.Branch("eth", &eth, "eth");
	bhabha_01.Branch("een", &een, "een");
	bhabha_01.Branch("ept", &ept, "ept");
	bhabha_01.Branch("exaxis", &exaxis, "exaxis");
	bhabha_01.Branch("eyaxis", &eyaxis, "eyaxis");
	bhabha_01.Branch("ezaxis", &ezaxis, "ezaxis");


	bhabha_01.Branch("pth", &pth, "pth");
	bhabha_01.Branch("pen", &pen, "pen");
	bhabha_01.Branch("ppt", &ppt, "ppt");
	bhabha_01.Branch("pxaxis", &pxaxis, "pxaxis");
	bhabha_01.Branch("pyaxis", &pyaxis, "pyaxis");
	bhabha_01.Branch("pzaxis", &pzaxis, "pzaxis");

	bhabha_01.Branch("NLlumi", &NLlumi, "NLlumi");
	bhabha_01.Branch("NRlumi", &NRlumi, "NRlumi");
	bhabha_01.Branch("NL1mm", &NL1mm, "NL1mm");
	bhabha_01.Branch("NL4mm", &NL4mm, "NL4mm");
	bhabha_01.Branch("NR1mm", &NR1mm, "NR1mm");
	bhabha_01.Branch("NR4mm", &NR4mm, "NR4mm");


	bhabha_01.Branch("NLfid", &NLfid, "NLfid");
	bhabha_01.Branch("NRfid", &NRfid, "NRfid");
	bhabha_01.Branch("NLfid1mm", &NLfid1mm, "NLfid1mm");
	bhabha_01.Branch("NLfid4mm", &NLfid4mm, "NLfid4mm");
	bhabha_01.Branch("NRfid1mm", &NRfid1mm, "NRfid1mm");
	bhabha_01.Branch("NRfid4mm", &NRfid4mm, "NRfid4mm");





	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

 	double z = 0.950; // 0.95 na CEPC
	Double_t r2_in = 0.028*0.028, r2_out = 0.01;
	Double_t r21_in = 0.038*0.038, r21_out = 0.08*0.08;


	Int_t broj2Dogadjaja = 0;
	Int_t broj2DogadjajaLumi = 0;


	//petlja koja iščitava .slcio podatke ukoliko ima više fajlova za jedan process
	for(UInt_t iJob=nFirstJob; iJob<=nLastJob; iJob++)
	{
		cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);
				try
				{
					lcReader->open(Form("%s%i.slcio", fName.Data(), iJob));
				}
				catch(lcio::IOException &ex)
				{
					cout << ". Could not open.\n"; // Exception " << ex.what() << endl;
					continue;
				}
				cout << ". Reading.\n";


		int brojDogadjaja = lcReader->getNumberOfEvents();	//Ukupan broj događaja
		cout << "Broj dogadjaja po fajlu je  : " << brojDogadjaja << endl;


		std::vector <TLorentzVector> allPhotons; //vektor koji prikuplja sve fotone
		Int_t brojac =0;

		// Prolazimo po svakom dogadjaju
		EVENT::LCEvent* evt = 0;
		while( (evt = lcReader->readNextEvent()) != 0 /*&& brojac < 1000*/)
		{
			brojac++;
		    vector <TLorentzVector> cestice,cesticeLumi, cestice50, cestice60, cestice70, cestice80, cestice90, cestice100, cestice110, cestice120, cestice130;
		    vector <TLorentzVector> cestice140, cestice150, cestice160, cestice170, cestice180, cestice190, cestice200, cestice210, cestice220, cestice230, cestice240, cestice250;
			vector<EVENT::MCParticle*> bhabhas;

		    int testC = 0;
		    Double_t xaxis, yaxis, zaxis, xeminus =0, xeplus=0, yeminus, yeplus, zeminus, zeplus, pdgeminus, pdgeplus, thetaeminus, thetaeplus;
		    bool electron = false;
		    bool positron = false;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			//std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < (int)colNames.size(); ++i)
			{
			//std::cout << colNames[i] << endl;
			}



				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCBhabha");
			//	IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCParticle");
				for (Int_t i = 0; i < mcParticles->getNumberOfElements() ; i++)
				{
					IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles->getElementAt(i);

					if (mcParticle->getGeneratorStatus() != 1) continue;
					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici

					const double *p = mcParticle->getMomentum(); // impuls čestice
					double e = mcParticle->getEnergy();	//energija čestice
					temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor

					Int_t particlePDG = mcParticle->getPDG();
				//	cout<<"particle PDG is: "<<particlePDG<<endl;
					Double_t Pt = temp.Pt();		//promenljiva koja nam daje Pt čestice

					//definisem unutrasnji i spoljasni prsten lumicala  - tacnije kvadrat

				if(particlePDG ==11){
					electron = true;
					bhabhas.push_back(mcParticle);
				}

				if(particlePDG ==-11){
					positron = true;
					bhabhas.push_back(mcParticle);
				}



	/*				Double_t xaxis, yaxis, zaxis;
					if (temp.Pz()<0) z = -z;
					double r = z * tan(temp.Theta());

					xaxis =r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
					yaxis = r* temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py());
					zaxis = z;
					Float_t fid, fid50;

					Double_t pprecnik2 = pow(xaxis,2) + pow(yaxis,2);

					bool el = true, pos=true;
					if (particlePDG == 11 && pprecnik2 < r2_in && pprecnik2>r2_out ){
						el = false;
					}
					if (particlePDG == -11 && pprecnik2 < r2_in && pprecnik2>r2_out ){
						el = false;
					}
					cout << "r elektrona je: " << pprecnik2<<endl;

					vl.pdgbaba = particlePDG;
					vl.thetaBhabha=temp.Theta();

					fid =  pow(xaxis,2)+ pow(yaxis,2);
					fid50= pow(xaxis+0.050,2)+ pow(yaxis,2);

					vl.xosa = r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
					vl.yosa = r * temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
					vl.pdgbaba = particlePDG;

				//	vl.fiducial = fid;
					vl.zosa = z;
					bhabha.Fill();*/




				}  // end of particle loop


				if (electron ==false || positron==false ) continue;
				float NL1 = 0, NLf=0; float NL4= 0;  float NR1= 0; float NR4 = 0, NRf=0, NLl = 0, NRl = 0, NLf1 = 0, NLf4 = 0, NRf1 = 0, NRf4 = 0;
				if (bhabhas.size()!= 2 )continue;
				//	cout<<"usao sam u pravi evt"<<endl;
					TLorentzVector temp; //četvorovektor u koji sakupljamo informacije o svakoj čestici
					for(int i = 0; i < (int) bhabhas.size(); i++){

						const double *p = bhabhas[i]->getMomentum(); // impuls čestice
						double e = bhabhas[i]->getEnergy();	//energija čestice
						temp.SetPxPyPzE(p[0], p[1], p[2], e);  	//zapisujemo vrednosti energije i impulsa u četvorovektor
						Int_t particlePDG = bhabhas[i]->getPDG();

				//		cout <<"pdg cestice: "<<particlePDG<<endl;
				//		cout<<"pre boosta px: "<<temp.Px()<<", py: "<< temp.Py()<<", pz: "<<temp.Pz()<<", E: "<<temp.E()<<endl;


				//	if (particlePDG==11)	temp.Boost(0 ,0 , 0.5);
					//if (particlePDG==-11)	temp.Boost(0 ,0 , -0.5);
					temp.Boost(0 ,0 , 0.0002);

				//	cout<<"posle  boosta px: "<<temp.Px()<<", py: "<< temp.Py()<<", pz: "<<temp.Pz()<<", E: "<<temp.E()<<endl;
					//cout<<"_________________________________"<<endl;


						double pt = temp.Pt();
						double theta = temp.Theta();

						if (temp.Pz()<0) z = -z;
						double r = z * tan(temp.Theta());

						xaxis =r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						yaxis = r* temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py());
						zaxis = z;

						if (particlePDG ==11){
							exaxis=xaxis;
							eyaxis=yaxis;
							ezaxis=z;
						//	ept= pt;
						//	eth = theta;
					//		een = e;

						}

						if (particlePDG ==-11){
							pxaxis=xaxis;
							pyaxis=yaxis;
							pzaxis=z;
					//		ppt = pt;
					//		pth = theta;
					//		pen = e;
						}
					//	cout <<"pdg čestice: "<<particlePDG<<", xaxis: "<< xaxis<<", "<<exaxis<<endl;

						//leva asimetrija
						if (z<0)
						{
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0295,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2)) NL1 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0325,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2)) NL4 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0285,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2) ) NLl =1;

							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.051,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2)) NLf1 = 1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.054,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2)) NLf4 = 1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0500,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2))  NLf = 1;

						}
						//desna asimetrija
						if (z>0)
						{
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0295,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2)) NR1 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0325,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2)) NR4 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0285,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.1,2)) NRl =1;

							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.051,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2)) NRf1 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.054,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2)) NRf4 =1;
							if (pow(xaxis,2)+ pow(yaxis,2)> pow(0.0500,2) && pow(xaxis,2)+ pow(yaxis,2) < pow(0.0752,2))  NRf = 1;

						}

					//	cout << "r elektrona je: " << pprecnik2<<endl;

						vl.pdgbaba = particlePDG;
						vl.thetaBhabha=temp.Theta();


						vl.xosa = r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						vl.yosa = r * temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						vl.pdgbaba = particlePDG;

					//	vl.fiducial = fid;
						vl.zosa = z;
					//	bhabha.Fill();
					//	if(pprecnik2 < r21_in || pprecnik2 > r21_out) continue;
						//cout <<  "r21 in: "<< 0.028*0.028<<", r cestice :"<<pow (xaxis,2)+ pow(yaxis,2)<<", r21_out: "<<0.01<<endl;
					//	bhabha.Fill();

					}
					NL1mm = NL1;
					NL4mm = NL4;
					NR1mm = NR1;
					NR4mm = NR4;
					NLlumi = NLl;
					NRlumi = NRl;

					NLfid1mm = NLf1;
					NLfid4mm = NLf4;
					NRfid1mm = NRf1;
					NRfid4mm = NRf4;
					NLfid = NLf;
					NRfid = NRf;



					bhabha_01.Fill(); // KG









		} // End of event loop


		lcReader->close();
		cout << "ukupan broj događaja po fajlu je: "<< brojac<<endl;


	} // End of file loop


//	cout << " broj događaja sa dve čestice u fiducialnoj zaprenini je: "<< broj2Dogadjaja<<endl;


	ofstream file ;
	file.open("proba.txt", ios_base::out | ios_base::app);

	file << "Lumi" << "\t" << broj2DogadjajaLumi<<endl;
/*	file << "0" << "\t" << broj2Dogadjaja << "\t"<< up<< "\t"<<down << "\t"<<up-down<<endl;
	file << "50" << "\t" << broj2Dogadjaja50 << "\t"<< up50<< "\t"<<down50 << "\t"<<up50-down50<< endl;
	*/
	file << "___________________________________ "<< endl;

	file.close();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
//	eventList.Write();
    bhabha.Write();
    bhabha_01.Write();

	rootFile.Close();



	return 0;
}



Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;
	UInt_t nFirstJob = 1;
	if(argc>iarg) nFirstJob = atoi(argv[iarg]); iarg++;
	UInt_t nLastJob = 10;
	if(argc>iarg) nLastJob   = atoi(argv[iarg]); iarg++;

	TString fName = "ee_rem_col_";
	if(argc>iarg) fName = argv[iarg]; iarg++;

	TString rfName = "bhabha.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}

