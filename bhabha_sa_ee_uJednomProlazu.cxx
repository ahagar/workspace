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


	TTree bhabha("treeBaba", "babe");
	bhabha.Branch("xosa", &(vl.xosa), "xosa");
	bhabha.Branch("yosa", &(vl.yosa), "yosa");
	bhabha.Branch("zosa", &(vl.zosa), "zosa");
	bhabha.Branch("thetaBhabha", &(vl.thetaBhabha), "thetaBhabha");
	bhabha.Branch("pdgbaba", &(vl.pdgbaba), "pdgbaba");
	bhabha.Branch("fiducial", &(vl.fiducial), "fiducial");




	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader() ;
	TString fName = fn;
	stringstream fNameStream;

 	double z = 950; // 0.95 na CEPS
	Int_t broj2Dogadjaja = 0;
	Int_t broj2DogadjajaLumi = 0;
	Int_t broj2Dogadjaja50 = 0;
	Int_t broj2Dogadjaja60 = 0;
	Int_t broj2Dogadjaja70 = 0;
	Int_t broj2Dogadjaja80 = 0;
	Int_t broj2Dogadjaja90 = 0;
	Int_t broj2Dogadjaja100 = 0;
	Int_t broj2Dogadjaja110 = 0;
	Int_t broj2Dogadjaja120 = 0;
	Int_t broj2Dogadjaja130 = 0;
	Int_t broj2Dogadjaja140 = 0;
	Int_t broj2Dogadjaja150 = 0;
	Int_t broj2Dogadjaja160 = 0;
	Int_t broj2Dogadjaja170 = 0;
	Int_t broj2Dogadjaja180 = 0;
	Int_t broj2Dogadjaja190 = 0;
	Int_t broj2Dogadjaja200 = 0;
	Int_t broj2Dogadjaja210 = 0;
	Int_t broj2Dogadjaja220 = 0;
	Int_t broj2Dogadjaja230 = 0;
	Int_t broj2Dogadjaja240 = 0;
	Int_t broj2Dogadjaja250 = 0;

	int up =0, up50 =0, up60=0, up70=0, up80=0, up90=0, up100=0, up110=0, up120=0, up130=0, up140=0, up150=0, up160=0, up170=0, up180=0, up190=0, up200=0, up210=0, up220=0, up230=0, up240=0, up250=0;
	int down=0, down50=0, down60=0, down70=0, down80=0, down90=0, down100=0, down110=0, down120=0, down130=0, down140=0, down150=0, down160=0, down170=0, down180=0, down190=0, down200=0, down210=0, down220=0, down230=0, down240=0, down250=0;





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

		    int testC = 0;



				std::vector<std::string> colNames = *evt->getCollectionNames();
			//std::cout << "\n\nCollection names: \n";
			for (int i = 0; i < (int)colNames.size(); ++i)
			{
			//std::cout << colNames[i] << endl;
			}



				IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt->getCollection("MCBhabha");
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
			//	cout << "energija cestice  "<<temp.E()<<endl;

					//cout << temp.Theta()<< ",   "<<atan(temp.Y()/temp.X())<<endl;

					if (particlePDG ==11 || particlePDG == -11)
					{

						if (temp.Pz()<0) z = -z;

						double r = z * tan(temp.Theta());
						Float_t fid, fid50, fid60, fid70, fid80, fid90, fid100, fid110, fid120, fid130, fid140, fid150, fid160, fid170, fid180, fid190, fid200, fid210, fid220, fid230, fid240, fid250  ;
						Double_t xaxis =r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						Double_t yaxis = r* temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py());


				    //	vl.xosa = r * cosFi;
					//	vl.yosa = r * sinFi;
				    	vl.xosa = r * temp.Px()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						vl.yosa = r * temp.Py()/sqrt(temp.Px()*temp.Px() + temp.Py()*temp.Py()) ;
						vl.pdgbaba = particlePDG;
						vl.thetaBhabha=temp.Theta();



						fid =  pow(xaxis,2)+ pow(yaxis,2);
						fid50= pow(xaxis+0.050,2)+ pow(yaxis,2);
						fid60= pow(xaxis+0.060,2)+ pow(yaxis,2);
						fid70= pow(xaxis+0.070,2)+ pow(yaxis,2);
						fid80= pow(xaxis+0.080,2)+ pow(yaxis,2);
						fid90= pow(xaxis+0.090,2)+ pow(yaxis,2);
						fid100= pow(xaxis+0.100,2)+ pow(yaxis,2);
						fid110= pow(xaxis+0.110,2)+ pow(yaxis,2);
						fid120= pow(xaxis+0.120,2)+ pow(yaxis,2);
						fid130= pow(xaxis+0.130,2)+ pow(yaxis,2);
						fid140= pow(xaxis+0.140,2)+ pow(yaxis,2);
						fid150= pow(xaxis+0.150,2)+ pow(yaxis,2);
						fid160= pow(xaxis+0.160,2)+ pow(yaxis,2);
						fid170= pow(xaxis+0.170,2)+ pow(yaxis,2);
						fid180= pow(xaxis+0.180,2)+ pow(yaxis,2);
						fid190= pow(xaxis+0.190,2)+ pow(yaxis,2);
						fid200= pow(xaxis+0.200,2)+ pow(yaxis,2);
						fid210= pow(xaxis+0.210,2)+ pow(yaxis,2);
						fid220= pow(xaxis+0.220,2)+ pow(yaxis,2);
						fid230= pow(xaxis+0.230,2)+ pow(yaxis,2);
						fid240= pow(xaxis+0.240,2)+ pow(yaxis,2);
						fid250= pow(xaxis+0.250,2)+ pow(yaxis,2);



						vl.fiducial = fid;
					//	if (temp.Pz()<0) z = -z;
						vl.zosa = z;

						if (fid > 2*25  &&fid< 100*100){
							cesticeLumi.push_back(temp);
						}

						if (fid > 50*50  &&fid< 75*75){
							cestice.push_back(temp);
						}


						if (fid50 > 50*50  &&fid50< 75*75){
							cestice50.push_back(temp);
							if (z>0 && xaxis>0 ) up50++;
							if (z>0 && xaxis<0 ) down50++;
						}else continue;
						if (fid60 > 50*50  &&fid60< 75*75){
							cestice60.push_back(temp);
							if (z>0 && xaxis>0 ) up60++;
							if (z>0 && xaxis<0 ) down60++;
						}
						if (fid70 > 50*50  &&fid70 < 75*75){
							cestice70.push_back(temp);
							if (z>0 && xaxis>0 ) up70++;
							if (z>0 && xaxis<0 ) down70++;
						}
						if (fid80 > 50*50  &&fid80 < 75*75){
							cestice80.push_back(temp);
							if (z>0 && xaxis>0 ) up80++;
							if (z>0 && xaxis<0 ) down80++;
						}
						if (fid90 > 50*50  &&fid90< 75*75){
							cestice90.push_back(temp);
							if (z>0 && xaxis>0 ) up90++;
							if (z>0 && xaxis<0 ) down90++;
						}
						if (fid100 > 50*50  &&fid100< 75*75){
							cestice100.push_back(temp);
							if (z>0 && xaxis>0 ) up100++;
							if (z>0 && xaxis<0 ) down100++;
						}
						if (fid110 > 50*50  &&fid110< 75*75){
							cestice110.push_back(temp);
							if (z>0 && xaxis>0 ) up110++;
							if (z>0 && xaxis<0 ) down110++;
						}
						if (fid120 > 50*50  &&fid120< 75*75) {
							cestice120.push_back(temp);
							if (z>0 && xaxis>0 ) up120++;
							if (z>0 && xaxis<0 ) down120++;
						}
						if (fid130 > 50*50  &&fid130< 75*75){
							cestice130.push_back(temp);
							if (z>0 && xaxis>0 ) up130++;
							if (z>0 && xaxis<0 ) down130++;
						}
						if (fid140 > 50*50  &&fid140< 75*75){
							cestice140.push_back(temp);
							if (z>0 && xaxis>0 ) up140++;
							if (z>0 && xaxis<0 ) down140++;
						}
						if (fid150 > 50*50  &&fid150< 75*75){
							cestice150.push_back(temp);
							if (z>0 && xaxis>0 ) up150++;
							if (z>0 && xaxis<0 ) down150++;
						}
						if (fid160 > 50*50  &&fid160< 75*75) {
							cestice160.push_back(temp);
							if (z>0 && xaxis>0 ) up160++;
							if (z>0 && xaxis<0 ) down160++;
						}
						if (fid170 > 50*50  &&fid170< 75*75) {
							cestice170.push_back(temp);
							if (z>0 && xaxis>0 ) up170++;
							if (z>0 && xaxis<0 ) down170++;
						}
						if (fid180 > 50*50  &&fid180< 75*75){
							cestice180.push_back(temp);
							if (z>0 && xaxis>0 ) up180++;
							if (z>0 && xaxis<0 ) down180++;
						}
						if (fid190 > 50*50  &&fid190< 75*75){
							cestice190.push_back(temp);
							if (z>0 && xaxis>0 ) up190++;
							if (z>0 && xaxis<0 ) down190++;
						}
						if (fid200 > 50*50  &&fid200< 75*75){
							cestice200.push_back(temp);
							if (z>0 && xaxis>0 ) up200++;
							if (z>0 && xaxis<0 ) down200++;
						}
						if (fid210 > 50*50  &&fid210< 75*75){
							cestice210.push_back(temp);
							if (z>0 && xaxis>0 ) up210++;
							if (z>0 && xaxis<0 ) down210++;
						}
						if (fid220 > 50*50  &&fid220< 75*75) {
							cestice220.push_back(temp);
							if (z>0 && xaxis>0 ) up220++;
							if (z>0 && xaxis<0 ) down220++;
						}
						if (fid230 > 50*50  &&fid230< 75*75){
							cestice230.push_back(temp);
							if (z>0 && xaxis>0 ) up230++;
							if (z>0 && xaxis<0 ) down230++;
						}
						if (fid240 > 50*50  &&fid240< 75*75) {
							cestice240.push_back(temp);
							if (z>0 && xaxis>0 ) up240++;
							if (z>0 && xaxis<0 ) down240++;
						}
						if (fid250 > 50*50 &&fid250< 75*75) {
							cestice250.push_back(temp);
							if (z>0 && xaxis>0 ) up250++;
							if (z>0 && xaxis<0 ) down250++;
						}



					//	cout <<"fid: "<<fid<<endl;
						bhabha.Fill();


					}



			//		cout << muplus.E()<<endl;

				}  // end of particle loop

				if (cesticeLumi.size() == 2) broj2DogadjajaLumi++;
				if (cestice.size() == 2) broj2Dogadjaja++;
				if (cestice50.size() == 2) broj2Dogadjaja50++;
				if (cestice60.size() == 2) broj2Dogadjaja60++;
				if (cestice70.size() == 2) broj2Dogadjaja70++;
				if (cestice80.size() == 2) broj2Dogadjaja80++;
				if (cestice90.size() == 2) broj2Dogadjaja90++;
				if (cestice100.size() == 2) broj2Dogadjaja100++;
				if (cestice110.size() == 2) broj2Dogadjaja110++;
				if (cestice120.size() == 2) broj2Dogadjaja120++;
				if (cestice130.size() == 2) broj2Dogadjaja130++;
				if (cestice140.size() == 2) broj2Dogadjaja140++;
				if (cestice150.size() == 2) broj2Dogadjaja150++;
				if (cestice160.size() == 2) broj2Dogadjaja160++;
				if (cestice170.size() == 2) broj2Dogadjaja170++;
				if (cestice180.size() == 2) broj2Dogadjaja180++;
				if (cestice190.size() == 2) broj2Dogadjaja190++;
				if (cestice200.size() == 2) broj2Dogadjaja200++;
				if (cestice210.size() == 2) broj2Dogadjaja210++;
				if (cestice220.size() == 2) broj2Dogadjaja220++;
				if (cestice230.size() == 2) broj2Dogadjaja230++;
				if (cestice240.size() == 2) broj2Dogadjaja240++;
				if (cestice250.size() == 2) broj2Dogadjaja250++;

		} // End of event loop


		lcReader->close();
		cout << "ukupan broj događaja po fajlu je: "<< brojac<<endl;


	} // End of file loop


	cout << " broj događaja sa dve čestice u fiducialnoj zaprenini je: "<< broj2Dogadjaja<<endl;


	ofstream file ;
	file.open("proba.txt", ios_base::out | ios_base::app);

	file << "Lumi" << "\t" << broj2DogadjajaLumi<<endl;
	file << "0" << "\t" << broj2Dogadjaja << "\t"<< up<< "\t"<<down << "\t"<<up-down<<endl;
	file << "50" << "\t" << broj2Dogadjaja50 << "\t"<< up50<< "\t"<<down50 << "\t"<<up50-down50<< endl;
	file << "60" << "\t" << broj2Dogadjaja60 << "\t"<< up60<< "\t"<<down60 << "\t"<<up60-down60<< endl;
	file << "70" << "\t" << broj2Dogadjaja70 << "\t"<< up70<< "\t"<<down70 << "\t"<<up70-down70<< endl;
	file << "80" << "\t" << broj2Dogadjaja80 << "\t"<< up80<< "\t"<<down80 << "\t"<<up80-down80<< endl;
	file << "90" << "\t" << broj2Dogadjaja90 << "\t"<< up90<< "\t"<<down90 << "\t"<<up90-down90<< endl;
	file << "100" << "\t" << broj2Dogadjaja100 << "\t"<< up100<< "\t"<<down100 << "\t"<<up100-down100<< endl;
	file << "110" << "\t" << broj2Dogadjaja110 << "\t"<< up110<< "\t"<<down110 << "\t"<<up110-down110<< endl;
	file << "120" << "\t" << broj2Dogadjaja120 <<"\t"<< up120<< "\t"<<down120 << "\t"<<up120-down120<< endl;
	file << "130" << "\t" << broj2Dogadjaja130 <<"\t"<< up130<< "\t"<<down130 << "\t"<<up130-down130<< endl;
	file << "140" << "\t" << broj2Dogadjaja140 <<"\t"<< up140<< "\t"<<down140 << "\t"<<up140-down140<< endl;
	file << "150" << "\t" << broj2Dogadjaja150 <<"\t"<< up150<< "\t"<<down150 << "\t"<<up150-down150 << endl;
	file << "160" << "\t" << broj2Dogadjaja160 <<"\t"<< up160<< "\t"<<down160 << "\t"<<up160-down160<< endl;
	file << "170" << "\t" << broj2Dogadjaja170<< "\t"<< up170<< "\t"<<down170 << "\t"<<up170-down170 << endl;
	file << "180" << "\t" << broj2Dogadjaja180<< "\t"<< up180<< "\t"<<down180 << "\t"<<up180-down180<< endl;
	file << "190" << "\t" << broj2Dogadjaja190<< "\t"<< up190<< "\t"<<down190 << "\t"<<up190-down190<< endl;
	file << "200" << "\t" << broj2Dogadjaja200<< "\t"<< up200<< "\t"<<down200 << "\t"<<up200-down200<< endl;
	file << "210" << "\t" << broj2Dogadjaja210 <<"\t"<< up210<< "\t"<<down210 << "\t"<<up210-down210<< endl;
	file << "220" << "\t" << broj2Dogadjaja220 <<"\t"<< up220<< "\t"<<down220 << "\t"<<up220-down220<< endl;
	file << "230" << "\t" << broj2Dogadjaja230<< "\t"<< up230<< "\t"<<down230 << "\t"<<up230-down230<< endl;
	file << "240" << "\t" << broj2Dogadjaja240<< "\t"<< up240<< "\t"<<down240 << "\t"<<up240-down240<< endl;
	file << "250" << "\t" << broj2Dogadjaja250 <<"\t"<< up250<< "\t"<<down250 << "\t"<<up250-down250<< endl;
	file << "___________________________________ "<< endl;

	file.close();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");
	eventList.Write();
    bhabha.Write();




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

	TString rfName = "anglemuons.root";
	if(argc>iarg) rfName = argv[iarg]; iarg++;

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}

