/*
 * Polarimetri za 350 GeV (Jinx)
 *
 *  Napravljeno: 7530. grozdober 27.
 *  Autori: KG, Gordana
 */

#ifndef __CINT__
	#include "TROOT.h"
	#include "TFile.h"
	#include "Riostream.h"
	#include "TFile.h"
	#include "TH1.h"
	#include "TH2.h"
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
	#include <algorithm>
	#include "TMatrixD.h"

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
	#include "TMatrixD.h"
#endif

//#include "EventShape.h"
//#include "jama_eig.h"  // Computes eigenvalues and eigenvectors of a real (non-complex) matrix.
//#include "tnt_math_utils.h"
#include "stdlib.h"
#include <sstream>
#include <iostream>
#include <iterator>
#include <fstream>
#include <vector>
#include <array>

using namespace std;

// --------------------------------------------------------------------------------------------------------------------

EVENT::MCParticle* getMCPion1(EVENT::LCEvent* evt)
{
	EVENT::MCParticle* pointer_to_MCpion1 = 0;

	std::vector<std::string> colNames = *evt -> getCollectionNames();

	IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");

	Int_t N_piona = 0;
	Int_t N_neutrina = 0;
	Int_t N_event = 0;

	for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
	{ // 1. Petlja za broj elemenata
		IMPL::MCParticleImpl* generated = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

		if (abs(generated -> getPDG() == 25))
		{// 1.1
			cout << "Nadjen je Higgs \n";

			if (abs (generated -> getDaughters()[0] -> getPDG() == 25))
			{// 1.2 Traze se tauoni
				const EVENT::MCParticleVec & Tauoni = generated -> getDaughters()[0] -> getDaughters();

				if (Tauoni.size() == 2)
				{// 2.1
					EVENT::MCParticle* prviT = (EVENT::MCParticle*) Tauoni[0];
					EVENT::MCParticle* drugiT = (EVENT::MCParticle*) Tauoni[1];

					if (prviT -> getPDG() == 15 && drugiT -> getPDG() == -15)
					{// 2.2 petlja po 1. i 2. tauonu
						cout << "Nadjeni su tauoni" << endl;

						if (prviT -> getDaughters().size() >= 2 && drugiT -> getDaughters().size() >= 2)
						{// 2.3
							EVENT::MCParticle* prviTpotomak1 = (EVENT::MCParticle*) prviT -> getDaughters()[0]; // prvi potomak T1
							EVENT::MCParticle* prviTpotomak2 = (EVENT::MCParticle*) prviT -> getDaughters()[1]; // drugi potomak T1
							EVENT::MCParticle* drugiTpotomak1 = (EVENT::MCParticle*) drugiT -> getDaughters()[0]; // prvi potomak T2
							EVENT::MCParticle* drugiTpotomak2 = (EVENT::MCParticle*) drugiT -> getDaughters()[1]; // drugi potomak T2

							if (abs (prviTpotomak1 -> getPDG()) == 211) N_piona++;
							if (abs (prviTpotomak2 -> getPDG()) == 16) N_neutrina++;
							if (abs (drugiTpotomak1 -> getPDG()) == 211) N_piona++;
							if (abs (drugiTpotomak2 -> getPDG()) == 16) N_neutrina++;

							if (N_piona == 2 && N_neutrina == 2)
							{// Uslov kada imamo 2 piona i 2 neutrina
								N_event++;
								cout << "Nadjen Tauon-Pion dogadjaj broj " << N_event << endl;

								if ((prviTpotomak1 -> getPDG() == -211 && prviTpotomak2 -> getPDG() == 16) || (prviTpotomak1 -> getPDG() == 16 && prviTpotomak2 -> getPDG() == -211))
								{// 2.4 Petlja po potomcima 1. tauona
									if (prviTpotomak1 -> getGeneratorStatus() == 1) // Ako je cestica finalna.
									{
										cout << "Nadjen 1. potomak 1. tauona " << prviTpotomak1 -> getPDG() << " cestica je finalna." << endl;
										pointer_to_MCpion1 = prviTpotomak1;
									}

									if (prviTpotomak2 -> getGeneratorStatus() == 1)
									{
										cout << "Nadjen je 2. potomak 1. tauona " << prviTpotomak2 -> getPDG() << " cestica je finalna." << endl;
										pointer_to_MCpion1 = prviTpotomak2;
									}

									if ((prviTpotomak1 -> getGeneratorStatus() != 1) && (prviTpotomak2 -> getGeneratorStatus() != 1))
									{// grananje za generatorske statuse 1. i 2. potomka 1. tauona

										if (prviTpotomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak nije konacni.
										{// 2.5 Petlja za 1. potomka 1. tauona

											for (Int_t i = 0; (int)i < prviTpotomak1 -> getDaughters().size(); i++)
											{
												EVENT::MCParticle* prviTpotomak1cestica94 = (EVENT::MCParticle*) prviTpotomak1 -> getDaughters()[i]; // Ako je 1. potomak 1. tauona "cestica" 94

												if (abs(prviTpotomak1cestica94 -> getPDG()) == 94)
												{// 2.5.1

													if (prviTpotomak1cestica94 -> getGeneratorStatus() != 1)
													{// 2.5.2 Petlja po potomcima pseudo-cestice 94
														EVENT::MCParticle* prviTpotomak1cestica94potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94 -> getDaughters()[0]; // 1. potomak cestice 94
														EVENT::MCParticle* prviTpotomak1cestica94potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94 -> getDaughters()[1]; // 2. potomak cestice 94

														if ((prviTpotomak1cestica94potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1 -> getGeneratorStatus() == 1))
														{// 2.5.3
															cout << "Nadjen je 1. potomak cestice 94 1. tauona" << prviTpotomak1cestica94potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak1cestica94potomak1;
														}// 2.5.3

														if ((prviTpotomak1cestica94potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2 -> getGeneratorStatus() == 1))
														{// 2.5.3
															cout << "Nadjen je 2. potomak cestice 94 1. tauona " << prviTpotomak1cestica94potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak1cestica94potomak2;
														}// 2.5.3

														if(prviTpotomak1cestica94potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak cestice 94 nije finalni.
														{//2.5.4 Petlja za 1. potomak cestice 94
															EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka cestice 94
															EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka cestice 94

															if ((prviTpotomak1cestica94potomak1potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1 -> getGeneratorStatus() == 1))
															{// 2.5.5 Petlja po potomcima 1. potomka cestice 94
																cout << "Nadjen je 1. potomak 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1;
															}// 2.5.5

															if ((prviTpotomak1cestica94potomak1potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak2 -> getPDG() == -16) && (prviTpotomak1cestica94potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.5 Petlja po potomcima 1. potomka cestice 94
																cout << "Nadjen je 2. potomak 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak2;
															}// 2.5.5

															if(prviTpotomak1cestica94potomak1potomak1 -> getGeneratorStatus() != 1)
															{// 2.5.6 Ako 1. potomak 1. potomka cestice 94 nije finalni.
																EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka cestice 94
																EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka cestice 94

																if ((prviTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																{// 2.5.7 Petlja po potomcima 1. potomka 1. potomka cestice 94
																	cout << "Nadjen je 1. potomak 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak1;
																}// 2.5.7

																if ((prviTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
																{// 2.5.7 Petlja po potomcima 1. potomka 1. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak2;
																}// 2.5.7

																if (prviTpotomak1cestica94potomak1potomak1potomak1 -> getGeneratorStatus() != 1)
																{//2.5.8
																	EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka 1. potomka cestice 94
																	EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka 1. potomka cestice 94

																	if ((prviTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																	{// 2.5.9
																		cout << "Nadjen je 1. potomak 1. potomka 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak1potomak1;
																	} // 2.5.9

																	if ((prviTpotomak1cestica94potomak1potomak1potomak1potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak1potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
																	{// 2.5.9
																		cout << "Nadjen je 2. potomak 1. potomka 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak1potomak2;
																	} // 2.5.9
																} // 2.5.8

																if (prviTpotomak1cestica94potomak1potomak1potomak2 -> getGeneratorStatus() != 1)
																{//2.5.8
																	EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 1. potomka 1. potomka cestice 94
																	EVENT::MCParticle* prviTpotomak1cestica94potomak1potomak1potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak1potomak1potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 1. potomka 1. potomka cestice 94

																	if ((prviTpotomak1cestica94potomak1potomak1potomak2potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak2potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak2potomak1 -> getGeneratorStatus() == 1))
																	{// 2.5.9
																		cout << "Nadjen je 1. potomak 2. potomka 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak2potomak1;
																	} // 2.5.9

																	if ((prviTpotomak1cestica94potomak1potomak1potomak2potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak1potomak1potomak2potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak1potomak1potomak2potomak2 -> getGeneratorStatus() == 1))
																	{// 2.5.9
																		cout << "Nadjen je 2. potomak 2. potomka 1. potomka 1. potomka cestice 94 1. tauona " << prviTpotomak1cestica94potomak1potomak1potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak1potomak1potomak2potomak2;
																	}// 2.5.9
																}// 2.5.8
															}// 2.5.6
														}// 2.5.4 Kraj petlje za 1. potomka cestice 94

														if(prviTpotomak1cestica94potomak2 -> getGeneratorStatus() != 1)  // Ako 2. potomak cestice 94 nije finalni.
														{// 2.6.4 Petlja za 2. potomak "cestice" 94
															EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka cestice 94
															EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka cestice 94

															if ((prviTpotomak1cestica94potomak2potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak1 -> getGeneratorStatus() == 1))
															{// 2.6.5 Petlja po potomcima 2. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak1;
															}// 2.6.5

															if ((prviTpotomak1cestica94potomak2potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2 -> getGeneratorStatus() == 1))
															{// 2.6.5 Petlja po potomcima 2. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2;
															}// 2.6.5

															if(prviTpotomak1cestica94potomak2potomak1 -> getGeneratorStatus() != 1)
															{// 2.6.6
																EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 2. potomka "cestice" 94
																EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 2. potomka "cestice" 94

																if ((prviTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
																{// 2.6.7 Petlja po potomcima 1. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 1. potomak 1. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak1potomak1;
																}// 2.6.7

																if ((prviTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
																{//2.6.7 Petlja po potomcima 1. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 1. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak1potomak2;
																}//2.6.7
															}// 2.6.6

															if (prviTpotomak1cestica94potomak2potomak2 -> getGeneratorStatus() != 1)
															{// 2.6.8
																EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 2. potomka "cestice" 94
																EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 2. potomka "cestice" 94

																if ((prviTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak1 -> getGeneratorStatus() == 1))
																{// 2.6.5 Petlja po potomcima 2. potomka "cestice" 94
																	cout << "Nadjen je 1. potomak 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak1;
																} // 2.6.5

																if ((prviTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak2 -> getGeneratorStatus() == 1))
																{// 2.6.5 Petlja po potomcima 2. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak2;
																}// 2.6.5

																if (prviTpotomak1cestica94potomak2potomak2potomak1 -> getGeneratorStatus() != 1)
																{// 2.6.6
																	EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2potomak1 -> getDaughters()[0];
																	EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2potomak1 -> getDaughters()[1];

																	if ((prviTpotomak1cestica94potomak2potomak2potomak1potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak1potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
																	{// 2.6.7
																		cout << "Nadjen je 1. potomak 1. potomka 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak1potomak1;
																	}// 2.6.7

																	if ((prviTpotomak1cestica94potomak2potomak2potomak1potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak1potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
																	{// 2.6.7
																		cout << "Nadjen je 2. potomak 1. potomka 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak1potomak2;
																	}// 2.6.7
																}// 2.6.6

																if (prviTpotomak1cestica94potomak2potomak2potomak2 -> getGeneratorStatus() != 1)
																{// 2.6.6
																	EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2potomak2 -> getDaughters()[0];
																	EVENT::MCParticle* prviTpotomak1cestica94potomak2potomak2potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak1cestica94potomak2potomak2potomak2 -> getDaughters()[1];

																	if ((prviTpotomak1cestica94potomak2potomak2potomak2potomak1 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak2potomak1 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak2potomak1 -> getGeneratorStatus() == 1))
																	{// 2.6.7
																		cout << "Nadjen je 1. potomak 2. potomka 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak2potomak1;
																	}// 2.6.7

																	if ((prviTpotomak1cestica94potomak2potomak2potomak2potomak2 -> getPDG() == -211 || prviTpotomak1cestica94potomak2potomak2potomak2potomak2 -> getPDG() == 16) && (prviTpotomak1cestica94potomak2potomak2potomak2potomak2 -> getGeneratorStatus() == 1))
																	{// 2.6.7
																		cout << "Nadjen je 2. potomak 2. potomka 2. potomka 2. potomka cestice 94 " << prviTpotomak1cestica94potomak2potomak2potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion1 = prviTpotomak1cestica94potomak2potomak2potomak2potomak2;
																	}// 2.6.7
																}// 2.6.6
															}// 2.6.8
														}// 2.6.4 Zatvorena petlja za 2. potomka "cestice" 94
													}// 2.5.2
												}// 2.5.1
											} // Kraj FOR petlje za potomke 1. potomka 1. tauona
										}// 2.5 Kraj IF za 1. potomka 1. tauona
									} // Kraj IF generatorskog statusa za 1. i 2. potomak 1. tauona

									// ----------------------------------------------------------------

									if (prviTpotomak1 -> getDaughters().size() == 0) // Ako 1. potomak 1. tauona nema potomaka.
									{
										if (prviTpotomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak nije krajnji.
										{// 2.5 Petlja za 2. potomka 1. tauona

											if (prviTpotomak2 -> getDaughters().size() == 0) continue;

											EVENT::MCParticle* prviTpotomak2cestica94 = (EVENT::MCParticle*) prviTpotomak2 -> getDaughters()[0]; // 1. potomak 2.potomka 1. tauona, "cestica" 94

											if (abs(prviTpotomak2cestica94 -> getPDG()) == 94)
											{// 2.5.1

												if (prviTpotomak2cestica94 -> getGeneratorStatus() != 1)
												{//2.5.2 Uslov po potomcima "cestice" 94
													EVENT::MCParticle* prviTpotomak2cestica94potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94 -> getDaughters()[0]; // 1. potomak "cestice" 94
													EVENT::MCParticle* prviTpotomak2cestica94potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94 -> getDaughters()[1]; // 2. potomak "cestice" 94

													if ((prviTpotomak2cestica94potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1 -> getGeneratorStatus() == 1))
													{// 2.5.3
														cout << "Nadjen je 1. potomak cestice 94 " << prviTpotomak2cestica94potomak1 -> getPDG() << " cestica je finalna." << endl;
														pointer_to_MCpion1 = prviTpotomak2cestica94potomak1;
													}// 2.5.3

													if ((prviTpotomak2cestica94potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak2 -> getGeneratorStatus() == 1))
													{// 2.5.3
														cout << "Nadjen je 2. potomak cestice 94 " << prviTpotomak2cestica94potomak2 -> getPDG() << " cestica je finalna." << endl;
														pointer_to_MCpion1 = prviTpotomak2cestica94potomak2;
													}// 2.5.3

													if(prviTpotomak2cestica94potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak "cestice" 94 nije konacan.
													{// 2.5.4 Uslov za potomke 1. potomka "cestice" 94
														EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka "cestice" 94
														EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka "cestice" 94

														if ((prviTpotomak2cestica94potomak1potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak1 -> getGeneratorStatus() == 1))
														{// 2.5.5 Uslov po potomcima 1. potomka "cestice" 94
															cout << "Nadjen je 1. potomak 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak1;
														}// 2.5.5

														if ((prviTpotomak2cestica94potomak1potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak2 -> getGeneratorStatus() == 1))
														{// 2.5.5 Uslov po potomcima 1. potomka "cestice" 94
															cout << "Nadjen je 2. potomak 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak2;
														}// 2.5.5

														if(prviTpotomak2cestica94potomak1potomak1 -> getGeneratorStatus() != 1)
														{// 2.5.6 Uslov po potomcima 1. potomka "cestice" 94
															EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka "cestice" 94
															EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka "cestice" 94

															if ((prviTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po potomcima 1. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 1. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak1potomak1;
															}// 2.5.7

															if ((prviTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po potomcima 1. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 1. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak1potomak2;
															}// 2.5.7

															if (prviTpotomak2cestica94potomak1potomak1potomak1 -> getGeneratorStatus() != 1)
															{// 2.5.8
																EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak1potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka 1. potomka "cestice" 94

																if ((prviTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 1. potomka 1. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak1potomak1potomak1;
																}
															}// 2.5.8
														}// 2.5.6 Kraj uslova po potomcima 1. potomka "cestice" 94

														if(prviTpotomak2cestica94potomak1potomak2 -> getGeneratorStatus() != 1)
														{// 2.5.6 Uslov po potomcima 2. potomka "cestice" 94
															EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 1. potomka "cestice" 94
															EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 1. potomka "cestice" 94

															if ((prviTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak2potomak1 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po cerkama cerke prve cerke 94
																cout << "Nadjen je 1. potomak 2. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak2potomak1;
															}// 2.5.7

															if ((prviTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak2potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po cerkama cerke prve cerke 94
																cout << "Nadjen je 2. potomak 2. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak2potomak2;
															}// 2.5.7

															if (prviTpotomak2cestica94potomak1potomak2potomak1 -> getGeneratorStatus() != 1)
															{// 2.5.8
																EVENT::MCParticle* prviTpotomak2cestica94potomak1potomak1potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak1potomak2potomak1 -> getDaughters()[0]; // 1. potomak 2. potomka 1. potomka s1. potomka "cestice" 94

																if ((prviTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 2. potomka 1. potomka 1. potomka cestice 94 " << prviTpotomak2cestica94potomak1potomak1potomak2potomak1 ->getPDG() << " cestica je finalna"<< endl;
																	pointer_to_MCpion1 = prviTpotomak2cestica94potomak1potomak1potomak2potomak1;
																}
															}// 2.5.8
														}// 2.5.6 Kraj uslova po potomcima 2. potomka "cestice" 94
													}// 2.5.4 Kraj uslova za potomke 1. potomka "cestice" 94

													if(prviTpotomak2cestica94potomak2 -> getGeneratorStatus() != 1)  // Ako 2. potomak "cestice" 94 nije konacan.
													{// 2.6.4 Uslov za potomke 2. potomka "cestice" 94
														EVENT::MCParticle* prviTpotomak2cestica94potomak2potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak2 -> getDaughters()[0]; // cerka druge cerke od 94
														EVENT::MCParticle* prviTpotomak2cestica94potomak2potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak2 -> getDaughters()[1]; // cerka druge cerke od 94

														if ((prviTpotomak2cestica94potomak2potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak2potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak2potomak1 -> getGeneratorStatus() == 1))
														{// 2.6.5 Uslov po potomcima druge cerke 94
															cout << "Nadjen je 1. potomak 2. potomka cestice 94 " << prviTpotomak2cestica94potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak2cestica94potomak2potomak1;
														} // 2.6.5

														if ((prviTpotomak2cestica94potomak2potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak2potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak2potomak2 -> getGeneratorStatus() == 1))
														{// 2.6.5 Uslov po potomcima 2. potomka "cestice" 94
															cout << "Nadjen je 2. potomak 2. potomka cestice 94 " << prviTpotomak2cestica94potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion1 = prviTpotomak2cestica94potomak2potomak2;
														}// 2.6.5

														if(prviTpotomak2cestica94potomak2potomak1 -> getGeneratorStatus() != 1)
														{// 2.6.6
															EVENT::MCParticle* prviTpotomak2cestica94potomak2potomak1potomak1 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak2potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 2. potomka "cestice" 94
															EVENT::MCParticle* prviTpotomak2cestica94potomak2potomak1potomak2 = (EVENT::MCParticle*) prviTpotomak2cestica94potomak2potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 2. potomka "cestice" 94

															if ((prviTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() == -211 || prviTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() == 16) && (prviTpotomak2cestica94potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
															{// 2.6.7 Petlja po potomcima 1. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 1. potomka 2. potomka cestice 94 " << prviTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak2potomak1potomak1;
															}// 2.6.7

															if ((prviTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() == -211 || prviTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() == 16) && (prviTpotomak2cestica94potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.6.7 Petlja po potomcima 1. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 1. potomka 2. potomka cestice 94 " << prviTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion1 = prviTpotomak2cestica94potomak2potomak1potomak2;
															}// 2.6.7
														}// 2.6.6
													}// 2.6.4 Kraj uslova za 2. potomka "cestice" 94
												}// 2.5.2
											}// 2.5.1
										}// 2.5 Kraj uslova za 2. potomka 1. tauona
									}// Kraj uslova ako 1. potomak 1. tauona nema potomaka
								}// 2.4 Kraj uslova po 1. tauonu
							}// Kraj uslova kada imamo 2 piona i 2 neutrina
						} // 2.3
					}// 2.2
				}// 2.1
			}// 1.2	Kraj petlje po tauonima

		}// 1.1 Kraj uslova za Higsa

	} // Kraj FOR petlje za broj elemenata

	return pointer_to_MCpion1;

} // Kraj funkcije getMCPion1

// ---------------------------------------------------------------------------------------------------------------------

EVENT::MCParticle* getMCPion2(EVENT::LCEvent* evt)
{
	EVENT::MCParticle* pointer_to_MCpion2 = 0;

	std::vector<std::string> colNames = *evt -> getCollectionNames();

	IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*)evt -> getCollection("MCParticlesSkimmed");

	Int_t N_piona = 0;
	Int_t N_neutrina = 0;
	Int_t N_event = 0;

	for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
	{// 1
		IMPL::MCParticleImpl* generated = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

		if (abs(generated -> getPDG() == 25))
		{// 1.1
			cout << "Nadjen je Higgs \n";

			if (abs (generated -> getDaughters()[0] -> getPDG() == 25))
			{// 1.2 traze se tauoni
				const EVENT::MCParticleVec & Tauoni = generated -> getDaughters()[0] -> getDaughters();

				if (Tauoni.size() == 2)
				{// 2.1
					EVENT::MCParticle* prviT = (EVENT::MCParticle*) Tauoni[0];
					EVENT::MCParticle* drugiT = (EVENT::MCParticle*) Tauoni[1];

					if (prviT -> getPDG() == 15 && drugiT -> getPDG() == -15)
					{// 2.2 petlja po T1 i T2
					cout << "Nadjeni su tauoni" << endl;

						if (prviT -> getDaughters().size() >= 2 && drugiT -> getDaughters().size() >= 2)
						{// 2.3
							EVENT::MCParticle* prviTpotomak1 = (EVENT::MCParticle*) prviT -> getDaughters()[0]; // 1. potomak T1
							EVENT::MCParticle* prviTpotomak2 = (EVENT::MCParticle*) prviT -> getDaughters()[1]; // 2. potomak T1
							EVENT::MCParticle* drugiTpotomak1 = (EVENT::MCParticle*) drugiT -> getDaughters()[0]; // 1. potomak T2
							EVENT::MCParticle* drugiTpotomak2 = (EVENT::MCParticle*) drugiT -> getDaughters()[1]; // 2. potomak T2

							if (abs (prviTpotomak1 -> getPDG()) == 211) N_piona++;
							if (abs (prviTpotomak2 -> getPDG()) == 16) N_neutrina++;
							if (abs (drugiTpotomak1 -> getPDG()) == 211) N_piona++;
							if (abs (drugiTpotomak2 -> getPDG()) == 16) N_neutrina++;

							if (N_piona == 2 && N_neutrina == 2)
							{// Uslov kada imamo 2 piona i 2 neutrina
								N_event++;

								cout << "Nadjen Tauon-Pion dogadjaj broj " << N_event << endl;

								if ((drugiTpotomak1 -> getPDG() == 211 && drugiTpotomak2 -> getPDG() == -16) || (drugiTpotomak1 -> getPDG() == -16 && drugiTpotomak2 -> getPDG() == 211))
								{// 2.4 Uslov po potomcima 2. tauona

									if (drugiTpotomak1 -> getGeneratorStatus() == 1)
									{
										cout << "Nadjen je 1. potomak 2. tauona " << drugiTpotomak1 -> getPDG() << " cestica je finalna." << endl;
										pointer_to_MCpion2 = drugiTpotomak1;
									}

									if (drugiTpotomak2 -> getGeneratorStatus() == 1)
									{
										cout << "Nadjen je 2. potomak 2. tauona " << drugiTpotomak2 -> getPDG() << " cestica je finalna." << endl;
										pointer_to_MCpion2 = drugiTpotomak2;
									}

									if ((drugiTpotomak1 -> getGeneratorStatus() != 1) && (drugiTpotomak2 -> getGeneratorStatus() != 1))
									{// Uslov za 1. i 2. potomak 2. tauona

										if (drugiTpotomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak nije konacan.
										{// 2.5 Uslov za 1. potomka 2. tauona

											for (Int_t i = 0;  (int)i < drugiTpotomak1 -> getDaughters().size(); i++)
											{
												EVENT::MCParticle* drugiTpotomak1cestica94 = (EVENT::MCParticle*) drugiTpotomak1 -> getDaughters()[0]; // Ako je 1. potomak 2. tauona "cestice" 94

												if (abs(drugiTpotomak1cestica94 -> getPDG()) == 94)
												{// 2.5.1 Uslov za "cesticu" 94

													if (drugiTpotomak1cestica94 -> getGeneratorStatus() != 1)
													{// 2.5.2 Uslov po potomcima "cestice" 94
														EVENT::MCParticle* drugiTpotomak1cestica94potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94 -> getDaughters()[0]; // 1. potomak "cestice" 94
														EVENT::MCParticle* drugiTpotomak1cestica94potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94 -> getDaughters()[1]; // 2. potomak "cestice" 94

														if ((drugiTpotomak1cestica94potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1 -> getGeneratorStatus() == 1))
														{// 2.5.3
															cout << "Nadjen je 1. potomak cestice 94 " << drugiTpotomak1cestica94potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1;
														}// 2.5.3

														if ((drugiTpotomak1cestica94potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2 -> getGeneratorStatus() == 1))
														{// 2.5.3
															cout << "Nadjen je 2. potomak cestice 94 " << drugiTpotomak1cestica94potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2;
														}// 2.5.3

														if(drugiTpotomak1cestica94potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak "cestice" 94 nije konacan
														{// 2.5.4 Uslov za 1. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak1cestica94potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka "cerke" 94
															EVENT::MCParticle* drugiTpotomak1cestica94potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka "cerke" 94

															if ((drugiTpotomak1cestica94potomak1potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1potomak1 -> getGeneratorStatus() == 1))
															{// 2.5.5 Uslov po potomcima 1. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 1. potomka cestice 94 " << drugiTpotomak1cestica94potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1potomak1;
															}// 2.5.5

															if ((drugiTpotomak1cestica94potomak1potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.5 Uslov po potomcima 1. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 1. potomka cestice 94 " << drugiTpotomak1cestica94potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1potomak2;
															}// 2.5.5

															if(drugiTpotomak1cestica94potomak1potomak1 -> getGeneratorStatus() != 1)
															{// 2.5.6
																EVENT::MCParticle* drugiTpotomak1cestica94potomak1potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak1cestica94potomak1potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka "cestice" 94

																if ((drugiTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																{// 2.5.7 Uslov po potomcima 1. potomka 1. potomka "cestice" 94
																	cout << "Nadjen je 1. potomak 1. potomka 1. potomka cestice 94 " << drugiTpotomak1cestica94potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1potomak1potomak1;
																}// 2.5.7

																if ((drugiTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
																{// 2.5.7 Uslov po potomcima 1. potomka 1. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 1. potomka 1. potomka cestice 94 " << drugiTpotomak1cestica94potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1potomak1potomak2;
																}// 2.5.7

																if (drugiTpotomak1cestica94potomak1potomak1potomak1 -> getGeneratorStatus() != 1)
																{// 2.5.8
																	EVENT::MCParticle* drugiTpotomak1cestica94potomak1potomak1potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak1potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka "cestice" 94

																	if ((drugiTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																	{
																		cout << "Nadjen je 1. potomak 1. potomka 1. potomka cestice 94 " << drugiTpotomak1cestica94potomak1potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																		pointer_to_MCpion2 = drugiTpotomak1cestica94potomak1potomak1potomak1potomak1;
																	}
																}// 2.5.8
															}// 2.5.6
														}// 2.5.4 Kraj uslova za 1. potomka "cestice" 94

														if(drugiTpotomak1cestica94potomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak "cestice" 94 nije konacan.
														{// 2.6.4 Uslov za 2. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka "cestice" 94

															if ((drugiTpotomak1cestica94potomak2potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak1 -> getGeneratorStatus() == 1))
															{// 2.6.5 Uslov po potomcima 2. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak1;
															}// 2.6.5

															if ((drugiTpotomak1cestica94potomak2potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak2 -> getGeneratorStatus() == 1))
															{// 2.6.5 Uslov po potomcima 2. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak2;
															} // 2.6.5

															if(drugiTpotomak1cestica94potomak2potomak1 -> getGeneratorStatus() != 1)
															{// 2.6.6
																EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 2. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 2. potomka "cestice" 94

																if ((drugiTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
																{// 2.6.7 Uslov po potomcima 1. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 1. potomak 1. potomka 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak1potomak1;
																}// 2.6.7

																if ((drugiTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
																{// 2.6.7 Uslov po potomcima 1. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 1. potomka 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak1potomak2;
																}// 2.6.7
															}// 2.6.6

															if(drugiTpotomak1cestica94potomak2potomak2 -> getGeneratorStatus() != 1)
															{// 2.6.6
																EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 2. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak1cestica94potomak2potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak1cestica94potomak2potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 2. potomka "cestice" 94

																if ((drugiTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak2potomak1 -> getGeneratorStatus() == 1))
																{// 2.6.7 Uslov po potomcima 2. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 1. potomak 2. potomka 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak2potomak1;
																}// 2.6.7

																if ((drugiTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() == 211 || drugiTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak1cestica94potomak2potomak2potomak2 -> getGeneratorStatus() == 1))
																{// 2.6.7 Uslov po potomcima 2. potomka 2. potomka "cestice" 94
																	cout << "Nadjen je 2. potomak 2. potomka 2. potomka cestice 94 " << drugiTpotomak1cestica94potomak2potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak1cestica94potomak2potomak2potomak2;
																}// 2.6.7
															}// 2.6.6
														}// 2.6.4 Kraj uslova za 2. potomka "cestice" 94
													}// 2.5.2
												}// 2.5.1 Kraj uslova za "cesticu" 94
											} // Kraj FOR petlje po potomcima 1. potomka 2. tauona
										}// 2.5 Kraj uslova za 1. potomka 2. tauona
									}// Kraj uslova za generatorski status 1. i 2. potomka 2. tauona

									///////////////////// ako 1.potomak 2. tauona nema potomke ////////////////////
									if (drugiTpotomak1 -> getDaughters().size() == 0)
									{
										if (drugiTpotomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak 2. tauona nije konacna cestica.
										{//2.5 Uslov za 2. potomka 2. tauona

											if (drugiTpotomak2 -> getDaughters().size() == 0) continue;

											EVENT::MCParticle* drugiTpotomak2cestica94 = (EVENT::MCParticle*)drugiTpotomak2 ->getDaughters()[0];//prva cerka prve cerke z2, 94

											if (abs(drugiTpotomak2cestica94 -> getPDG()) == 94)
											{// 2.5.1

												if (drugiTpotomak2cestica94 -> getGeneratorStatus() != 1)
												{// 2.5.2 Uslov po potomcima "cestice" 94
													EVENT::MCParticle* drugiTpotomak2cestica94potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94 -> getDaughters()[0]; // 1. potomak "cestice" 94
													EVENT::MCParticle* drugiTpotomak2cestica94potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94 -> getDaughters()[1]; // 2. potomak "cestice" 94

													if ((drugiTpotomak2cestica94potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1 -> getGeneratorStatus() == 1))
													{// 2.5.3
														cout << "Nadjen je 1. potomak cestice 94 2. tauona " << drugiTpotomak2cestica94potomak1 -> getPDG() << " cestica je finalna." << endl;
														pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1;
													}// 2.5.3

													if ((drugiTpotomak2cestica94potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2 -> getGeneratorStatus() == 1))
													{// 2.5.3
														cout << "Nadjen je 2. potomak cestice 94 2. tauona " << drugiTpotomak2cestica94potomak2 -> getPDG() << " cestica je finalna." << endl;
														pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2;
													}// 2.5.3

													if(drugiTpotomak2cestica94potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak "cestice" 94 nije konacan.
													{// 2.5.4 Uslov za 1. potomka "cestice" 94
														EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka "cestice" 94
														EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka "cestice" 94

														if ((drugiTpotomak2cestica94potomak1potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1 -> getGeneratorStatus() == 1))
														{// 2.5.5 Uslov po potomcima 1. potomka "cestice" 94
															cout << "Nadjen je 1. potomak 1. potomka cestice 94 2. tauona " << drugiTpotomak2cestica94potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1;
														}// 2.5.5

														if ((drugiTpotomak2cestica94potomak1potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak2 -> getGeneratorStatus() == 1))
														{// 2.5.5 Uslov po potomcima 2. potomka "cestice" 94
															cout << "Nadjen je 2. potomak 1. potomka cestice 94 2. tauona " << drugiTpotomak2cestica94potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak2;
														}// 2.5.5

														if(drugiTpotomak2cestica94potomak1potomak1 -> getGeneratorStatus() != 1)
														{//2.5.6
															EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka "cestice" 94

															if ((drugiTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
															{//2.5.7 Uslov po potomcima 1. potomka 1. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 1. potomka 1. potomka cestice 94 2. tauona " << drugiTpotomak2cestica94potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak1;
															}// 2.5.7

															if ((drugiTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po potomcima 1. potomka 1. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 1. potomka 1. potomka cestice 94 2. tauona " << drugiTpotomak2cestica94potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak2;
															}// 2.5.7

															if (drugiTpotomak2cestica94potomak1potomak1potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak 1. potomka 1. potomka nije konacan
															{// 2.5.8
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 1. potomka 1. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 1. potomka 1. potomka "cestice" 94

																if ((drugiTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 1. potomka 1. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak1potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak1potomak1;
																}

																if ((drugiTpotomak2cestica94potomak1potomak1potomak1potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak1potomak2 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 1. potomka 1. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak1potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak1potomak2;
																}
															}// 2.5.8

															if (drugiTpotomak2cestica94potomak1potomak1potomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak 1. potomka 1. potomka nije konacan
															{// 2.5.8
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 1. potomka 1. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak1potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak1potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 1. potomka 1. potomka "cestice" 94

																if ((drugiTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 2. potomka 1. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak1potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak2potomak1;
																}

																if ((drugiTpotomak2cestica94potomak1potomak1potomak2potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak1potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak1potomak2potomak2 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 2. potomak 2. potomka 1. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak1potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak1potomak2potomak2;
																}
															}// 2.5.8

														}//2.5.6 drugiTpotomak2cestica94potomak1potomak1

														if(drugiTpotomak2cestica94potomak1potomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak 1. potomka "cestice" 94 nije konacan.
														{// 2.5.6
															EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 1. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 1. potomka "cestice" 94

															if ((drugiTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak2potomak1 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov po potomcima 2. potomka 1. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 2. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak2potomak1;
															}// 2.5.7

															if ((drugiTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak2potomak2 -> getGeneratorStatus() == 1))
															{// 2.5.7 Uslov potomcima 2. potomka 1. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 2. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak2potomak2;
															}// 2.5.7

															if (drugiTpotomak2cestica94potomak1potomak2potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak 2. potomka 1. potomka "cestice" 94 nije konacan
															{// 2.5.8
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak2potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak2potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 2. potomka 1. potomka "cestice" 94
																EVENT::MCParticle* drugiTpotomak2cestica94potomak1potomak2potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak1potomak2potomak1 -> getDaughters()[1]; // 1. potomak 1. potomka 2. potomka 1. potomka "cestice" 94

																if ((drugiTpotomak2cestica94potomak1potomak2potomak1potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak2potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 1. potomak 1. potomka 2. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak2potomak1potomak1;
																}

																if ((drugiTpotomak2cestica94potomak1potomak2potomak1potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak1potomak2potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak1potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
																{
																	cout << "Nadjen je 2. potomak 1. potomka 2. potomka 1. potomka cestice 94 " << drugiTpotomak2cestica94potomak1potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																	pointer_to_MCpion2 = drugiTpotomak2cestica94potomak1potomak2potomak1potomak2;
																}
															}// 2.5.8

														}// 2.5.6 Kraj uslova ako 2. potomak 1. potomka "cestice" 94 nije konacan.
													}// 2.5.4 Kraj uslova ako 1. potomak "cestice" 94 nije konacan.

													if(drugiTpotomak2cestica94potomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak "cestice" 94 nije konacan.
													{// 2.6.4 Uslov za 2. potomka "cestice" 94
														EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka "cestice" 94
														EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka "cestice" 94

														if ((drugiTpotomak2cestica94potomak2potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak1 -> getGeneratorStatus() == 1))
														{// 2.6.5 Uslov po potomcima 2. potomka "cestice" 94
															cout << "Nadjen je 1. potomak 2. potomka cestice 94 " << drugiTpotomak2cestica94potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak1;
														}// 2.6.5

														if ((drugiTpotomak2cestica94potomak2potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak2 -> getGeneratorStatus() == 1))
														{// 2.6.5 Uslov po potomcima 2. potomka "cestice" 94
															cout << "Nadjen je 2. potomak 2. potomka cestice 94 " << drugiTpotomak2cestica94potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
															pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak2;
														}// 2.6.5

														if(drugiTpotomak2cestica94potomak2potomak1 -> getGeneratorStatus() != 1) // Ako 1. potomak 2. potomka "cestice" 94 nije konacan.
														{// 2.6.6
															EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak1potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2potomak1 -> getDaughters()[0]; // 1. potomak 1. potomka 2. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak1potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2potomak1 -> getDaughters()[1]; // 2. potomak 1. potomka 2. potomka "cestice" 94

															if ((drugiTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak1potomak1 -> getGeneratorStatus() == 1))
															{// 2.6.7 Uslov po potomcima 1. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 1. potomka 2. potomka cestice 94 " << drugiTpotomak2cestica94potomak2potomak1potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak1potomak1;
															}// 2.6.7

															if ((drugiTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak1potomak2 -> getGeneratorStatus() == 1))
															{// 2.6.7 Uslov po potomcima 1. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 1. potomka 2. potomka cestice 94" << drugiTpotomak2cestica94potomak2potomak1potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak1potomak2;
															}// 2.6.7
														}// 2.6.6

														if(drugiTpotomak2cestica94potomak2potomak2 -> getGeneratorStatus() != 1) // Ako 2. potomak 2. potomka "cestice" 94 nije konacan.
														{// 2.6.6
															EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak2potomak1 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2potomak2 -> getDaughters()[0]; // 1. potomak 2. potomka 2. potomka "cestice" 94
															EVENT::MCParticle* drugiTpotomak2cestica94potomak2potomak2potomak2 = (EVENT::MCParticle*) drugiTpotomak2cestica94potomak2potomak2 -> getDaughters()[1]; // 2. potomak 2. potomka 2. potomka "cestice" 94

															if ((drugiTpotomak2cestica94potomak2potomak2potomak1 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak2potomak1 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak2potomak1 -> getGeneratorStatus() == 1))
															{// 2.6.7 Uslov po potomcima 2. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 1. potomak 2. potomka 2. potomka cestice 94 " << drugiTpotomak2cestica94potomak2potomak2potomak1 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak2potomak1;
															}// 2.6.7

															if ((drugiTpotomak2cestica94potomak2potomak2potomak2 -> getPDG() == 211 || drugiTpotomak2cestica94potomak2potomak2potomak2 -> getPDG() == -16) && (drugiTpotomak2cestica94potomak2potomak2potomak2 -> getGeneratorStatus() == 1))
															{// 2.6.7 Uslov po potomcima 2. potomka 2. potomka "cestice" 94
																cout << "Nadjen je 2. potomak 2. potomka 2. potomka cestice 94" << drugiTpotomak2cestica94potomak2potomak2potomak2 -> getPDG() << " cestica je finalna." << endl;
																pointer_to_MCpion2 = drugiTpotomak2cestica94potomak2potomak2potomak2;
															}// 2.6.7
														}// 2.6.6
													}// 2.6.4 Kraj uslova za 2. potomka "cestice" 94
												}// 2.5.2 Kraj uslova po potomcima "cestice" 94
											}// 2.5.1
										}// 2.5 Kraj uslova za 2. potomka 2. tauona
									}// Kraj uslova za potomke 1. potomka 2. tauona
								}// 2.4 Kraj uslova po potomcima 2. tauona
							}// Kraj uslova kada imamo 2 piona i 2 neutrina
						} // 2.3
					}// 2.2
				}// 2.1
			}// 1.2	Kraj petlje po tauonima

		}// 1.1 Kraj uslova za Higsa

	} // 1 Kraj FOR petlje

	return pointer_to_MCpion2;

} // Kraj funkcije getMCPion2

// --------------------------------------------------------------------------------------------------------------------

Int_t slcio2appTree(UInt_t nFirstJob, UInt_t nLastJob, const char * fn, const char * rfn)
{
	#ifdef __CINT__
		gSystem -> Load("${LCIO}/lib/liblcio.so");
		gSystem -> Load("${LCIO}/lib/liblcioDict.so");
	#endif

	// TFile rootFile("Polarimetar_Jinx_MC.root","RECREATE", "Polarimetar_Jinx_MC", 1);
	TTree JinxTree ("JinxTree", "Generator particle tree");

	Float_t Pt_Tauon_minus;
	Float_t Pt_Tauon_plus;
	Float_t Pt_Pion_minus;
	Float_t Pt_Pion_plus;
//	Float_t Polarimetar_6_minus;
//	Float_t Polarimetar_6_plus;
//	Float_t Delta_Fi;

	JinxTree.Branch("Pt_Tauon_minus", &Pt_Tauon_minus, "Pt_Tauon_minus");
	JinxTree.Branch("Pt_Tauon_plus", &Pt_Tauon_plus, "Pt_Tauon_plus");
	JinxTree.Branch("Pt_Pion_minus", &Pt_Pion_minus, "Pt_Pion_minus");
	JinxTree.Branch("Pt_Pion_plus", &Pt_Pion_plus, "Pt_Pion_plus");
//	JinxTree.Branch("Polarimetar_6_minus", &Polarimetar_6_minus, "Polarimetar_6_minus");
//	JinxTree.Branch("Polarimetar_6_plus", &Polarimetar_6_plus, "Polarimetar_6_plus");
//	JinxTree.Branch("Delta_Fi", &Delta_Fi, "Delta_Fi");

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance() -> createLCReader();
	TString fName = fn;
	stringstream fNameStream;

	// Int_t N_event_ukupno = 0;
	// Int_t N_isolep = 0;
	// Int_t N_signal = 0;
	// const Double_t mH_theory = 125.0;
	// bool right_evt = false;

	// N_event = 0;

	// Petlja koja iščitava .slcio fajlove
	for(UInt_t iJob = nFirstJob; iJob <= nLastJob; iJob++)
	{
		cout << "Opening " << Form("%s%i.slcio", fName.Data(), iJob);

		try
		{
			lcReader -> open(Form("%s%i.slcio", fName.Data(), iJob));
		}

		catch(lcio::IOException &ex)
		{
			cout << ". Ne mere se otvorit.\n";
			continue;
		}

		cout << ". Učitavanje.\n";

		Int_t N_event = 0;

		// Petlja po dogadjajima
		EVENT::LCEvent* evt = 0;

		while( (evt = lcReader -> readNextEvent()) != 0)
		{
			N_event++;

			std::vector<std::string> colNames = *evt -> getCollectionNames();

			IMPL::LCCollectionVec* mcParticles = (IMPL::LCCollectionVec*) evt -> getCollection("MCParticlesSkimmed");
			// IMPL::LCCollectionVec* pfos = (IMPL::LCCollectionVec*) evt -> getCollection("PandoraPFANewPFOs");
			// IMPL::LCCollectionVec* colJet = (IMPL::LCCollectionVec*) evt -> getCollection("twoRefJetsZep");
			// IMPL::LCCollectionVec* jets2 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_2");
			// IMPL::LCCollectionVec* jets4 = (IMPL::LCCollectionVec*) evt -> getCollection("FJ_Jets_4");
			// IMPL::LCCollectionVec* isolep = (IMPL::LCCollectionVec*) evt -> getCollection("Isolep_Selected");

			vector <TLorentzVector> Tauon_minus;
			vector <TLorentzVector> Tauon_plus;
			// vector <TLorentzVector> Pion_minus_1;
			// vector <TLorentzVector> Pion_minus_2;
			// vector <TLorentzVector> Pion_plus_1;
			// vector <TLorentzVector> Pion_plus_2;

			Int_t N_tauona_plus = 0;
			Int_t N_tauona_minus = 0;
			Int_t N_piona_plus = 0;
			Int_t N_piona_minus = 0;

			Int_t N_tauona_plus_ukupno = 0;
			Int_t N_tauona_minus_ukupno = 0;

			cout << "------------" << endl;
			cout << N_event << ". događaj" << endl;
			cout << "------------" << endl;

			for (Int_t i = 0; i < mcParticles -> getNumberOfElements(); i++)
			{
				IMPL::MCParticleImpl* mcParticle = (IMPL::MCParticleImpl*) mcParticles -> getElementAt(i);

				if ( mcParticle -> getPDG() == 25 && mcParticle -> getDaughters()[0] -> getPDG() == 15 && mcParticle -> getDaughters()[1] -> getPDG() == -15 )
				{
//				const EVENT::MCParticleVec & parent = mcParticle -> getParents();
//				const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

					TLorentzVector tau1( TVector3(mcParticle -> getDaughters()[0] -> getMomentum()), mcParticle -> getDaughters()[0] -> getEnergy() );
					TLorentzVector tau2( TVector3(mcParticle -> getDaughters()[1] -> getMomentum()), mcParticle -> getDaughters()[1] -> getEnergy() );

//				TLorentzVector tau1( mcParticle -> getDaughters()[0] -> getSpin() );
//				TLorentzVector tau2( mcParticle -> getDaughters()[1] -> getSpin() );

					vector <TLorentzVector> Pion_minus;
					vector <TLorentzVector> Pion_plus;

					if ( mcParticle -> getPDG() == 15 )
					{
						if ( mcParticle -> getPDG() == 15 && mcParticle -> getDaughters()[0] -> getPDG() == -211 && mcParticle -> getDaughters()[1] -> getPDG() == 16 )
						{
//						const EVENT::MCParticleVec & parent = mcParticle -> getParents();
//						const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

							EVENT::MCParticle* pion1 = getMCpion1(evt);

							cout << "Događaj sa pionom minus " << endl;

							N_piona_minus++;

							TLorentzVector pion1a( TVector3(mcParticle -> getDaughters()[0] -> getMomentum()), mcParticle -> getDaughters()[0] -> getEnergy() );

							Pion_minus.push_back(pion1a);
						}

						Pt_Pion_minus = Pion_minus.Pt();

					}

					if ( mcParticle -> getPDG() == -15 )
					{
						if ( mcParticle -> getPDG() == -15 && mcParticle -> getDaughters()[0] -> getPDG() == 211 && mcParticle -> getDaughters()[1] -> getPDG() == -16 )
						{
//						const EVENT::MCParticleVec & parent = mcParticle -> getParents();
//						const EVENT::MCParticleVec & daughter = mcParticle -> getDaughters();

							EVENT::MCParticle* pion2 = getMCpion2(evt);

							cout << "Događaj sa pionom plus " << endl;

							N_piona_plus++;

							TLorentzVector pion1c( TVector3(mcParticle -> getDaughters()[0] -> getMomentum()), mcParticle -> getDaughters()[0] -> getEnergy() );

							Pion_plus.push_back(pion1c);
						}

						Pt_Pion_plus = Pion_plus.Pt();

					}

					Tauon_minus.push_back(tau1);
					Tauon_plus.push_back(tau2);

					N_tauona_plus++;
					N_tauona_minus++;

				} // Kraj IF uslova za tauone

				N_tauona_plus_ukupno++;
				N_tauona_minus_ukupno++;

			} // Kraj FOR petlje po broju elemenata

			Pt_Tauon_minus = Tauon_minus[0].Pt();
			Pt_Tauon_plus = Tauon_plus[0].Pt();

			cout << "Broj tauona minus " << N_tauona_minus << endl;
			cout << "Broj tauona plus " << N_tauona_plus << endl;
			cout << "Broj tauona minus ukupno " << N_tauona_minus_ukupno << endl;
			cout << "Broj tauona plus ukupno " << N_tauona_plus_ukupno << endl;
			cout << endl;

			cout << "Broj piona minus " << N_piona_minus << endl;
			cout << "Broj piona plus " << N_piona_plus << endl;
			cout << endl;

		} // Kraj WHILE petlje

		JinxTree.Fill();

		lcReader -> close();

	} // Kraj FOR petlje koja iščitava .slcio fajlove

//	JinxTree.Fill();

	TString tfName(rfn);
	if(!tfName.EndsWith(".root")) tfName.Append(".root");
	TFile rootFile(tfName.Data(),"RECREATE");

	JinxTree.Write();
	// evtTree.Write();
	rootFile.Write();
	rootFile.Close(); //

	return 0;

} // Kraj funkcije

// ------------------------------------------------------------------

Int_t main(int argc, char* argv[])
{
	Int_t iarg = 1;

	UInt_t nFirstJob = 1;

	if(argc > iarg) {nFirstJob = atoi(argv[iarg]); iarg++;}

	UInt_t nLastJob = 50;

	if(argc > iarg) {nLastJob  = atoi(argv[iarg]); iarg++;}

	TString fName = "test_";
	if(argc > iarg) fName = argv[iarg]; //iarg++;

	TString rfName = "Polarimetar_Jinx_MC.root";
	if(argc > iarg) rfName = argv[0];

	return slcio2appTree(nFirstJob, nLastJob, fName.Data(), rfName.Data());
}
