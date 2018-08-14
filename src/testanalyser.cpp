#include <vector>
#include "interface/testanalyser.h"
#include "../interface/helpers.h"
#include "../interface/ttbar_solver.h"
#include "../interface/pzCalculator.h"
#include "TMath.h"
using namespace std;

double M_mu =  0.10566;
double M_e = 0.511e-3;
double M_T = 172.5;

void testanalyser::analyze(size_t childid /* this info can be used for printouts */){

	d_ana::dBranchHandler<Electron>		elecs(tree(),"Electron");
	d_ana::dBranchHandler<HepMCEvent>	event(tree(),"Event");
	d_ana::dBranchHandler<Jet>		genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>		jet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>		muontight(tree(),"Muon");
	d_ana::dBranchHandler<MissingET>	met(tree(),"MissingET");
	d_ana::dBranchHandler<Jet>		ParticleFlowJet04(tree(),"ParticleFlowJet04");

	TH1* neutrinoMomentumHisto = addPlot(new TH1D("neutrinoMomentumHisto", "nu histo",100,0,200), "pT","N");
	TH1* histoSL=addPlot(new TH1D("tquark_SL","t-quark m_{inv} mass in semileptonic decay",100,0,500),"M [GeV]","N");
	TH1* topMHistoNonBoosted=addPlot(new TH1D("topMHistoBoostNonBoosted","t-quark m_{inv} mass",100,0,400),"M [GeV]","N");
	TH1* tbarMHistoNonBoosted=addPlot(new TH1D("tbarMHistoBoostNonBoosted","anti t-quark m_{inv} mass",100,0,400),"M [GeV]","N");
	TH1* massDiffHistoNonBoosted = addPlot(new TH1D("massDiffHistoNonBoosted", "m_t - m_{#bar{t}}", 50, -50, 50), "#Delta m", "N_{events}");
	TH1* chi2Histo= addPlot(new TH1D("chi2Histo", "#Chi^2", 100, 0, 100), "#Chi^2", "N_{events}");
	TH1* chi2VsMassDiffHisto= addPlot(new TH1D(" chi2VsMassDiffHisto", "#Chi^2 vs mass difference", 50, -3, 3), "#Delta m", "#Chi^2");
	TH1* chi2VsMassHisto= addPlot(new TH1D(" chi2VsMassHisto", "#Chi^2 vs top mass", 50, 0, 300), "m_t", "#Chi^2");
//	TH1* chi2Histo2D= addPlot(new TH2D("chi2Histo2D", "#Chi^2 for m_t and m_{#bar{t}}", 100, 0, 200, 100, 0,200), "top", "antitop", "Chi^2");
	TH1* hadronicTopMassHisto = addPlot(new TH1D("hadronicTopMassHisto", "Top mass from hadronic decay", 50, 0, 400), "m_t (Gev)", "N");
	TH1* hadronicTbarMassHisto = addPlot(new TH1D("hadronicTbarMassHisto", "Anti top mass from hadronic decay", 50, 0, 400), "m_t (Gev)", "N");
	
	size_t nevents=tree()->entries();

	if(isTestMode())
		nevents/=100;

	for(size_t eventno=0; eventno<nevents; eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);

		//if (ParticleFlowJet04.size()<4) continue;
		bool acceptableEvent = true;

		// an even must contain at least 2 jets, 2 of which must be b-tagged 
		size_t bjets=0;
		size_t ljets=0;
		for(size_t i=0; i<jet.size(); i++){
			if(jet.at(i)->BTag){
				bjets++;

			}else{
				ljets++;
			}
		}

		if(bjets<2 || ljets<2) continue;

		// remaining jets are required to satisfy certain conditions
		for(size_t i=0; i<jet.size(); i++){
			if(!acceptLJet(jet.at(i))
			&& !acceptBJet(jet.at(i))){
				acceptableEvent = false;
			}
		}

		if (!acceptableEvent) continue;

		//define the lower bound for lepton transverse momentum according to the paper
		double muonPtLimit = 20;
		double electronEtLimit = 25;
		
		//determine the number of leptons with PT>ptLimit
		size_t nrOfHighPtLeps=0;
		size_t nrOfHighPtElecs=0;
		size_t nrOfHighPtMus=0;
		size_t lepIndex=0; // hold the value corresponding to the index of the high PT lepton
		double M_lepton = M_mu;
		for(size_t i=0; i<elecs.size(); i++){
			if (((elecs.at(i)->PT) > electronEtLimit) && (TMath::Abs(elecs.at(i)->Eta) < 2.1)){
				nrOfHighPtLeps++;
				nrOfHighPtElecs++;
			}
		}	
		for(size_t i=0; i<muontight.size(); i++){
			if (((muontight.at(i)->PT) > muonPtLimit) && (TMath::Abs(muontight.at(i)->Eta) < 2.1)){
				nrOfHighPtLeps++;
				nrOfHighPtMus++;
			}
		}	

		// exactly one charged lepton is required
		if (nrOfHighPtLeps != 1) continue;
		if (nrOfHighPtElecs==1){
			for(size_t i=0; i<elecs.size(); i++){
				if (elecs.at(i)->PT > electronEtLimit){
					lepIndex = i;
				}
			}
		}else{
			for(size_t i=0; i<muontight.size(); i++){
				if (muontight.at(i)->PT > muonPtLimit){
					lepIndex = i;
				}
			}
		}

		TLorentzVector lepVec(0.,0.,0.,0.);
		TLorentzVector elecVec(0.,0.,0.,0.);
		TLorentzVector muVec(0.,0.,0.,0.);
		bool antiLepton = false;
		double leptonPT, leptonEta, leptonPhi;

		//met must be larger than 20 for muonic decay
		double minimumMET = 20;
		TLorentzVector missingET(0.,0.,0.,0.);
		missingET=(met.at(0)->P4());

		if (nrOfHighPtElecs==1){
			//met must be larger than 30 for electronic decay
			minimumMET = 30;
			if (met.at(0)->MET < minimumMET){
				acceptableEvent = false;
			}
			M_lepton = M_e;
			leptonPT = elecs.at(lepIndex)->PT;
			leptonEta = elecs.at(lepIndex)->Eta;
			leptonPhi = elecs.at(lepIndex)->Phi;
			lepVec.SetPtEtaPhiM(leptonPT, leptonEta, leptonPhi , M_lepton);
			if (elecs.at(lepIndex)->Charge == 1){
				antiLepton = true;
			}
		} else{
			if (met.at(0)->MET< minimumMET){
				acceptableEvent = false;
			}
			leptonPT = muontight.at(lepIndex)->PT;
			leptonEta = muontight.at(lepIndex)->Eta;
			leptonPhi = muontight.at(lepIndex)->Phi;
			lepVec.SetPtEtaPhiM(leptonPT, leptonEta, leptonPhi , M_lepton);
			if (muontight.at(lepIndex)->Charge == 1){
				antiLepton = true;
			}
		}
		if(!acceptableEvent) continue;

		//finding neutrino mass from missing transverse energy
		pzCalculator calculator;
		calculator.setLepton(lepVec);
		calculator.setMET(missingET);
		calculator.setLepMass(M_lepton);
		double pznu = calculator.getPz();
		TLorentzVector neutrinoP4(missingET.Px(),
					missingET.Py(),
					pznu,
					TMath::Sqrt(TMath::Power(missingET.Pt(),2)+TMath::Power(pznu,2)));
		neutrinoMomentumHisto->Fill(neutrinoP4.Pt());

		TLorentzVector t(0., 0., 0., 0.); //top quark
		TLorentzVector tbar(0., 0., 0., 0.); // top antiquark
		TLorentzVector bJetVec(0.,0.,0.,0.);
		TLorentzVector empty(0.,0.,0.,0.);// an Lorentz vector that stays empty
		
		// minimum Chi^2
		ttbar_solver solver;	

		double topmass = 0;
		double tbarmass = 0;
		solver.setLepton(lepVec);
		solver.setNeutrino(neutrinoP4);
		double bestchi = 100000;
		solver.setMtopLeptonic(M_T);
		for(double m = 125; m<250; m=m+5){
			solver.setMtopHadronic(m);
			for(size_t i=0; i<jet.size(); i++){
				//bjet from blv
				if(acceptBJet(jet.at(i))){
					TLorentzVector bjet1(0.,0.,0.,0.);
					TLorentzVector bjet2(0.,0.,0.,0.);
					TLorentzVector ljet1(0.,0.,0.,0.);
					TLorentzVector ljet2(0.,0.,0.,0.);
					bjet1 = makeTLorentzVector(jet.at(i));
					solver.setBJetA(bjet1);
					for(size_t j=0; j<jet.size(); j++){
						//bjet from bqq'
						if(acceptBJet(jet.at(j)) && j!=i){
							solver.setBJetB(bjet2);
							for(size_t k=0;k<jet.size();k++){
								//first light jet
								if(acceptLJet(jet.at(k))){
									ljet1=makeTLorentzVector(jet.at(k));
									solver.setLightJetA(ljet1);
									for(size_t l=k;l<jet.size();l++){
										//second light jet
										if(acceptLJet(jet.at(l)) && l!=k){
											ljet2=makeTLorentzVector(jet.at(l));
											solver.setLightJetB(ljet2);
											// vary top mass,
											// leave antitop fixed
											double chi2=solver.getChi2();
											if(chi2<bestchi){
												bestchi=chi2;	
												if(antiLepton){
													if (solver.getIsChiA()==true){
														tbar = bjet1 +ljet1 + ljet2;
														tbarmass = m;
													} else{
														tbar = bjet2 + ljet1 + ljet2;
														tbarmass = m;
													}
												} else{
													if (solver.getIsChiA()==true){
														t = bjet1 +ljet1+ljet2;
														topmass = m;
													}else{
														t = bjet2 + ljet1 + ljet2;
														topmass = m;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}//i<jet.size()
		}
		//if (t!=empty) hadronicTopMassHisto->Fill(t.M());
		//if (tbar!=empty) hadronicTbarMassHisto->Fill(tbar.M());
		if (antiLepton && (tbarmass > 0)) hadronicTbarMassHisto->Fill(tbarmass);
		if (topmass>0) hadronicTopMassHisto->Fill(topmass);

	} // for event
	processEndFunction();
}//void testanalyser::analyze


void testanalyser::postProcess(){
	/*
	 * This function can be used to analyse the output histograms, e.g. extract a signal contribution etc.
	 * The function can also be called directly on an output file with the histograms, if
	 * RunOnOutputOnly = true
	 * is set in the analyser's config file
	 *
	 * This function also represents an example of how the output of the analyser can be
	 * read-back in an external program.
	 * Just include the sampleCollection.h header and follow the procedure below
	 *
	 */

	/*
	 * Here, the input file to the extraction of parameters from the histograms is the output file
	 * of the parallelised analysis.
	 * The sampleCollection class can also be used externally for accessing the output consistently
	 */
	d_ana::sampleCollection samplecoll;
	samplecoll.readFromFile(getOutPath());

	std::vector<TString> alllegends = samplecoll.listAllLegends();
	
	if(alllegends.size()>0){
		d_ana::histoCollection histos=samplecoll.getHistos(alllegends.at(0));

		/*
		 * here, the histogram created in the analyze() function is selected and evaluated
		 * The histoCollection maintains ownership (you don't need to delete the histogram)
		 */
		const TH1* myplot=histos.getHisto("neutrinoMomentumHisto");

		std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;
		/*
		 * If the histogram is subject to changes, please clone it and take ownership
		 */

		TH1* myplot2=histos.cloneHisto("neutrinoMomentumHisto");
		
		/*
		 * do something with the histogram
		 */

		delete myplot2;
	}

	/*
	 * do the extraction here.
	 */
}

