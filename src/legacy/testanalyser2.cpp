#include <vector>
#include "interface/testanalyser.h"
#include "../interface/helpers.h"
#include "../interface/ttbar_solver.h"
#include "TMath.h"
using namespace std;

Double_t W_BOSON_MASS = 80.385;
Double_t TOP_MASS = 172.44;

void testanalyser::analyze(size_t childid /* this info can be used for printouts */){

	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"Muon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");

	TH1* topQuarkHisto=addPlot(new TH1D("topQuarkHisto","topHistoTitle",100,0,500),"m [GeV]","N_{e}");
	TH1* wBosonHistoL=addPlot(new TH1D("wBosonHistoL","W",100,79,130),"m [GeV]","N_{e}");
	TH1* neutrinoMomentumHisto = addPlot(new TH1D("neutrinoMomentumHisto", "nu histo",100,0,200), "pT","N");
	TH1* leadingJetHistoAll = addPlot(new TH1D("leadingJetHistoAll", "Leading jet masses (all events)",100,0,50), "m_jet","N_e");
	TH1* leadingJetHistoLeptonic = addPlot(new TH1D("leadingJetHistoLeptonic", "Leading jet masses in single lepton events",100,0,50), "m_jet","N_e");
	TH1* quarkHistoHadronicBoosted= addPlot(new TH1D("quarkHistoHadronicBoosted", "Top quark mass from full hadronic, boosted",100,0,400), "m_3jets","N_e");
	TH1* histoSL=addPlot(new TH1D("tquark_SL","t-quark m_{inv} mass in semileptonic decay",50,100,400),"M [GeV]","N");
	
	TH1* topMHistoBoost=addPlot(new TH1D("topMHistoBoost","t-quark m_{inv} mass",100,0,400),"M [GeV]","N");

	TH1* tbarMHistoBoost=addPlot(new TH1D("tbarMHistoBoost","antit-quark m_{inv} mass",100,0,400),"M [GeV]","N");

	TH1* massDiffHisto = addPlot(new TH1D("massDiffHisto", "m_t - m_{bar{t}}", 100, -50, 50), "#Delta m", "N_{events}");


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

		//At least 2 b-jets
		size_t bjets=0;
		size_t ljets=0;
		for(size_t i=0;i<jet.size();i++){
			if(fabs(jet.at(i)->Eta)<2.5){
				if(jet.at(i)->BTag){
					bjets++;
				}else{ljets++;}
			}
		}
		if(bjets<2){continue;}
		if(ljets<2){continue;}


		TLorentzVector leadingJetVectorA(0.,0.,0.,0.); //the leading jet of an event
		Double_t leadingJetPT=0;
		TLorentzVector leadingJetVectorLeptonic(0.,0.,0.,0.);
		Double_t leadingJetPTLeptonic=0;
		TLorentzVector lepVec(0.,0.,0.,0.);

		for (size_t i = 0; i<jet.size(); i++){
			if ((jet.at(i)->PT) > leadingJetPT){
				leadingJetVectorA.SetPtEtaPhiM(jet.at(i)->PT, jet.at(i)->Eta, jet.at(i)->Phi, jet.at(i)->Mass);
			}
		}

		//determine the number of leptons with pt>30
		size_t nrOfHighPtLeps=0;
		size_t nrOfHighPtElecs=0;
		size_t nrOfHighPtMus=0;
		size_t lepIndex=0; // hold the value corresponding to the index of the high PT lepton
		for(size_t i=0; i<elecs.size(); i++){
			if (((elecs.at(i)->PT) > 30) && (TMath::Abs(elecs.at(i)->Eta) < 2.1)){
				nrOfHighPtLeps=nrOfHighPtLeps+1;
				nrOfHighPtElecs=nrOfHighPtElecs+1;
			}
		}	
		for(size_t i=0; i<muontight.size(); i++){
			if (((muontight.at(i)->PT) > 30) && (TMath::Abs(muontight.at(i)->Eta) < 2.1)){
				nrOfHighPtLeps=nrOfHighPtLeps+1;
				nrOfHighPtMus=nrOfHighPtMus+1;
			}
		}	

		if (nrOfHighPtLeps != 1) continue;
		if (nrOfHighPtElecs==1){
			for(size_t i=0; i<elecs.size(); i++){
				if (elecs.at(i)->PT > 30){
					lepIndex = i;
				}
			}
		}else{
			for(size_t i=0; i<muontight.size(); i++){
				if (muontight.at(i)->PT > 30){
					lepIndex = i;
				}
			}
		}

		for (size_t i = 0; i<jet.size(); i++){
			if ((jet.at(i)->PT) > leadingJetPTLeptonic){
				leadingJetVectorLeptonic.SetPtEtaPhiM(jet.at(i)->PT, jet.at(i)->Eta, jet.at(i)->Phi, jet.at(i)->Mass);
			}
		}
		int type = 1;
		double M_W  = W_BOSON_MASS;
		double M_mu =  0.10566;
		double M_e = 0.511e-3;
		double M_lepton = M_mu;
		double emu=0;
		double pxmu=0;
		double pymu=0;
		double pzmu=0; 
		double leptonPT=0;
		double leptonEta=0;
		double leptonPhi=0;
		TLorentzVector elecVec(0.,0.,0.,0.);
		TLorentzVector muVec(0.,0.,0.,0.);
		bool antiLepton = false;

		if (nrOfHighPtElecs==1){
			M_lepton = M_e;
			leptonPT = elecs.at(lepIndex)->PT;
			leptonEta = elecs.at(lepIndex)->Eta;
			leptonPhi = elecs.at(lepIndex)->Phi;
			elecVec.SetPtEtaPhiM(leptonPT, leptonEta, leptonPhi , M_lepton);
			pymu = elecVec.Py();
			pzmu = elecVec.Pz();
			lepVec = elecVec;
			if (elecs.at(lepIndex)->Charge == 1){
				antiLepton = true;
			}
		} else{
			leptonPT = muontight.at(lepIndex)->PT;
			leptonEta = muontight.at(lepIndex)->Eta;
			leptonPhi = muontight.at(lepIndex)->Phi;
			muVec.SetPtEtaPhiM(leptonPT, leptonEta, leptonPhi , M_lepton);
			emu = muVec.E();
			pxmu = muVec.Px();
			pymu = muVec.Py();
			pzmu = muVec.Pz();
			lepVec = muVec;
			if (muontight.at(lepIndex)->Charge == 1){
				antiLepton = true;
			}
		}
		TLorentzVector missingET(0.,0.,0.,0.);
		missingET=(met.at(0)->P4());
		double pxnu = missingET.Px();
		double pynu =missingET.Py();
		double pznu = 0.;
		// use pznu = - B/2*A +/- sqrt(B*B-4*A*C)/(2*A)

		double a = M_W*M_W - M_lepton*M_lepton + 2.0*(pxmu*pxnu + pymu*pynu);
		double A = 4.0*(emu*emu - pzmu*pzmu);
		double B = -4.0*a*pzmu;
		double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;

		double tmproot = B*B - 4.0*A*C;

		if (tmproot<0) {
			bool isComplex_= true;
			pznu = - B/(2*A); // take real part of complex roots
		}
		else {
			bool isComplex_ = false;
			double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
			double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

			if (type == 0 ) {
				// two real roots, pick the one closest to pz of muon
				if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
				else pznu = tmpsol1;
				// if pznu is > 300 pick the most central root
				if ( pznu > 300. ) {
					if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
					else pznu = tmpsol2;
				}
			}//if type==0
			if (type == 1 ) {
				// two real roots, pick the one closest to pz of muon
				if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
				else pznu = tmpsol1;
			}//if type==1
			if (type == 2 ) {
				// pick the most central root.
				if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
				else pznu = tmpsol2;
			}//if type==2
			if (type == 3 ) {
				// pick the largest value of the cosine
				TVector3 p3w, p3mu;
				p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol1);
				p3mu.SetXYZ(pxmu, pymu, pzmu );
				double sinthcm1 = 2.*(p3mu.Perp(p3w))/M_W;
				p3w.SetXYZ(pxmu+pxnu, pymu+pynu, pzmu+ tmpsol2);
				double sinthcm2 = 2.*(p3mu.Perp(p3w))/M_W;

				double costhcm1 = TMath::Sqrt(1. - sinthcm1*sinthcm1);
				double costhcm2 = TMath::Sqrt(1. - sinthcm2*sinthcm2);

				if ( costhcm1 > costhcm2 ) pznu = tmpsol1;
				else pznu = tmpsol2;
			}//if type==3
		}//else
		
		TLorentzVector neutrinoP4(missingET.Px(),
							missingET.Py(),
							pznu,
							TMath::Sqrt(TMath::Power(missingET.Pt(),2)+TMath::Power(pznu,2)));

		neutrinoMomentumHisto->Fill(neutrinoP4.Pt());
		wBosonHistoL->Fill((lepVec+neutrinoP4).M());
		



		TLorentzVector t(0., 0., 0., 0.); //top quark
		TLorentzVector tbar(0., 0., 0., 0.); // top antiquark
		TLorentzVector wBosLep = lepVec + neutrinoP4; // W boson of the leptonic decay
		TLorentzVector wBosHad(0.,0.,0.,0.); // W boson of the hadronic decay
		TLorentzVector bJetVec(0.,0.,0.,0.);
		TLorentzVector empty(0.,0.,0.,0.);// an Lorentz vector that stays empty

		for (size_t i = 0; i<jet.size(); i++){
			if (acceptBJet(jet.at(i))){
				double DelR=0;
				if (nrOfHighPtElecs==1){
					DelR = deltaR(jet.at(i), elecs.at(lepIndex));
				}
				if (nrOfHighPtMus==1){
					DelR = deltaR(jet.at(i), muontight.at(lepIndex));
				}
				if (DelR<1.2){ // getting top mass from leptonic decay
					bJetVec = makeTLorentzVector(jet.at(i));
					if(antiLepton){t = bJetVec+lepVec+neutrinoP4;}
					else{tbar = bJetVec+lepVec+neutrinoP4;}
					double JetR=0;
					JetR = TMath::Sqrt(TMath::Power(jet.at(i)->DeltaEta,2)+TMath::Power(jet.at(i)->DeltaPhi,2));
					if(DelR<JetR && antiLepton){
						t = t - lepVec;
					//	topMHistoBoost->Fill(t.M());
						}
					if(DelR<JetR && !antiLepton){
						tbar = tbar - lepVec;
					//	tbarMHistoBoost->Fill(tbar.M());
					}
					//topQuarkHisto->Fill((quarkVec).M());


				}else{ // getting top mass from hadronic decay (copied from Aleksei) 
					if (t==empty && tbar==empty) continue;
					for(size_t j=0;j<jet.size();j++){
						if (acceptLJet(jet.at(j))){
							DelR = deltaR(jet.at(i),jet.at(j));
							if (DelR<1.2){
								for(size_t k=j;k<jet.size();k++){
									if(acceptLJet(jet.at(k))){
										size_t leading=0;
										if(jet.at(i)->PT > jet.at(j)->PT){
											leading=i;
										}else{leading=j;}

										if(jet.at(leading)->PT < jet.at(k)->PT){
											leading=k;
										}
										if(leading!=k){
											DelR=deltaR(jet.at(leading),jet.at(k));
										}else{
											DelR=deltaR(jet.at(i),jet.at(k));
											double DelR_2=0;
											DelR = deltaR(jet.at(j),jet.at(k));
											if(DelR<DelR_2){
												DelR=DelR_2;
											}
										}
										if(DelR<1.2){
											TLorentzVector jet1(0.,0.,0.,0.);
											jet1=makeTLorentzVector(jet.at(i));
											TLorentzVector jet2(0.,0.,0.,0.);
											jet2=makeTLorentzVector(jet.at(j));
											TLorentzVector jet3(0.,0.,0.,0.);
											jet3=makeTLorentzVector(jet.at(k));

											TLorentzVector topQ(0.,0.,0.,0.);
											if(antiLepton){
												tbar=jet1+jet2+jet3;
											}else{t=jet1+jet2+jet3;}
											//quarkHistoHadronicBoosted->Fill(topQ.M());
										}
									}
								}
							}
						}
					}
				
				}//else statement (corresponding to DelR>=1.2)
			}//if (acceptBJet) 
		}// jet loop
		if (t.M()>0 && tbar.M()>0){
		
			massDiffHisto->Fill(t.M() - tbar.M());
			topMHistoBoost->Fill(t.M());
			tbarMHistoBoost->Fill(tbar.M());
		}
		leadingJetHistoAll->Fill(leadingJetVectorA.M()); // seeing what the leading jet looks like
		if (leadingJetVectorLeptonic.M()>0)
			//the leading jet mass corresponding to leptonic decay
			leadingJetHistoLeptonic->Fill(leadingJetVectorLeptonic.M());
		
		// minimum chi^2
		ttbar_solver solver;	
		double topmass1 = 0;
		double topmass2 = 0;
		solver.setLepton(lepVec);
		solver.setNeutrino(neutrinoP4);
		double bestchi = 1000;
		for(size_t i=0; i<jet.size(); i++){
			//bjet from blv
			if(acceptBJet(jet.at(i))){
				TLorentzVector bjet1(0.,0.,0.,0.);
				TLorentzVector bjet2(0.,0.,0.,0.);
				TLorentzVector ljet1(0.,0.,0.,0.);
				TLorentzVector ljet2(0.,0.,0.,0.);
				bjet1=makeTLorentzVector(jet.at(i));
				solver.setBJetA(bjet1);
				for(size_t j=0; j<jet.size(); j++){
					//bjet from bqq'
					if(acceptBJet(jet.at(j)) && j!=i){
						bjet2=makeTLorentzVector(jet.at(j));
						solver.setBJetB(bjet2);
						for(size_t k=0;k<jet.size();k++){
							//first light jet
							if(acceptLJet(jet.at(k))){
								ljet1=makeTLorentzVector(jet.at(k));
								solver.setLightJetA(ljet1);
								for(size_t l=k;l<jet.size();l++){
									//second lightjet
									if(acceptLJet(jet.at(l)) && l!=k){
										ljet2=makeTLorentzVector(jet.at(l));
										solver.setLightJetB(ljet2);
										//comparing chi2 with previous combination
										double chi2=solver.getChi2();
										if(chi2<bestchi){
											bestchi=chi2;	
											topmass1=(bjet1+neutrinoP4+lepVec).M();
											topmass2=(bjet2+ljet1+ljet2).M();
										}
									}
								}
							}
						}
					}
				}
			}
		}//chi2 outer loop
		//histogram of semileptonic decay
		histoSL->Fill(topmass1); //from blv
		histoSL->Fill(topmass2); //from bqq'

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

	/*
	 * Example how to process the output.
	 * Usually, one would define the legendname of the histogram to be used here
	 * by hand, e.g. "signal" or "background".
	 * To make this example run in any case, I am using alllegends.at(0), which
	 * could e.g. be the signal legend.
	 *
	 * So in practise, the following would more look like
	 * samplecoll.getHistos("signal");
	 */
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

