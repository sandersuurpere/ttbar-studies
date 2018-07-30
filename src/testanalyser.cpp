#include <vector>
#include "interface/testanalyser.h"
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
 		TLorentzVector leadingJetVectorA; //the leading jet of an event
		Double_t leadingJetPT;
		TLorentzVector leadingJetVectorLeptonic;
		Double_t leadingJetPTLeptonic;
		TLorentzVector lepVec;
		for (size_t i = 0; i<jet.size(); i++){
			if ((jet.at(i)->PT) > leadingJetPT){
				leadingJetVectorA.SetPtEtaPhiM(jet.at(i)->PT, jet.at(i)->Eta, jet.at(i)->Phi, jet.at(i)->Mass);
			}
		}
		if((elecs.size()==1 && ((elecs.at(0)->PT) > 30) && (TMath::Abs(elecs.at(0)->Eta) < 2.1) && muontight.size()==0) ||
 			(muontight.size()==1 && ((muontight.at(0)->PT) > 30) && (TMath::Abs(muontight.at(0)->Eta) < 2.1) && elecs.size()==0)){ 
			//count n b-jets
			int nrOfBJets;	
			for(size_t ij=0; ij<jet.size(); ij++) 
				{
				if(jet.at(ij)->BTag) nrOfBJets+=1;
				}
			//if (nrOfBJets<2) continue;
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

			double emu; 
			double pxmu;
			double pymu;
			double pzmu;

			if (elecs.size()==1) M_lepton = M_e;
			TLorentzVector elecVec;
			if (elecs.size()==1){
				elecVec.SetPtEtaPhiM(elecs.at(0)->PT, elecs.at(0)->Eta, elecs.at(0)->Phi, M_e);
				emu = elecVec.E();
				pxmu = elecVec.Px();
				pymu = elecVec.Py();
				pzmu = elecVec.Pz();
				lepVec = elecVec;
			}
			TLorentzVector muVec;
			if (elecs.size()==0){
				muVec.SetPtEtaPhiM(muontight.at(0)->PT, muontight.at(0)->Eta, muontight.at(0)->Phi, M_mu);
				emu = muVec.E();
				pxmu = muVec.Px();
				pymu = muVec.Py();
				pzmu = muVec.Pz();
				lepVec = muVec;
			}
			TLorentzVector missingET;
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
			TLorentzVector bJetVec;
			for (size_t i = 0; i<jet.size(); i++){
				if (jet.at(i)->BTag && (TMath::Abs(jet.at(i)->Eta) < 2.5)){
					//TLorentzVector temporaryVec;
					//temporaryVec.SetPtEtaPhiM(jet.at(i)->PT, jet.at(i)->Eta, jet.at(i)->Phi, jet.at(i)->Mass);
					//if (TMath::Abs(temporaryVec.M()-TOP_MASS)<(TMath::Abs(bJetVec.M()-TOP_MASS)))
					//	bJetVec = temporaryVec;
					bJetVec.SetPtEtaPhiM(jet.at(i)->PT, jet.at(i)->Eta, jet.at(i)->Phi, jet.at(i)->Mass);

					topQuarkHisto->Fill((bJetVec+lepVec+neutrinoP4).M());
				}
			}
		}
		leadingJetHistoAll->Fill(leadingJetVectorA.M());
		if (leadingJetVectorLeptonic.M()>0)
			leadingJetHistoLeptonic->Fill(leadingJetVectorLeptonic.M());
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

