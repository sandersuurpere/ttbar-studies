// This file contains the old code that we currently don't use anymore


#include <vector>
#include "interface/testanalyser.h"
using namespace std;

Double_t W_BOSON_MASS = 80.385;
Double_t TOP_MASS = 172.44;

void testanalyser::analyze(size_t childid /* this info can be used for printouts */){

	/*
	 * This skeleton analyser runs directly on the Delphes output.
	 * It can be used to create histograms directly or a skim.
	 * If a skim is created, a new input configuration will be written automatically
	 * and stored in the output directory together with the ntuples.
	 * The skim can contain delphes objects again or can be flat. This is up
	 * to the user.
	 * Examples for both are given here.
	 *
	 * The same skeleton can be used to read the skim. Please refer to the comments
	 * marked with "==SKIM=="
	 *
	 * These parts are commented, since the code is supposed to work as an example without
	 * modifications on Delphes output directly.
	 */



	/*
	 * Define the branches that are to be considered for the analysis
	 * This branch handler (notice the "d")
	 * is used to run directly in Delphes output.
	 * For skimmed ntuples, see below
	 */
	d_ana::dBranchHandler<Electron> elecs(tree(),"Electron");
	/*
	 * Other branches might be the following
	 * (for a full list, please inspect the Delphes sample root file with root)
	 * For the Delphes class description, see $DELPHES_PATH/classes/DelphesClasses.h
	 */
	//
	d_ana::dBranchHandler<HepMCEvent>  event(tree(),"Event");
	d_ana::dBranchHandler<Jet>         genjet(tree(),"GenJet");
	d_ana::dBranchHandler<Jet>         jet(tree(),"Jet");
	d_ana::dBranchHandler<Muon>        muontight(tree(),"Muon");
	d_ana::dBranchHandler<MissingET>   met(tree(),"MissingET");


	/* ==SKIM==
	 *
	 * If a skim of the Delphes outout was created in a way indicated
	 * further below, use the tBranchHandler (please notive the "t")
	 * to access vectors of objects...
	 *
	 */
	// d_ana::tBranchHandler<std::vector<Electron> > electrons(tree(),"Electrons");

	/*==SKIM==
	 *
	 * Or an object directly
	 *
	 */
	//d_ana::tBranchHandler<MissingET> met(tree(),"MET");



	/*
	 * Always use this function to add a new histogram (can also be 2D)!
	 * Histograms created this way are automatically added to the output file
	 */
	TH1* histo=addPlot(new TH1D("elecHistoName1","elecHistoTitle1",100,0,100),"p_{T} [GeV]","N_{e}");
	TH1* etaHisto=addPlot(new TH1D("elecEtaHisto","elecEtaHistoTitle",100,0,3),"#eta","N_{e}");

	TH1* jetHisto=addPlot(new TH1D("jetHistoName1","JetHistoTitle1",100,20,100),"p_{T} [GeV]","N_{e}");
	TH1* jetEtaHisto=addPlot(new TH1D("jetEtaHisto","jetEtaHistoTitle",100,0,5),"#eta","N_{e}");
	
		

	TH1* jetMassHisto=addPlot(new TH1D("jetMassHisto","Jet mass",250,0,1000),"m_{jet}","N_{e}");
	
	TH1* jetSumHisto2=addPlot(new TH1D("jetSumHisto2","Jet sum mass",450,0,1000),"m_{2jets}","N_{e}");

	TH1* jetSumHisto3=addPlot(new TH1D("jetSumHisto3","Jet sum mass",450,0,1000),"m_{3jets}","N_{e}");
	
	TH1* jetMassOfEventHisto =addPlot(new TH1D("jetMassOfEventHisto","#sum_{event} m_{jet}",450,0,1000),"m_{jet, total}","N_{e}");

	TH1* wBosonHisto=addPlot(new TH1D("wBosonHisto","W boson",100,79,81),"m_{2jets}","N_{e}");

	TH1* wBosonHisto2=addPlot(new TH1D("wBosonHisto2","W boson",100,0,300),"m_{2jets}","N_{e}");

	TH1* topQuarkHisto=addPlot(new TH1D("topQuarkHisto","top",100,0,500),"m_{3jets}","N_{e}");
	/*
	 * If (optionally) a skim or a flat ntuple is to be created, please use the following function to initialize
	 * the tree.
	 * The output files will be written automatically, and a config file will be created.
	 */
	//TTree* myskim=addTree();
	/*
	 * Add a simple branch to the skim
	 */
	//Double_t elecPt=0;
	//myskim->Branch("elecPt", &elecPt);
	/*
	 * Or store a vector of objects (also possible to store only one object)
	 */
	//std::vector<Electron> skimmedelecs;
	//myskim->Branch("Electrons",&skimmedelecs);



	size_t nevents=tree()->entries();
	if(isTestMode())
		nevents/=100;
	for(size_t eventno=0;eventno<nevents;eventno++){
		/*
		 * The following two lines report the status and set the event link
		 * Do not remove!
		 */
		reportStatus(eventno,nevents);
		tree()->setEntry(eventno);

		TLorentzVector v1;

		/*
		 * Begin the event-by-event analysis
		 */
		for(size_t i=0;i<elecs.size();i++){
			histo->Fill(elecs.at(i)->PT);
			etaHisto->Fill(elecs.at(i)->Eta);	
		}

		for(size_t i=0; i<jet.size(); i++){
			jetHisto->Fill(jet.at(i)->PT);
			jetEtaHisto->Fill(jet.at(i)->Eta);
		}

		TLorentzVector totalLorentzVectorOfEvent;
		for(size_t i=0; i<jet.size(); i++){
			v1.SetPtEtaPhiM(jet.at(i)->PT,
					jet.at(i)->Eta,
					jet.at(i)->Phi,
					jet.at(i)->Mass);
			jetMassHisto->Fill(v1.M());
			totalLorentzVectorOfEvent +=v1;
			for(size_t j=i; j<jet.size(); j++){
				if (i != j){
					TLorentzVector v2;
					v2.SetPtEtaPhiM(jet.at(j)->PT,
							jet.at(j)->Eta,
							jet.at(j)->Phi,
							jet.at(j)->Mass);
					jetSumHisto2->Fill((v2+v1).M());
				
					for(size_t k=j; k<jet.size(); k++){
					
						if (i != j && i!=j && j!=k){
							TLorentzVector v3;
							v3.SetPtEtaPhiM(jet.at(k)->PT,
									jet.at(k)->Eta,
									jet.at(k)->Phi,
									jet.at(k)->Mass);
							jetSumHisto3->Fill((v1+v2+v3).M());
						}//if statement
					}// for k
				}//if statement
			}//for j
		}//for i
		
		
		//Method 1 (the inceorrect method) for combining two jets to get W 
		TLorentzVector jetVector1;
		Double_t suggestedMassOfW; // 2jet mass combination that is closest to the W boson mass
		for(size_t i=0; i<jet.size(); i++){
			jetVector1.SetPtEtaPhiM(jet.at(i)->PT,
						jet.at(i)->Eta,
						jet.at(i)->Phi,
						jet.at(i)->Mass);
			for(size_t j=i; j<jet.size(); j++){
				if (i != j){
					TLorentzVector jetVector2;
					jetVector2.SetPtEtaPhiM(jet.at(j)->PT,
							jet.at(j)->Eta,
							jet.at(j)->Phi,
							jet.at(j)->Mass);
					Double_t candidateMassOfW = (jetVector1+jetVector2).M();
					if (TMath::Abs((candidateMassOfW-W_BOSON_MASS))<TMath::Abs((suggestedMassOfW-W_BOSON_MASS))){
						suggestedMassOfW = candidateMassOfW;
					}
				}//if (i!=j)
			}//for j
		}//for i
	
		wBosonHisto->Fill(suggestedMassOfW);


		//method 2 for combining two jets to get W (basically Aleksei's method, aka the correct method :) )
		TLorentzVector jetVec01;
		int suggestedIndex=0; // this index holds element number corresponding to the best fit
		// a vector that keeps track of which i and j indeces are already used:
		vector<bool> usedOrNot(jet.size()); // true = used index, false = unused index
		for(size_t i=0; i<jet.size(); i++){
			if((usedOrNot[i]==false) && !(jet.at(i)->BTag)){
				usedOrNot.at(i)=true;
				jetVec01.SetPtEtaPhiM(	jet.at(i)->PT,
							jet.at(i)->Eta,
							jet.at(i)->Phi,
							jet.at(i)->Mass);
				TLorentzVector suggestedW; // 2jet combination that is closest to the W boson
				for(size_t j=i; j<jet.size(); j++){
					if ((i != j) && (usedOrNot.at(j)==false)){
						TLorentzVector jetVec02;
						jetVec02.SetPtEtaPhiM(	jet.at(j)->PT,
									jet.at(j)->Eta,
									jet.at(j)->Phi,
									jet.at(j)->Mass);
						TLorentzVector candidateOfW = jetVec01+jetVec02;
						if (TMath::Abs((candidateOfW.M()-W_BOSON_MASS))<TMath::Abs((suggestedW.M()-W_BOSON_MASS))){
							suggestedW = candidateOfW;
							suggestedIndex = j;
						}
					}//if (i!=j)
				}//for j
				if (suggestedW.M()>5){
					wBosonHisto2->Fill(suggestedW.M());
					usedOrNot.at(suggestedIndex)=true;
				}
			}/*if(usedOrNot[i]==false)*/ 
		}//for i
	

		// Combining 3 jets to get the histogram for top quark mass
		TLorentzVector jetV1;
		int suggestedIndexJ=0; // this index holds element number corresponding to the best fit
		int suggestedIndexK=0; // this index holds element number corresponding to the best fit
		// a vector that keeps track of which i,j,k indeces are already used:
		vector<bool> used(jet.size()); // true = used index, false = unused index
		for(size_t i=0; i<jet.size(); i++){
			if((!used.at(i)) && (jet.at(i)->BTag)){
				used.at(i)=true;
				jetV1.SetPtEtaPhiM(	jet.at(i)->PT,
							jet.at(i)->Eta,
							jet.at(i)->Phi,
							jet.at(i)->Mass);
				TLorentzVector suggestedTop;// 3jet combination that is closest to the top quark
				TLorentzVector wSuggested;// 2jet combo that is closest to W 
				for(size_t j=i; j<jet.size(); j++){
					if ((i != j) && (!used.at(j)) && !(jet.at(j)->BTag)){
						TLorentzVector jetV2;
						jetV2.SetPtEtaPhiM(	jet.at(j)->PT,
									jet.at(j)->Eta,
									jet.at(j)->Phi,
									jet.at(j)->Mass);
						for(size_t k=j; k<jet.size(); k++){
							if (i != j && i!=j && j!=k && !used.at(k) && !(jet.at(j)->BTag)){
								TLorentzVector jetV3;
								jetV3.SetPtEtaPhiM(	jet.at(k)->PT,
											jet.at(k)->Eta,
											jet.at(k)->Phi,
											jet.at(k)->Mass);
								TLorentzVector wCandidate = jetV2+jetV3;
								TLorentzVector topCandidate = jetV1 + jetV2+jetV3;
								if ((TMath::Abs(wCandidate.M()-W_BOSON_MASS)
									<TMath::Abs(wSuggested.M()-W_BOSON_MASS))
									&& (TMath::Abs(topCandidate.M()-TOP_MASS)
									<TMath::Abs(suggestedTop.M()-TOP_MASS))){
									suggestedTop = topCandidate;
									suggestedIndexK = k;
									suggestedIndexJ = j;
								}
							}//if !used.at(k) 
						}// for k
					}//if (i!=j)
				}//for j
				if(suggestedTop.M()){
					topQuarkHisto->Fill(suggestedTop.M());
				}
				used.at(suggestedIndexK)=true;
				used.at(suggestedIndexJ)=true;
			}/*if(usedOrNot[i]==false)*/ 
		}//for i
	
      if(bJets.size()>=2 && lightJets.size()>=2)
        {
          //visible system
          TLorentzVector visSystem(leptons[0].p4()+bJets[0].p4()+bJets[1].p4()+lightJets[0].p4()+lightJets[1].p4());
          //determine the neutrino kinematics
          TLorentzVector met(0, 0, 0, 0);
          met.SetPtEtaPhiM(ev.met_pt[0], 0, ev.met_phi[0], 0.);
          neutrinoPzComputer.SetMET(met);
          neutrinoPzComputer.SetLepton(leptons[0].p4());
          float nupz=neutrinoPzComputer.Calculate();
          TLorentzVector neutrinoP4(met.Px(),met.Py(),nupz ,TMath::Sqrt(TMath::Power(met.Pt(),2)+TMath::Power(nupz,2)));
          //ttbar system
          TLorentzVector ttbarSystem(visSystem+neutrinoP4);
          ht.fill("ht",         scalarht, plotwgts);          
          ht.fill("mttbar_cen", ttbarSystem.M(),   plotwgts);
        }

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
		const TH1* myplot=histos.getHisto("elecHistoName1");

		const TH1* etaPlot=histos.getHisto("elecEtaHisto");

		
		const TH1* myJetPlot=histos.getHisto("jetHistoName1");

		const TH1* jetEtaPlot=histos.getHisto("jetEtaHisto");

		const TH1* jetMassPlot = histos.getHisto("jetMassHisto");

		const TH1* jetSumPlot2 = histos.getHisto("jetSumHisto2");

		const TH1* jetSumPlot3 = histos.getHisto("jetSumHisto3");

		std::cout << "(example output): the integral is " << myplot->Integral() <<std::endl;

		/*
		 * If the histogram is subject to changes, please clone it and take ownership
		 */

		TH1* myplot2=histos.cloneHisto("elecHistoName1");
		
		
		/*
		 * do something with the histogram
		 */

		delete myplot2;
	}

	/*
	 * do the extraction here.
	 */



}

