#include "../interface/ttbar_solver.h"

double ttbar_solver::getChi2()const{
	double M_W = 80.385;
	double M_T = 172.44;
	
	//resolutions
	double sigmaljet = 5.0;
	double sigmabjet = 5.0;
	double sigmalepton= 1.0;
	double sigmaneutrino = 5.0;

	double chi2 = 0;

	if(!fh){ // if not full hadronic decay
		TLorentzVector w1 = neutrino+lepton;
		TLorentzVector top1 = w1+bjeta;
		TLorentzVector w2 = ljeta+ljetb;
		TLorentzVector top2 = w2+bjetb;

		chi2 = (M_W-w1.M())*(M_W-w1.M())/((sigmalepton+sigmaneutrino)*(sigmalepton+sigmaneutrino))
		+ (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet))
		+ (mtop_lept_-top1.M())*(mtop_lept_-top1.M())/(TMath::Power((sigmalepton+sigmaneutrino+sigmabjet) ,2))
		+ (mtop_hadr_-top2.M())*(mtop_hadr_-top2.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet));
	}else{ // if full hadronic decay
		TLorentzVector w1=ljeta+ljetb;
		TLorentzVector top1=w1+bjeta;
		TLorentzVector w2=ljetc+ljetd;
		TLorentzVector top2=w2+bjetb;
		chi2 = (M_W-w1.M())*(M_W-w1.M())/((2*sigmaljet)*(2*sigmaljet))
		+ (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet))
		+ (mtop_hadr_-top1.M())*(mtop_hadr_-top1.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet)) 
		+ (mtop_hadr_-top2.M())*(mtop_hadr_-top2.M())/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet));
	}

	return chi2;
}
