#include "../interface/ttbar_solver.h"

double ttbar_solver::getChi2(){
	double M_W = 80.385;
	double M_T = 172.5;
	
	//resolutions
/*	double sigmaljet = 5.0;
	double sigmabjet = 5.0;
	double sigmalepton= 1.0;
	double sigmaneutrino = 5.0;
*/
	double sigmaljet = 1.0;
	double sigmabjet = 1.0;
	double sigmalepton= 1.0;
	double sigmaneutrino = 1.0;

	double chi2 = 0;

	TLorentzVector w1 = neutrino+lepton;
	TLorentzVector w2 = ljeta+ljetb;

	double chi2_a = (TMath::Power(((ljeta+ljetb+bjeta).M()-mtop_hadr_), 2))/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet))
		+ (TMath::Power(((lepton+neutrino+bjetb).M()-mtop_lept_), 2))/(TMath::Power((sigmalepton+sigmaneutrino+sigmabjet), 2));

	double chi2_b = (TMath::Power(((ljeta+ljetb+bjetb).M()-mtop_hadr_), 2))/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet))
		+ (TMath::Power(((lepton+neutrino+bjeta).M()-mtop_lept_), 2))/(TMath::Power((sigmalepton+sigmaneutrino+sigmabjet), 2));
/*
	double chi2_1 = (M_W-w1.M())*(M_W-w1.M())/((sigmalepton+sigmaneutrino)*(sigmalepton+sigmaneutrino))
		+ (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet))
		+ (TMath::Power(((ljeta+ljetb+bjeta).M()-mtop_hadr_), 2))/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet))
		+ (TMath::Power(((lepton+neutrino+bjetb).M()-mtop_lept_), 2))/(TMath::Power((sigmalepton+sigmaneutrino+sigmabjet) ,2));

	double chi2_2 = (M_W-w1.M())*(M_W-w1.M())/((sigmalepton+sigmaneutrino)*(sigmalepton+sigmaneutrino))
		+ (M_W-w2.M())*(M_W-w2.M())/((2*sigmaljet)*(2*sigmaljet))
		+ (TMath::Power(((ljeta+ljetb+bjetb).M()-mtop_hadr_), 2))/((2*sigmaljet+sigmabjet)*(2*sigmaljet+sigmabjet))
		+ (TMath::Power(((lepton+neutrino+bjeta).M()-mtop_lept_), 2))/(TMath::Power((sigmalepton+sigmaneutrino+sigmabjet) ,2));
*/
	if (chi2_a < chi2_b){
		chi2 = chi2_a;
		setIsChiA(true);// check if using the first equation
	}else{
		chi2 = chi2_b;
		setIsChiA(false);
	}

/*
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
*/
	return chi2;
}
