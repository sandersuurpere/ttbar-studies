#include "TLorentzVector.h"
#include "TMath.h"

double deltaPhi(const double& a, const double& b);

template<class T>
TLorentzVector makeTLorentzVector(T* c){
	TLorentzVector out;
	out.SetPtEtaPhiM(c->PT, c->Eta, c->Phi, c->Mass);
	return out;
}

template<class T, class V>
double deltaR(T* a, V* b){
	double dphi = deltaPhi(a->Phi, b->Phi);
	double deta = a->Eta - b->Eta;
	return TMath::Sqrt(dphi*dphi + deta*deta);
}

template<class T>
bool acceptLJet(T* jet){
	bool acceptable = true;
	if ((jet->BTag) || (TMath::Abs(jet->Eta)>2.5 || (jet->PT <= 25))){
		acceptable = false;
	}
	return acceptable;
}

template<class T>
bool acceptBJet(T* jet){
	bool acceptable = true;
	if ((!(jet->BTag)) || (TMath::Abs(jet->Eta)>2.5) || (jet->PT <= 25)){
		acceptable = false;
	}
	return acceptable;
}
