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
/*
double deltaR(double aPhi, double bPhi, double aEta, double bEta){
	double dphi = deltaPhi(aPhi, bPhi);
	double deta = aEta - bEta;
	return TMath::Sqrt(dphi*dphi + deta*deta);
}
*/
