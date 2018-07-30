#include "TLorentzVector.h"
#include "Math.h"

double deltaPhi(const double& a, const double& b);

template<class T>
TLorentzVector makeTLorentzVector(T* c){
	TLorentzVector out;
	out.SetPtEtaPhiM(c->PT, c->Eta, c->Phi, c->Mass);
	return out;
}

template<class T>
double deltaR(T* a, T* b){
	double dphi = deltaPhi(a->Phi, b->Phi);
	double deta = a->Eta - b->Eta;
	return TMath::Sqrt(phi*dphi + deta*deta);
}
