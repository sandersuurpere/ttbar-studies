#include "../interface/pzCalculator.h"

double M_W = 80.385;

double pzCalculator::getPz()const{

	double emu=lepton.E();
	double pxmu=lepton.Px();
	double pymu=lepton.Py();
	double pzmu=lepton.Pz();
	double pxnu=met.Px();
	double pynu=met.Py();
	double pznu=0.;

	double a = M_W*M_W - M_lepton*M_lepton + 2.0*(pxmu*pxnu + pymu*pynu);
	double A = 4.0*(emu*emu - pzmu*pzmu);
	double B = -4.0*a*pzmu; 
	double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;
	double tmproot = B*B - 4.0*A*C;

	if (tmproot<0) {
		pznu = - B/(2*A); // take real part of complex roots
	}else{
		double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
		double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);

		// two real roots, pick the one closest to pz of muon
		if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) {
			pznu = tmpsol2;
		}else{pznu = tmpsol1;}
	
	} 
	return pznu;

}
