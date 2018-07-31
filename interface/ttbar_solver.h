#include "TLorentzVector.h"

class ttbar_solver{
	public:
		void setLightJetA(const TLorentzVector& jet){
		TLorentzVector ljeta=jet;
		}
		void setLightJetB(TLorentzVector jet){
		TLorentzVector ljetb=jet;
		}
		void setBJetA(TLorentzVector jet){
		TLorentzVector bjeta=jet;
		}
		void setBJetB(TLorentzVector jet){
		TLorentzVector bjetb=jet;
		}
		void setLepton(TLorentzVector lep){
		TLorentzVector lepton=lep;
		}
		void setNeutrino(TLorentzVector neu){
		TLorentzVector neutrino=neu;
		}
		double getChi2()const;

	private:
		TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);
		TLorentzVector ljeta, ljetb, bjeta, bjetb, leptona, neutrino;

};
