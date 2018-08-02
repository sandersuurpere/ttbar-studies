#include "TLorentzVector.h"

class ttbar_solver{
	public:
		void setLightJetA(const TLorentzVector& jet){
			ljeta=jet;
		}
		void setLightJetB(TLorentzVector jet){
			ljetb=jet;
		}
		void setLightJetC(TLorentzVector jet){
			ljetc=jet;
			fh = true;
		}
		void setLightJetD(TLorentzVector jet){
			ljetc=jet;
			fh = true;
		}
		void setBJetA(TLorentzVector jet){
			bjeta=jet;
		}
		void setBJetB(TLorentzVector jet){
			bjetb=jet;
		}
		void setLepton(TLorentzVector lep){
			lepton=lep;
			fh = false;
		}
		void setNeutrino(TLorentzVector neu){
			neutrino=neu;
		}
		double getChi2()const;

	private:
		TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);
		TLorentzVector ljeta, ljetb, ljetc, ljetd, bjeta, bjetb, lepton, neutrino;
		bool fh = false;
};
