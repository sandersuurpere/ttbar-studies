#include "TLorentzVector.h"

class ttbar_solver{
	public:

		ttbar_solver():mtop_hadr_(172.5), mtop_lept_(172.5){}

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

		void setMtopHadronic(const double& m){
			mtop_hadr_ = m;
		}
		void setMtopLeptonic(const double& m){
			mtop_lept_ = m;
		}

	private:
		TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);
		TLorentzVector ljeta, ljetb, ljetc, ljetd, bjeta, bjetb, lepton, neutrino;
		bool fh = false;
		double mtop_hadr_, mtop_lept_;
};
