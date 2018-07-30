class ttbar_solver{
	public:
	    void setLightJetA(const TLorentzVector& jet){ljeta=jet;}
	    void setLightJetB(TLorentzVector jet){ljetb=jet;}

	    void setBJetA(TLorentzVector jet){bjeta=jet;}
	    void setBJetB(TLorentzVector jet){bjetb=jet;}

	    void setLepton(TLorentzVector lep){lepton=lep;}
	    void setNeutrino(TLorentzVector neu){neutrino=neu;}

	    double getChi2()const;

	private:
	    TLorentzVector getNeutrinoVector(TLorentzVector lepton, TLorentzVector met);
	    TLorentzVector ljeta, ljetb, bjeta, bjetb, leptona, neutrino;

}

