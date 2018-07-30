/*
 * testanalyser.h
 *
 *  Created on: 24 Aug 2016
 *      Author: jkiesele
 */

#ifndef testanalyser_H_
#define testanalyser_H_

#include "interface/basicAnalyzer.h"
#include "interface/sampleCollection.h"
#include "classes/DelphesClasses.h"


class testanalyser: public d_ana::basicAnalyzer{
public:
	testanalyser():d_ana::basicAnalyzer(){}
	~testanalyser(){}


private:
	void analyze(size_t id);

	void postProcess();
};





#endif /* testanalyser_H_ */
