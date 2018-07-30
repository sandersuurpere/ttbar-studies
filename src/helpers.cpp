#include "../interface/helpers.h"

double deltaPhi(const double& a, const double& b){
	const double pi = 3.14159265358979323846;
	double delta = (a - b);
	while (delta >= pi)  delta-= 2* pi;
	while (delta < -pi)  delta+= 2* pi;
	return delta;
}

