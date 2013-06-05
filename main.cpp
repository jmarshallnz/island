/*	Simulation study for "source". 

*/
#include "myutils.h"
#include "cluster.h"

using namespace myutils;

int main(const int argc, const char* argv[]) {
	if(argc!=6 && argc!=7) error("SYNTAX: in_file out_file niter thinning alpha [seed]");
	const char* train_file = argv[1];
	const char* out_file = argv[2];
	const int niter = /*100000;*/atoi(argv[3]);
	const int thinning = /*50;*/atoi(argv[4]);
	const double alpha = atof(argv[5]);
	Random ran;
	if(argc==7) ran.setseed(atoi(argv[6]));
	cout << "Seed set to " << ran.getseed() << endl;

	Cluster clust;
	clust.open_all(train_file);
	double HS = 0., HT = 0.;
	cout << "FST = " << clust.FST(HS,HT);
	cout << "\tHS = " << HS << "\tHT = " << HT << endl;

#if defined(_MODEL6)
	cout << "Linked model" << endl;
	clust.mcmc6f(alpha,1.0,1.0,ran,niter,thinning,out_file);
#else
	error("This program was compiled under an unrecognized model");
#endif
	return 0;
}

