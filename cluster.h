#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include "myutils.h"
#include "tsv.h"
#include "mydouble.h"
//#include "powell.h"
#include <fstream>
#include <time.h>

using namespace myutils;

class Cluster;

//class ML : public PowellFunction {
//public:
//  Cluster &clust;
//  /* Constructor */
//  ML(Cluster &clust_in) : clust(clust_in) {}
//  /* Call f(theta,rho,rhostar) */
//  double f(const vector<double> &x);
//};

class Cluster {
public:
	int ng;							// # groups
	Vector< Matrix<int> > MLST;		// haps for each
	Vector<int> size;				// size of each
	Vector<mydouble> SIZE;			// size of each
	Vector<int> nST;				// # unique STs in each
	Matrix<int> nalleles;			// nalleles[i][j] # unique alleles in group i at locus j
	Vector< Vector<double> > FREQ;		// freq of STs in each group
	Vector< Vector<double> > ABUN;		// abundance of STs in each group
	Vector< Vector<mydouble> > freq;	// freq of STs in each group
	Vector< Vector<mydouble> > abun;	// abundance of STs in each group
	Matrix< Vector<double> > acount;	// acount[i][j][k] gives the count, in pop i, and locus j, of allele k
	int nloc;						// # loci
	bool init;
	double A;

	Matrix<int> human;				// those sampled from humans
	Vector<int> htime;	///< mapping of human time
	int         ntime;	///< total times

	Vector< Matrix<mydouble> > ksl2;

	mydouble punique;
	Vector<mydouble> puniq,psame,pdiff;
	Matrix<bool> human_unique;
	Vector< Matrix<bool> > beast_unique;
	bool ****same;
	mydouble one;
	bool *****ksame;
	Vector<int> hid;	///< previous identical human for saving likelihood computation
	Vector<int> htid;	///< previous identical human and time for saving likelihood computation
	Vector<int> G;
	Matrix<mydouble> simLIK;
	Vector<int> identicals;
	Matrix<int> simMLST;

public:
	Cluster() {
		init = false;
		nloc = 7;
		A = 1.0e-6;
		one = 1.0;
	}
	void open(const char* filename);
	void open_human(const char* filename);
	void open_all(const char* filesources, const char *filehuman);
	void mcmc(const double th0, const double rh0, const double rhs0, const double alpha, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc(const double th0, const double rh0, const double rhs0, const double alpha, Random &ran, const int niter, const int thin, ofstream& out);
	void mcmc2(const double th0, const double rh0, const double rhs0, const double alpha, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc3(const double th0, const double rh0, const double mi0, const double et0, const double alpha, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc3(const double th0, const double rh0, const double mi0, const double et0, const double alpha, Random &ran, const int niter, const int thin, ofstream& out);
	// mcmc4 uses the unlinked model
	void mcmc4(const double alpha, const double beta, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc4(const double alpha, const double beta, Random &ran, const int niter, const int thin, ofstream &out);
	void mcmc4b(const double alpha, const double beta, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc4b(const double alpha, const double beta, Random &ran, const int niter, const int thin, ofstream &out);
	void mcmc4c(const double alpha, const double beta, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc4c(const double alpha, const double beta, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	// mcmc5 ???
	void mcmc5(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc5(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out);
	void mcmc5b(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc5b(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out);
	void mcmc5c(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc5c(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out);
	// mcmc6 use the linkage model
	// mcmc6b uses sequences of known and unknown origin to jointly infer F,M,R but terminally fails to mix
	void mcmc6b(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6b(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	// mcmc6c uses only sequences of known origin to infer M and R. its method of inferring F though is incorrect - the F chain spends much time in 'burn-in'
	void mcmc6c(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6c(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	// mcmc6d uses sequences of known and unknown origin to jointly infer F,M,R using a different likelihood ratio calc_lik6d(), but this is incorrect and fails to converge
	void mcmc6d(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6d(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	// mcmc6e infers M and R from seqs of known origin, and draws G assuming a constant, uniform F
	void mcmc6e(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6e(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	// mcmc6f infers M and R from seqs of known origin, and runs 100 side-chains to infer F given M and R
	void mcmc6f(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6f(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, ofstream &out, ofstream &o2, ofstream &o3);
	// mcmc6g is the same as mcmc6f, but also performs posterior predictive goodness-of-fit simulations
	void mcmc6g(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6g(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, ofstream &out, ofstream &o2, ofstream &o3, ofstream &o4);
	// mcmc6h is the same as mcmc6f but F is fixed and uniform
	void mcmc6h(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6h(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, ofstream &out, ofstream &o2, ofstream &o3);
	// mcmc7 is uses the linkage model with group-and-locus-specific mutation and recombination rates
	void mcmc7c(const double alpha, const double beta, const double gamma, const double delta, const double epsilon, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc7c(const double alpha, const double beta, const double gamma, const double delta, const double epsilon, Random &ran, const int niter, const int thin, ofstream &out, const char* gfilename);
	~Cluster() {
		if(init) {
			int i,j;
			for(i=0;i<ng;i++) {
				for(j=0;j<nloc;j++) {
					acount[i][j].resize(0);
				}
				MLST[i].resize(0,0);
			}
		}
	}

	int multinom(Vector<double> &p, Random &ran);
	mydouble likHi(const int id, const int g, const double theta, const double rho, const double rhostar);
	mydouble likHi(const int id, const int g, const double theta, const double rho, const double M, const double eta);
	mydouble likHi3(const int id, const int g, const double theta, const double rho, const double M, const double eta);
	mydouble likHi4(const int id, const int g, Matrix<double> &a, Matrix< Vector<double> > &b);
	mydouble likHi5(const int id, const int g, Matrix<double> &a, Matrix< Vector<double> > &b, Vector<double> &r);
	mydouble likHi5(const int id, const int g, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
//	mydouble likHi5(const int id, const int i, Matrix<mydouble> &a, Matrix< Vector<double> > &b, Matrix<mydouble> &r);
	mydouble likHi6(const int id, const int i, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
	mydouble likHi7(const int id, const int i, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
	void gamma(Vector<int> &g, Vector<double> &gam, const double alpha) {
		int i;
		for(i=0;i<gam.size();i++) gam[i] = alpha;
		for(i=0;i<g.size();i++) ++gam[g[i]];
	}
	mydouble lik(const double theta, const double rho, const double rhostar);
	mydouble likprint(const double theta, const double rho, const double rhostar);
	mydouble lik(const int id, const double theta, const double rho, const double rhostar);
	double gammln(const double xx);
	mydouble full_lik(Vector<double> &gam, const double alpha, Matrix< Vector<int> > &ass, Vector<double> &f);
	mydouble known_source_lik(const double theta, const double rho, const double rhostar, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik(const double theta, const double rho, const double M, const double eta, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik2(const double theta, const double rho, const double M, const double eta, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik3(const double theta, const double rho, const double M, const double eta, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik_composite(const double theta, const double rho, const double rhostar, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik3_composite(const double theta, const double rho, const double M, const double eta, Vector<double> &SIZE, Matrix< Vector<double> > &ACOUNT);
	mydouble known_source_lik4_composite(Matrix<double> &a, Matrix< Vector<double> > &b);
	mydouble known_source_lik5_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Vector<double> &r);
	mydouble known_source_lik5_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
//	mydouble known_source_lik5_composite(Matrix<mydouble> &a, Matrix< Vector<double> > &b, Matrix<mydouble> &r);
	mydouble known_source_lik6_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
	mydouble known_source_lik7_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Vector< Matrix<double> > &P);
	double FST();
	double FST(double &HS, double &HT);
	void recalc_b(Matrix<double> &a, Matrix< Vector<double> > &b);
//	void recalc_b(Matrix<mydouble> &a, Matrix< Vector<double> > &b);
	double sigma(const double x);
	void calc_sum(mydouble *v, const int n);
	void calc_A(Matrix<mydouble> &a, Matrix<double> &A);
	void calc_Ai(Matrix<mydouble> &a, Matrix<double> &A, const int i);
	void calc_F(Vector<mydouble> &f, Vector<mydouble> &F);
	void calc_R(Matrix<mydouble> &r, Matrix<double> &R);
	void calc_Ri(Matrix<mydouble> &r, Matrix<double> &R, const int i);
	void precalc();
	mydouble calc_lik6(Matrix<mydouble> &LIKHI, Matrix<double> &A, Matrix<mydouble> &a, Matrix< Vector<double> > &b, Matrix<double> &R, Vector<mydouble> &F);
	mydouble calc_lik6(Matrix<mydouble> &LIKHI_use, Matrix<mydouble> &LIKHI_notuse, Vector<mydouble> &F_prime);
	mydouble calc_lik6d(Matrix<mydouble> &LIKHI, Matrix<double> &A, Matrix<mydouble> &a, Matrix< Vector<double> > &b, Matrix<double> &R, Vector<mydouble> &F);
	mydouble calc_lik7(Matrix<mydouble> &LIKHI, Matrix<double> &A, Matrix<mydouble> &a, Matrix< Vector<double> > &b, Vector< Matrix<double> > &P, Vector<mydouble> &F);
	mydouble calc_lik7(Matrix<mydouble> &LIKHI_use, Matrix<mydouble> &LIKHI_notuse, Vector<mydouble> &F_prime);
//	void calc_rho(Vector<mydouble> &rho, Vector<mydouble> &RHO);
//	void calc_theta(Vector<mydouble> &theta, Vector<mydouble> &THETA);
	void recalc_P(Matrix<mydouble> &r, Vector<mydouble> &theta, Vector<mydouble> &rho, Vector< Matrix<double> > &P);
	void recalc_Pi(Matrix<mydouble> &r, Vector<mydouble> &theta, Vector<mydouble> &rho, Vector< Matrix<double> > &P, const int popid);
	template<typename T> void pSWAP(Vector<T> &a, Vector<T> &b);
	template<typename T> void pSWAP(Matrix<T> &a, Matrix<T> &b);
	void simHi6(Vector<mydouble> &F, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r, Matrix<mydouble> &LIKHI, ofstream &o4, Random &ran);
	void posteriorP6g(const char* param_file, Random &ran, const char* out_file);
};

#endif//_CLUSTER_H_
