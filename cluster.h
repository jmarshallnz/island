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

	mydouble punique;
	Vector<mydouble> puniq,psame,pdiff;
	Matrix<bool> human_unique;
	Vector< Matrix<bool> > beast_unique;
	bool ****same;
	mydouble one;
	bool *****ksame;

public:
	Cluster() {
		init = false;
		nloc = 7;
		A = 1.0e-6;
		one = 1.0;
	}
	void open(const char* filename);
	void open_human(const char* filename);
	void open_all(const char* filename);

	// mcmc6f infers M and R from seqs of known origin, and runs 100 side-chains to infer F given M and R
	void mcmc6f(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6f(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, ofstream &out, ofstream &o3, const std::string &filename);

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
	mydouble likHi6(const int id, const int i, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);
	mydouble known_source_lik6_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r);

	void recalc_b(Matrix<double> &a, Matrix< Vector<double> > &b);
	void calc_A(Matrix<mydouble> &a, Matrix<double> &A);
	void calc_Ai(Matrix<mydouble> &a, Matrix<double> &A, const int i);
	void calc_R(Matrix<mydouble> &r, Matrix<double> &R);
	void calc_Ri(Matrix<mydouble> &r, Matrix<double> &R, const int i);
	void precalc();
	template<typename T> void pSWAP(Vector<T> &a, Vector<T> &b);
	template<typename T> void pSWAP(Matrix<T> &a, Matrix<T> &b);
};

#endif//_CLUSTER_H_
