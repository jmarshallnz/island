#ifndef _CLUSTER_H_
#define _CLUSTER_H_

#include "myutils.h"
#include "mydouble.h"
#include <fstream>
#include <map>

class Cluster {
	int ng;							// # groups
  myutils::Vector< myutils::Matrix<int> > MLST;		// haps for each
  myutils::Vector<int> size;				// size of each
  myutils::Vector<mydouble> SIZE;			// size of each
  myutils::Vector<int> nST;				// # unique STs in each
  myutils::Matrix<int> nalleles;			// nalleles[i][j] # unique alleles in group i at locus j
  myutils::Vector< myutils::Vector<double> > FREQ;		// freq of STs in each group
  myutils::Vector< myutils::Vector<double> > ABUN;		// abundance of STs in each group
  myutils::Vector< myutils::Vector<mydouble> > freq;	// freq of STs in each group
  myutils::Vector< myutils::Vector<mydouble> > abun;	// abundance of STs in each group
  myutils::Matrix< myutils::Vector<double> > acount;	// acount[i][j][k] gives the count, in pop i, and locus j, of allele k
  mydouble one;

  int nloc;						// # loci
	bool init;

	myutils::Matrix<int> human;				// those sampled from humans

	mydouble punique;
	myutils::Vector<mydouble> puniq,psame,pdiff;
	myutils::Matrix<bool> human_unique;
	myutils::Vector< myutils::Matrix<bool> > beast_unique;
	bool ****same;
	bool *****ksame;

public:
	Cluster() {
		init = false;
		nloc = 7;
		one = 1.0;
		same = NULL;
		ksame = NULL;
	}
	void open_all(const char* filename);

	// mcmc6f infers M and R from seqs of known origin, and runs 100 side-chains to infer F given M and R
	void mcmc6f(const double alpha, const double beta, const double gamma, myutils::Random &ran, const int niter, const int thin, const char* filename);
	void mcmc6f(const double alpha, const double beta, const double gamma_, myutils::Random &ran, const int niter, const int thin, myutils::Matrix<double> &traces, std::map<int, myutils::Matrix<double> > phi_out);

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
		/* free memory */
		for(int i = 0; i < human.nrows(); i++) {
		  for(int ii = 0; ii < ng; ii++) {
		    for(int jj = 0; jj < nST[ii]; jj++) {
		      delete[] same[i][ii][jj];
		    }
		    delete[] same[i][ii];
		  }
		  delete[] same[i];
		}
		delete[] same;

		for(int i = 0; i < ng; i++) {
		  for(int j = 0; j < nST[i]; j++) {
		    for(int ii = 0; ii < ng; ii++) {
		      for(int jj = 0; jj < nST[ii]; jj++) {
		        delete[] ksame[i][j][ii][jj];
		      }
		      delete[] ksame[i][j][ii];
		    }
		    delete[] ksame[i][j];
		  }
		  delete[] ksame[i];
		}
		delete[] ksame;
	}

	int multinom(myutils::Vector<double> &p, myutils::Random &ran);
	mydouble likHi6(const int id, const int i, myutils::Matrix<double> &a, myutils::Matrix< myutils::Vector<double> > &b, myutils::Matrix<double> &r);
	mydouble known_source_lik6_composite(myutils::Matrix<double> &a, myutils::Matrix< myutils::Vector<double> > &b, myutils::Matrix<double> &r);

	void recalc_b(myutils::Matrix<double> &a, myutils::Matrix< myutils::Vector<double> > &b);
	void calc_A(myutils::Matrix<mydouble> &a, myutils::Matrix<double> &A);
	void calc_Ai(myutils::Matrix<mydouble> &a, myutils::Matrix<double> &A, const int i);
	void calc_R(myutils::Matrix<mydouble> &r, myutils::Matrix<double> &R);
	void calc_Ri(myutils::Matrix<mydouble> &r, myutils::Matrix<double> &R, const int i);
	void precalc();
	template<typename T> void pSWAP(myutils::Vector<T> &a, myutils::Vector<T> &b);
	template<typename T> void pSWAP(myutils::Matrix<T> &a, myutils::Matrix<T> &b);
};

#endif//_CLUSTER_H_
