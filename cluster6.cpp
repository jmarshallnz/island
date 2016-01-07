#include "cluster.h"

#include <sstream> // for the stream stuff
#include <iostream>

using namespace myutils;

mydouble Cluster::known_source_lik6_composite(Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r) {
#if defined(_MODEL4)
	return known_source_lik4_composite(a,b);
#elif defined(_FLAT_LIKELIHOOD)
	return mydouble(1.0);
#endif
	int i,j,ii,jj,l;
	mydouble lik = 1.0;
	/* Cycle through each unique ST in each group, taking account of abundance of the STs */
	for(i=0;i<ng;i++) {
		punique = mydouble(a[i][ng]);
		for(j=0;j<nST[i];j++) {
			mydouble ncopiesj = ABUN[i][j];
			mydouble l_j(0.0);
			for(l=0;l<nloc;l++) {
				int allele = MLST[i][j][l];
				double ac = acount[i][l][allele];
				double ac_ = (ac*(double)size[i]-1.0)/(double)(size[i]-1);
				double bk = b[i][l][allele] - a[i][i]*ac + a[i][i]*ac_;
				double b_ = r[i][0] * bk + r[i][1] * (1.0-a[i][ng]);// if no rec then must be same as CF (barring mutation)
				if(fabs(b_)<1.0e-7) {
					b_ = 0.0;
				}
				psame[l] = mydouble(b_);
				b_ = r[i][0] * bk;	// different so must have been recombination
				if(fabs(b_)<1.0e-7) {
					b_ = 0.0;
				}
				pdiff[l] = mydouble(b_);
			}
			for(ii=0;ii<ng;ii++) {						//	Cycle through source of the clonal frame
				mydouble l_ii(0.0);
				mydouble mii(a[i][ii]/(1.0-a[i][ng]));
				for(jj=0;jj<nST[ii];jj++) {				//	Cycle through each ST from that source
					mydouble ncopiesjj = (i==ii && j==jj) ? abun[ii][jj]-MIN(abun[ii][jj],one)
						: abun[ii][jj];
					mydouble l_jj = mii;
					bool *BEAST_UNIQUE = beast_unique[i][j];
					bool *SAME = ksame[i][j][ii][jj];
					mydouble *PDIFF = pdiff.element;
					mydouble *PSAME = psame.element;
					for(l=0;l<nloc;l++,BEAST_UNIQUE++,SAME++,PDIFF++,PSAME++) {
						if(*BEAST_UNIQUE) {				// new allele (allow some rounding error)
							l_jj *= punique;
						}
						else if(*SAME) {				// previously observed and same as CF
							l_jj *= *PSAME;
						}
						else {							// previously observed but different to CF
							l_jj *= *PDIFF;
						}
					}
					l_ii += l_jj * ncopiesjj;
				}
				l_j += l_ii / SIZE[ii];
			}
			lik *= l_j^ncopiesj;
		}
	}
	return lik;
}

mydouble Cluster::likHi6(const int id, const int i, Matrix<double> &a, Matrix< Vector<double> > &b, Matrix<double> &r) {
#if defined(_MODEL4)
	return likHi4(id,i,a,b);
#elif defined(_FLAT_LIKELIHOOD)
	return mydouble(1.0);
#endif
	int ii,jj,l;
/// NOTE: Little a in this function is A everywhere else!!!

//	punique = mydouble(a[i][ng]);						// MAKE SURE THIS IS SET BEFORE CALLING likHi6()
	for(l=0;l<nloc;l++) {
		int human_allele = human[id][l];
		pdiff[l] = mydouble(MAX(r[i][0] * b[i][l][human_allele],0.0));
		psame[l] = mydouble(MAX(r[i][0] * b[i][l][human_allele] + r[i][1] * (1.0-a[i][ng]),0.0));
	}
	mydouble lik(0.0);
	for(ii=0;ii<ng;ii++) {								// Cycle through source of the clonal frame
		mydouble mii(a[i][ii]/(1.0-a[i][ng]));
		mydouble l_ii(0.0);
		for(jj=0;jj<nST[ii];jj++) {
			mydouble l_jj = mii;						//	Cycle through each ST from that source
			bool* HUMAN_UNIQUE = human_unique[id];
			bool* SAME = same[id][ii][jj];
			mydouble* PSAME = psame.element;
			mydouble* PDIFF = pdiff.element;
			for(l=0;l<nloc;l++,HUMAN_UNIQUE++,SAME++,PSAME++,PDIFF++) {
				if(*HUMAN_UNIQUE) {						// new allele (allow some rounding error)
					l_jj *= punique;
				}
				else if(*SAME) {						// previously observed and same as CF
					l_jj *= *PSAME;
				}
				else {									// previously observed but different to CF
					l_jj *= *PDIFF;
				}
			}
			mydouble &ncopiesjj = abun[ii][jj];
			l_ii += l_jj * ncopiesjj;
		}
		lik += l_ii / SIZE[ii];
	}
	return lik;
}

void Cluster::precalc() {

	human_unique = Matrix<bool>(human.nrows(),nloc);
	beast_unique = Vector< Matrix<bool> >(ng);
	for (int i = 0; i < ng; i++)
	  beast_unique[i] = Matrix<bool>(nST[i],nloc);

	same = new bool***[human.nrows()];
	for(int i = 0; i < human.nrows(); i++) {
		same[i] = new bool**[ng];
		for(int ii = 0; ii < ng; ii++) {
			same[i][ii] = new bool*[nST[ii]];
			for(int jj = 0; jj < nST[ii]; jj++) {
				same[i][ii][jj] = new bool[nloc];
				for(int l = 0; l < nloc; l++) {
					same[i][ii][jj][l] = (human[i][l]==MLST[ii][jj][l]);
				}
			}
		}
		for(int l = 0; l < nloc; l++) {
			int human_allele = human[i][l];
		  human_unique[i][l] = (human_allele>=acount[ng][l].size()
                            || acount[ng][l][human_allele]==0);
		}
	}
	puniq = Vector<mydouble>(nloc);
	psame = Vector<mydouble>(nloc);
	pdiff = Vector<mydouble>(nloc);

	ksame = new bool****[ng];
	for(int i = 0; i < ng; i++) {
		ksame[i] = new bool***[nST[i]];
		for(int j = 0; j < nST[i]; j++) {
			ksame[i][j] = new bool**[ng];
			for(int ii = 0; ii < ng; ii++) {
				ksame[i][j][ii] = new bool*[nST[ii]];
				for(int jj = 0; jj < nST[ii]; jj++) {
					ksame[i][j][ii][jj] = new bool[nloc];
					for(int l = 0; l < nloc; l++) {
						ksame[i][j][ii][jj][l] = (MLST[i][j][l] == MLST[ii][jj][l]);
					}
				}
			}
			for(int l = 0; l < nloc; l++) {
				int allele = MLST[i][j][l];
				double num = acount[ng][l][allele] * (double)size[ng];
				beast_unique[i][j][l] = (num < 1.1);
			}
		}
	}
}

/* This version uses the clonal frame version of the likelihood */
void Cluster::mcmc6f(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename) {

	/* Create output matrices */
	Matrix<double> out_mat(niter/thin+1, 1+ng*(ng+1)+ng+4); // iter, A, r, loglik, loglik2, logalpha, move
  std::map<int, Matrix<double> > phi_out;

  /* Run model */
	mcmc6f(alpha,beta,gamma,ran,niter,thin,out_mat,phi_out,filename);

	/* Dump output to text file(s) */

	/* Open the file */
	std::ofstream out(filename);
	out << "iter";
	const char tab = '\t';
	for (int i = 0; i < ng; i++) for(int j = 0; j < ng+1; j++) out << tab << "A[" << i << "," << j << "]";
	for (int i = 0; i < ng; i++) out << tab << "r" << i;
	out << tab << "loglik";
	out << tab << "loglik2";
	out << tab << "logalpha";
	out << tab << "move";
	out << std::endl;

	for (int i = 0; i < out_mat.nrows(); i++) {
	  out << out_mat[i][0];
	  for (int j = 1; j < out_mat.ncols(); j++) {
	    out << tab << out_mat[i][j];
	  }
	  out << std::endl;
	}
	out.close();

	/* Dump it to the file */
	for (std::map< int, Matrix<double> >::iterator i = phi_out.begin(); i != phi_out.end(); ++i) {
  	std::stringstream s; s << i->first;
  	std::string file = "phi_" + s.str() + "_" + filename;
  	std::ofstream o(file.c_str());
  	char tab = '\t';
  	o << "Isolate";
  	for (int i = 0; i < human.ncols(); i++)
  	  o << tab << "Loci" << i;
  	for (int i = 0; i < ng; i++)
  	  o << tab << "Source" << i;
  	o << "\n";
  	for (int h = 0; h < human.nrows(); h++) {
  	  o << h;
  	  for (int i = 0; i < human.ncols(); i++)
  	    o << tab << human[h][i];
  	  for (int j = 0; j < ng; j++)
  	    o << tab << i->second[h][j];
  	  o << "\n";
  	}
  	o.close();
	}
}

void append_traces(int iter, Matrix<double> &A, Matrix<double> &R, double lik1, double lik2, double logalpha, int move, Matrix<double> &traces, int trace_row) {
  int col = 0;
  traces[trace_row][col++] = iter;
  for (int i = 0; i < A.nrows(); i++) {
    for (int j = 0; j < A.ncols(); j++)
      traces[trace_row][col++] = A[i][j];
  }
  for (int i = 0; i < R.nrows(); i++)
    traces[trace_row][col++] = R[i][0];
  traces[trace_row][col++] = lik1;
  traces[trace_row][col++] = lik2;
  traces[trace_row][col++] = logalpha;
  traces[trace_row][col++] = move;
}

/* This version uses the clonal frame version of the likelihood */
void Cluster::mcmc6f(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, Matrix<double> &traces, std::map<int, Matrix<double> > phi_out, const std::string &filename) {
	precalc();
	int i,j;
	/* Initialize the Markov chain */
	int use = 0; int notuse = (int)!use;
	Vector<double> BETA(ng+1,beta);			///<	Dirichlet hyperparameters of migration matrix (beta==1)
	Vector< Matrix<mydouble> > a(2);
	a[use] = Matrix<mydouble>(ng,ng+2);		///<	Reparametrised migration matrix. a[,ng+1]=sum(a[,1:ng]), a[,ng]=mutation rate
	a[notuse] = Matrix<mydouble>(ng,ng+2);
	Vector< Matrix<double> > A(2);
	A[use] = Matrix<double>(ng,ng+1);		///<	Migration matrix M, M[ng] = mutation rates?
	A[notuse] = Matrix<double>(ng,ng+1);
	const bool a_constraint = false;
	for(i=0;i<ng;i++) {
		while(true) {
			for(j=0;j<ng+1;j++) {
				a[use][i][j] = ran.gamma(1.,BETA[j]);
			}

			if(!a_constraint) break;
			mydouble amax = a[use][i][0];
			for(j=1;j<ng;j++) if(a[use][i][j]>amax) amax = a[use][i][j];
			if(a[use][i][i]==amax) break;
		}
	}
	calc_A(a[use],A[use]);	///< Does transformation to M

	Vector< Matrix< Vector<double> > > b(2);	///< b[use][grp][loc][i] = sum_j(freq[grp][loc][i] * M[grp][j])
	b[use] = Matrix< Vector<double> >(ng,nloc);
	b[notuse] = Matrix< Vector<double> >(ng,nloc);
	for(i=0;i<ng;i++) {
		for(j=0;j<nloc;j++) {
			b[use][i][j] = Vector<double>(acount[i][j].size());	///< number of alleles on group i locus j
			b[notuse][i][j] = Vector<double>(acount[i][j].size());	///< number of alleles on group i locus j
		}
	}
	recalc_b(A[use],b[use]);
	Vector< Matrix<mydouble> > r(2);	///< Reparameterised per-group recombination rates
	r[use] = Matrix<mydouble>(ng,3);	///< r[u,grp,3] = sum(r[u,grp,1:2])
	r[notuse] = Matrix<mydouble>(ng,3);
	Vector< Matrix<double> > R(2);
	R[use] = Matrix<double>(ng,2);		///< R[u,grp,1:2] = r[u,grp,1:2]/r[u,grp,3]
	R[notuse] = Matrix<double>(ng,2);
	Vector<double> GAMMA_(ng,gamma_);
	for(i=0;i<ng;i++) {
		//r[i] = ran.beta(gamma_,gamma_);
		for(j=0;j<2;j++) {
			r[use][i][j] = ran.gamma(1.,gamma_);
		}
	}
	calc_R(r[use],R[use]);

	/* Storage for likelihoods */
	mydouble likelihood = known_source_lik6_composite(A[use],b[use],R[use]);

	/* Proposal probabilities */
	Vector<double> proprob(6,0.0);
	proprob[0] = 42./5.;							//	Update A: switching proposal
	proprob[1] = 42.;							//	Update A: log-normal proposal
	proprob[4] = 12./5.;							//	Update r: switching proposal
	proprob[5] = 12.;							//	Update r: log-normal proposal
	double tot = 0.0;
	for(i=0;i<proprob.size();i++) tot += proprob[i];
	for(i=0;i<proprob.size();i++) proprob[i] /= tot;

	double sigma_a = 0.5;							//	factor for normal proposal in MH change of a (case 1)
	double sigma_f = 0.5;							//	factor for normal proposal in MH change of f (case 3)
	double sigma_r = 0.5;							//	factor for normal proposal in MH change of r (case 5)

	/* Trace output matrix */
	int trace_row = 0;
	append_traces(0, A[use], R[use], likelihood.LOG(), likelihood.LOG(), 0, NAN, traces, trace_row++);

	clock_t start = clock(), current;
	clock_t next = start + (clock_t)CLOCKS_PER_SEC;
	std::cout << "Done 0 of " << niter << " iterations";

	mydouble newlik, logalpha;
	int iter, fiter, move, ctr = 0;
	const int fniter = 20000;
	const int fburnin = 10000;
	const int inc = MAX((int)floor((double)niter*.9/100.),1);
	const int burnin = (int)floor((double)niter*.1);
	for(iter=0;iter<niter;iter++) {
		if(iter>=burnin && (iter-burnin)%inc==0) {

			/* Now dump the likelihoods for this iteration.
			   Weird that this has nothing to do with the thinning */

			/* Compute likelihood of human isolate from each source */
			Matrix<mydouble> phi(human.nrows(), ng);
			{
			  int j;
			  for(int h = 0; h < human.nrows(); h++) {
          // calculate the likelihood
          for (int j = 0; j < ng; j++) {
            punique = a[use][j][ng];              // NOTE USE of little a here!!!
            phi[h][j] = likHi6(h,j,A[use],b[use],R[use]);
          }
			  }
			}

			/* Dump it to a matrix */
			Matrix<double> phi_out_i(human.nrows(), ng);
			for (int h = 0; h < human.nrows(); h++) {
			  for (int j = 0; j < ng; j++)
			    phi_out_i[h][j] = phi[h][j].LOG();
			}
			phi_out[iter] = phi_out_i;
		}
		else {
			newlik = likelihood;
			logalpha = 1;
			move = multinom(proprob,ran);			//	random sweep for proposing moves
			switch(move) {
				case 0:	{// update A: switching proposal
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "mig" matrix
					int id1 = ran.discrete(0,ng);		// Principal elements of mig matrix to change
					int id2 = ran.discrete(0,ng-1);
					if(id2==id1) id2 = ng;
					if(a_constraint && (id1==popid || id2==popid)) break;
					a[notuse] = a[use];
					A[notuse] = A[use];
					SWAP(a[notuse][popid][id1],a[notuse][popid][id2]);
					calc_Ai(a[notuse],A[notuse],popid);
					logalpha = 1.0;
					// Prior ratio equals 1 because prior is symmetric
					// Hastings ratio equals 1 because proposal is symmetric
					// Likelihood ratio
					recalc_b(A[notuse],b[notuse]);
					mydouble oldlik = likelihood;
					newlik = known_source_lik6_composite(A[notuse],b[notuse],R[use]);

					logalpha *= newlik / oldlik;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						r[notuse] = r[use];
						R[notuse] = R[use];
						SWAP(use,notuse);
						likelihood = newlik;
					}
					else { // reject
					}
					break;
				}
				case 1:	{// update A: log-normal proposal
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "mig" matrix
					int id = ran.discrete(0,ng);		// Principal element of mig matrix to change
					a[notuse] = a[use];
					A[notuse] = A[use];
					mydouble *ap = a[use][popid], *ap_prime = a[notuse][popid];
					ap_prime[id].setlog(ran.normal(ap[id].LOG(),sigma_a));
					bool reject = false;
					if(a_constraint) {
						mydouble ap_primemax = ap_prime[0];
						for(j=1;j<ng;j++) if(ap_prime[j]>ap_primemax) ap_primemax = ap_prime[j];
						if(ap_prime[popid]!=ap_primemax) reject = true;
					}
					if(reject) break;
					calc_Ai(a[notuse],A[notuse],popid);
					// Prior-Hastings ratio
					logalpha.setlog(ap[id].todouble()-ap_prime[id].todouble());
					logalpha *= (ap_prime[id]/ap[id])^(beta);
					// Likelihood ratio
					recalc_b(A[notuse],b[notuse]);
					mydouble oldlik = likelihood;
					newlik = known_source_lik6_composite(A[notuse],b[notuse],R[use]);

					logalpha *= newlik / oldlik;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						r[notuse] = r[use];
						R[notuse] = R[use];
						SWAP(use,notuse);
						likelihood = newlik;
					}
					else { // reject
					}
					break;
				}
				case 4: {// update r (switching move)
					int popid = ran.discrete(0,ng-1);
					r[notuse] = r[use];
					R[notuse] = R[use];
					SWAP(r[notuse][popid][0],r[notuse][popid][1]);
					calc_Ri(r[notuse],R[notuse],popid);
					logalpha = 1.0;
					// Prior ratio equals 1 because prior is symmetric
					// Symmetric proposal so Hastings ratio equals 1
					// Likelihood ratio
					mydouble oldlik = likelihood;
					newlik = known_source_lik6_composite(A[use],b[use],R[notuse]);

					logalpha *= newlik / oldlik;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						a[notuse] = a[use];
						A[notuse] = A[use];
						b[notuse] = b[use];
						SWAP(use,notuse);
						likelihood = newlik;
					}
					else { // reject
					}
					break;
				}
				case 5:	{// update r (log-normal move)
					int popid = ran.discrete(0,ng-1);	// Choose the source for which to change the "rec" parameter
					int id = ran.discrete(0,1);			// Change one or other of the gamma components
					r[notuse] = r[use];
					R[notuse] = R[use];
					mydouble *rp = r[use][popid], *rp_prime = r[notuse][popid];
					rp_prime[id].setlog(ran.normal(rp[id].LOG(),sigma_r));
					calc_Ri(r[notuse],R[notuse],popid);
					// Prior-Hastings ratio
					logalpha.setlog(rp[id].todouble()-rp_prime[id].todouble());
					logalpha *= (rp_prime[id]/rp[id])^(gamma_);
					// Likelihood ratio
					mydouble oldlik = likelihood;
					newlik = known_source_lik6_composite(A[use],b[use],R[notuse]);

					logalpha *= newlik / oldlik;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						a[notuse] = a[use];
						A[notuse] = A[use];
						b[notuse] = b[use];
						SWAP(use,notuse);
						likelihood = newlik;
					}
					else { // reject
					}
					break;
				}
				default: {
					error("Move not recognised");
				}
			}
		}

		/* output thinned traces of island model fit */
		if((iter+1)%thin==0)
		  append_traces(iter+1, A[use], R[use], likelihood.LOG(), newlik.LOG(), logalpha.LOG(), move, traces, trace_row++);

		if((current=clock())>next) {
		  std::cout << "\rDone " << (iter+1) << " of " << niter << " iterations in " << (double)(current-start)/CLOCKS_PER_SEC << " s " << std::flush;
			next = current + (clock_t)CLOCKS_PER_SEC;
		}
	}
	std::cout << std::endl;
}

// helper stuff below here

// Assumes A is correctly sized
void Cluster::calc_A(Matrix<mydouble> &a, Matrix<double> &A) {
  int i,j;
  const int n = a.ncols()-1;
  for(i=0;i<a.nrows();i++) {
    a[i][n] = 0.0;
    for(j=0;j<n;j++) a[i][n] += a[i][j];
    for(j=0;j<n;j++) A[i][j] = (a[i][j]/a[i][n]).todouble();
  }
}

// Assumes A is correctly sized
void Cluster::calc_Ai(Matrix<mydouble> &a, Matrix<double> &A, const int i) {
  int j;
  const int n = a.ncols()-1;
  a[i][n] = 0.0;
  for(j=0;j<n;j++) a[i][n] += a[i][j];
  for(j=0;j<n;j++) A[i][j] = (a[i][j]/a[i][n]).todouble();
}

// Assumes R is correctly sized
void Cluster::calc_R(Matrix<mydouble> &r, Matrix<double> &R) {
  int i,j;
  const int n = r.ncols()-1;
  for(i=0;i<r.nrows();i++) {
    r[i][n] = 0.0;
    for(j=0;j<n;j++) r[i][n] += r[i][j];
    for(j=0;j<n;j++) R[i][j] = (r[i][j]/r[i][n]).todouble();
  }
}

// Assumes R is correctly sized
void Cluster::calc_Ri(Matrix<mydouble> &r, Matrix<double> &R, const int i) {
  int j;
  const int n = r.ncols()-1;
  r[i][n] = 0.0;
  for(j=0;j<n;j++) r[i][n] += r[i][j];
  for(j=0;j<n;j++) R[i][j] = (r[i][j]/r[i][n]).todouble();
}

void Cluster::recalc_b(Matrix<double> &a, Matrix< Vector<double> > &b) {
  int i,j,k,l;
  for(i=0;i<ng;i++) {
    for(j=0;j<nloc;j++) {
      for(k=0;k<acount[i][j].size();k++) {
        b[i][j][k] = 0.0;
        for(l=0;l<ng;l++) {
          b[i][j][k] += acount[l][j][k] * a[i][l];
          //					bk[i][j][k] += (acount[l][j][k]*(double)size[l]-1.0)/(double)(size[l]-1) * a[i][l];
        }
      }
    }
  }
}
