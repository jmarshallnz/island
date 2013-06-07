#include "cluster.h"

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
	int i,j,ii,jj,l;
	const int h = human.nrows();

	human_unique = Matrix<bool>(h,nloc);
	beast_unique = Vector< Matrix<bool> >(ng);
	for(i=0;i<ng;i++) beast_unique[i] = Matrix<bool>(nST[i],nloc);

	same = new bool***[h];
	for(i=0;i<h;i++) {
		same[i] = new bool**[ng];
		for(ii=0;ii<ng;ii++) {
			same[i][ii] = new bool*[nST[ii]];
			for(jj=0;jj<nST[ii];jj++) {
				same[i][ii][jj] = new bool[nloc];
				for(l=0;l<nloc;l++) {
					same[i][ii][jj][l] = (human[i][l]==MLST[ii][jj][l]);
				}
			}
		}
		for(l=0;l<nloc;l++) {
			int human_allele = human[i][l];
			if(human_allele>=acount[ng][l].size()
				|| acount[ng][l][human_allele]==0) human_unique[i][l] = true;
			else human_unique[i][l] = false;
		}
	}
	puniq = Vector<mydouble>(nloc);
	psame = Vector<mydouble>(nloc);
	pdiff = Vector<mydouble>(nloc);

	ksame = new bool****[ng];
	for(i=0;i<ng;i++) {
		ksame[i] = new bool***[nST[i]];
		for(j=0;j<nST[i];j++) {
			ksame[i][j] = new bool**[ng];
			for(ii=0;ii<ng;ii++) {
				ksame[i][j][ii] = new bool*[nST[ii]];
				for(jj=0;jj<nST[ii];jj++) {
					ksame[i][j][ii][jj] = new bool[nloc];
					for(l=0;l<nloc;l++) {
						ksame[i][j][ii][jj][l] = (MLST[i][j][l]==MLST[ii][jj][l]);
					}
				}
			}
			for(l=0;l<nloc;l++) {
				int allele = MLST[i][j][l];
				double num = acount[ng][l][allele] * (double)size[ng];
				if(num<1.1) beast_unique[i][j][l] = true;
				else beast_unique[i][j][l] = false;
			}
		}
	}

	// Identifies haplotypes that are identical, to save recalculating the likelihoods
	hid = Vector<int>(h);
	for(i=0;i<h;i++) {
		for(j=0;j<i;j++) {
			bool same = true;
			for(l=0;l<nloc;l++) {
				if(human[i][l]!=human[j][l]) {
					same = false;
					break;
				}
			}
			if(same) break;
		}
		hid[i] = j;
	}
	G = Vector<int>(h);
	simLIK = Matrix<mydouble>(h,ng);
	identicals = Vector<int>(h);
	simMLST = Matrix<int>(h,nloc);
}

mydouble Cluster::calc_lik6(Matrix<mydouble> &LIKHI, Matrix<double> &A, Matrix<mydouble> &a, Matrix< Vector<double> > &b, Matrix<double> &R, Vector<mydouble> &F) {
#if defined(_FLAT_LIKELIHOOD)
	return mydouble(1.0);
#endif
	int i,j;
	const int h = LIKHI.nrows();
	mydouble lik = 1.0;
	for(i=0;i<h;i++) {
		if(hid[i]<i) { ///< have an identical haplotype, so copy the likelihood over
			const int ii = hid[i];
			for(j=0;j<ng;j++) LIKHI[i][j] = LIKHI[ii][j];
			LIKHI[i][ng] = LIKHI[ii][ng];
		}
		else {	// calculate the likelihood
			LIKHI[i][ng] = 0.0;
			for(j=0;j<ng;j++) {
				punique = a[j][ng];			// NOTE USE of little a here!!!
				LIKHI[i][j] = likHi6(i,j,A,b,R);
				LIKHI[i][ng] += F[j] * LIKHI[i][j];
			}
		}
		lik *= LIKHI[i][ng];
	}
	return lik;
}

/*! \brief used to determine the new likelihood after a swap in F
  i.e. the individual within-group likelihoods will be the same, so just multiply
       up by the group likelihood (F)
 */
mydouble Cluster::calc_lik6(Matrix<mydouble> &LIKHI_use, Matrix<mydouble> &LIKHI_notuse, Vector<mydouble> &F_prime) {
#if defined(_FLAT_LIKELIHOOD)
	return mydouble(1.0);
#endif
	int i,j;
	const int h = LIKHI_use.nrows();
	mydouble lik = 1.0;
	for(i=0;i<h;i++) {
		if(hid[i]<i) {
			const int ii = hid[i];
			LIKHI_notuse[i][ng] = LIKHI_notuse[ii][ng];
		}
		else {
			LIKHI_notuse[i][ng] = 0.0;
			for(j=0;j<ng;j++) {
				LIKHI_notuse[i][ng] += F_prime[j] * LIKHI_use[i][j];
			}
		}
		lik *= LIKHI_notuse[i][ng] / LIKHI_use[i][ng];
	}
	return lik;
}

void calc_logit_F(const Vector<double> f, Vector<mydouble> &F)
{
  	/* Assumes F is one larger than f */
	double fs = 0;
	for (int i = 0; i < f.size(); i++)
	{
		fs += exp(f[i]);
		F[i].setlog(f[i]);		///< equivalent to F[i] = exp(f[i])
  	}
	F[f.size()].setlog(0); 			///< equivalent to F[f.size()] = 1
	fs += 1;
	for (int i = 0; i < F.size(); i++)
 		F[i] /= fs;
}

/* This version uses the clonal frame version of the likelihood */
void Cluster::mcmc6f(const double alpha, const double beta, const double gamma, Random &ran, const int niter, const int thin, const char* filename) {
	int i, j, h = human.nrows();
	/* Open the file */
	ofstream out(filename);
	string gfilename = string("g_") + string(filename);
	ofstream o2(gfilename.c_str());
	string ffilename = string("f_") + string(filename);
	ofstream o3(ffilename.c_str());
	char tab = '\t';
	out << "iter";
	for(i=0;i<ng;i++) for(j=0;j<ng+1;j++) out << tab << "A[" << i << "," << j << "]";
	for(i=0;i<ng;i++) out << tab << "r" << i;
	out << tab << "loglik";
	out << tab << "loglik2";
	out << tab << "logalpha";
	out << tab << "move";
	out << endl;
	o3 << "iter";
	for(i=0;i<ng;i++) o3 << tab << "f" << i;
	o3 << endl;
	return mcmc6f(alpha,beta,gamma,ran,niter,thin,out,o2,o3);
}

/* This version uses the clonal frame version of the likelihood */
void Cluster::mcmc6f(const double alpha, const double beta, const double gamma_, Random &ran, const int niter, const int thin, ofstream &out, ofstream &o2, ofstream &o3) {
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

	int h = human.nrows();

	/* Storage for likelihoods */
	Vector< Matrix<mydouble> > LIKHI(2);
	LIKHI[use] = Matrix<mydouble>(h,ng+1);		///< likelihood of each human isolate
	LIKHI[notuse] = Matrix<mydouble>(h,ng+1);	///< likelihood of each human isolate
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

	/* Source probabilities */
	Vector<double> ALPHA(ng-1, 0);			///< Mean of assignment proportion (logit scale)
	Vector<double> TAU(ng-1, 1);			///< Std dev of assignment proportion (logit scale)
	Vector<double> f(ng-1,0), f_prime(ng-1);	///< F values on the logit scale
	Vector<mydouble> F(ng), F_prime(ng);		///< Probability of source

	/* Output to file */
	char tab = '\t';
	out << 0;
	for(i=0;i<ng;i++) {
		for(j=0;j<ng+1;j++) out << tab << A[use][i][j];
	}
	for(i=0;i<ng;i++) out << tab << R[use][i][0];
	out << tab << likelihood.LOG();
	out << tab << likelihood.LOG();
	out << tab << "0";
	out << tab << "NA";
	out << endl;

	clock_t start = clock(), current;
	clock_t next = start + (clock_t)CLOCKS_PER_SEC;
	cout << "Done 0 of " << niter << " iterations";

	Matrix<mydouble> GLIK(h,ng,mydouble(0.0));

	mydouble newlik, logalpha;
	int iter, fiter, move, ctr = 0;
	const int fniter = 20000;
	const int fburnin = 10000;
	const int inc = MAX((int)floor((double)niter*.9/100.),1);
	const int burnin = (int)floor((double)niter*.1);
	for(iter=0;iter<niter;iter++) {
		if(iter>=burnin && (iter-burnin)%inc==0) {
			/* Draw our hyper-prior parameters */
			for(j=0; j < f.size(); j++)
				f[j] = ran.normal(ALPHA[j], TAU[j]);
			/* Calculate our F's */
			calc_logit_F(f,F);

			/* Compute likelihood */
			mydouble flik = calc_lik6(LIKHI[use],A[use],a[use],b[use],R[use],F);

			/* Side chain */
			for(fiter=0;fiter<fniter;fiter++) {
				move = (ran.U()<.05) ? 2 : 3; /* Switching likelihood */
				switch(move) {
				case 2: {
					/* update f_logit (switching move)
					   choose two values to switch, ensuring they're different */
					int id1 = ran.discrete(0,f.size()-1);
					int id2 = ran.discrete(0,f.size()-2);
					if(id2==id1) id2 = f.size()-1;
					f_prime = f;
					SWAP(f_prime[id1],f_prime[id2]);
					calc_logit_F(f_prime,F_prime);
					logalpha = 1.0;
					// Prior ratio equals 1 because prior is symmetric
					// Symmetric proposal so Hastings ratio equals 1
					// Likelihood ratio
					mydouble lik_ratio = calc_lik6(LIKHI[use],LIKHI[notuse],F_prime);

					logalpha *= lik_ratio;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						f = f_prime;
						F = F_prime;
						for(i=0;i<h;i++) LIKHI[use][i][ng] = LIKHI[notuse][i][ng];
						flik *= lik_ratio;
					}
					else { // reject
					}
					break;
				}
				case 3: {// update f (II. Using a Metropolis-Hastings step with normal proposal)
					int id = ran.discrete(0,f.size()-1);
					f_prime = f;
					f_prime[id] = ran.normal(f[id],sigma_f);
					calc_logit_F(f_prime,F_prime);
					// Prior-Hastings ratio = Proposal(f,f')/Proposal(f',f) * Prior(f')/Prior(f)
					// Proposal is symmetric, so this drops down to the prior. Prior is
					// exp((f^2-f'^2)/(2*tau^2))
					logalpha.setlog((f[id]*f[id]-f_prime[id]*f_prime[id])/(2*TAU[id]*TAU[id]));
					// Likelihood ratio
					mydouble lik_ratio = calc_lik6(LIKHI[use],LIKHI[notuse],F_prime);

					logalpha *= lik_ratio;
					if(logalpha.LOG()>=0.0 || ran.U()<logalpha.todouble()) {	// accept
						f = f_prime;
						F = F_prime;
						for(i=0;i<h;i++) LIKHI[use][i][ng] = LIKHI[notuse][i][ng];
						flik *= lik_ratio;
					}
					else { // reject
					}
					break;
				}
				default:
					error("Cluster::mcmc6f(): Undefined move in f chain");
				}
				if(fiter%100==0) {
					o3 << fiter;
					for(i=0;i<ng;i++) o3 << tab << F[i].todouble();
					o3 << endl;
					if(fiter>=fburnin) {
						++ctr;
						for(i=0;i<h;i++) {
							for(j=0;j<ng;j++) {
								GLIK[i][j] += F[j] * LIKHI[use][i][j] / LIKHI[use][i][ng];
							}
						}
					}
				}
			}
			cout << "\rDone f chain ";
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
		if((iter+1)%thin==0) {
			out << (iter+1);
			for(i=0;i<ng;i++) {
				for(j=0;j<ng+1;j++) out << tab << A[use][i][j];
			}
			for(i=0;i<ng;i++) out << tab << R[use][i][0];
			out << tab << likelihood.LOG();
			/* Check */
			//mydouble newlik = known_source_lik6_composite(a[use],b[use],r[use]);
			out << tab << newlik.LOG();
			out << tab << logalpha.LOG();
			out << tab << move;
			out << endl;
		}
		if((current=clock())>next) {
			cout << "\rDone " << (iter+1) << " of " << niter << " iterations in " << (double)(current-start)/CLOCKS_PER_SEC << " s " << flush;
			next = current + (clock_t)CLOCKS_PER_SEC;
		}
	}
	cout << endl;
	out.close();
	o3.close();

	// Re-normalise the posterior probabilities of source
	for(i=0;i<h;i++) {
		for(j=0;j<ng;j++) {
			GLIK[i][j] /= (double)ctr;
		}
	}
	//Vector<int> jmax(h);
	//for(i=0;i<h;i++) {
	//	jmax[i] = 0;
	//	mydouble likmax = GLIK[i][0];
	//	for(j=1;j<ng;j++) {
	//		if(GLIK[i][j]>likmax) {
	//			jmax[i] = j;
	//			likmax = GLIK[i][j];
	//		}
	//	}
	//}
	//o2 << jmax[0];
	//for(i=1;i<h;i++) o2 << tab << jmax[i];
	//for(i=0;i<h;i++) o2 << tab << GLIK[i][jmax[i]].todouble();
	for(i=0;i<h;i++) {
		for(j=0;j<ng;j++) {
			if(!(i==0 && j==0)) o2 << tab;
			o2 << GLIK[i][j].todouble();
		}
	}
	o2 << endl;
	o2.close();
}
