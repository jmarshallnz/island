#include "cluster.h"

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

mydouble Cluster::known_source_lik4_composite(Matrix<double> &a, Matrix< Vector<double> > &b) {
	int i,j,l;
	mydouble lik = 1.0;
	/* Cycle through each unique ST in each group, taking account of abundance of the STs */
	for(i=0;i<ng;i++) {
		for(j=0;j<nST[i];j++) {
			double ncopies = ceil(freq[i][j].todouble()*(double)size[i]-0.5);
			for(l=0;l<nloc;l++) {
				int allele = MLST[i][j][l];
				double num = acount[ng][l][allele] * (double)size[ng];
				if(num<1.1) { // new allele (allow some rounding error)
					lik *= mydouble(a[i][ng])^ncopies;
				}
				else {	// previously observed
//					double bk = b[i][l][allele] + a[i][i]/(double)(size[i]-1)*(acount[i][l][allele] - 1.0);
					double ac = acount[i][l][allele];
					double bk = b[i][l][allele] - a[i][i]*ac + a[i][i]*(ac*(double)size[i]-1.0)/(double)(size[i]-1);
					if(fabs(bk)<1.0e-7) {
						bk = 0.0;
						warning("Factor of zero in likelihood");
					}
					lik *= mydouble(bk)^ncopies;
				}
			}
		}
	}
	return lik;
}

mydouble Cluster::likHi4(const int id, const int i, Matrix<double> &a, Matrix< Vector<double> > &b) {
//	return mydouble(1.0);
	int k;
	mydouble l_i(1.0);
	/* Assume free recombination between loci */
	for(k=0;k<nloc;k++) {
		int human_allele = human[id][k];
		if(human_allele>=acount[ng][k].size() || acount[ng][k][human_allele]==0) {	// unique
			l_i *= mydouble(a[i][ng]);
		}
		else {	// not unique
			l_i *= mydouble(b[i][k][human_allele]);
		}
	}
	return l_i;
//	return l_i^0.1;
}

