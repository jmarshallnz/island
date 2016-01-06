#include "cluster.h"

int Cluster::multinom(Vector<double> &p, Random &ran) {
	double U = ran.U();
	int i;
	for(i=0;i<p.size();i++) {
		if((U-=p[i])<=0.0) break;
	}
	if(U>0.0) error("Problem in multinom");
	return i;
}

void Cluster::open_all(const char* filename) {
	init = true;
	// Read in the raw data file
	tsv TSV;
	TSV.read(filename);

	// Count the number of groups and loci
	nloc = TSV.data.ncols()-2;
	if(!(nloc>0))
		error("There must be a positive number of loci");
	ng = TSV.n_values()[nloc+1]-1;
	if(!(ng>=2)) error("There must be at least one target and two source pops");

	// Convert to integers
	Matrix<int> isolate(TSV.data.nrows(),TSV.data.ncols());
	int i,j;
	for(i=0;i<TSV.data.nrows();i++) {
		for(j=0;j<TSV.data.ncols();j++) {
			isolate[i][j] = atoi(TSV.data[i][j].c_str());
		}
	}

	// Do some preliminary passes
	int maxST = -1;			///< maximum ST
	Vector<int> maxallele(nloc,-1);	///< maximum allele number for each locus
	int nhuman = 0; 		///< number of human isolates
	size = Vector<int>(ng+1,0);	///< size[0:(ng-1)] is total in each source group, size[ng] is total in source groups
	for(i=0;i<isolate.nrows();i++) {
		if(!(isolate[i][0]>=0)) {
			cout << "Sequence " << i+1 << " column 1 reads " << TSV.data[i][0] << endl;
			error("Sequences should be labelled with a non-negative integer");
		}
		if(isolate[i][0]>10*isolate.nrows()) {
//			cout << "Sequence " << i+1 << " column 1 reads " << TSV.data[i][0] << endl;
//			cout << "WARNING: this integer (" << isolate[i][0] << ") is large - possible mistake?" << endl;
		}
		if(isolate[i][0]>maxST) maxST = isolate[i][0];
		for(j=1;j<=nloc;j++) {
			if(!(isolate[i][j]>=0)) {
				cout << "Sequence " << i+1 << " column " << j+1 << " reads " << TSV.data[i][j] << endl;
				error("Alleles must be given non-negative integer labels");
			}
			if(isolate[i][j]>10*isolate.nrows()) {
	//			cout << "Sequence " << i+1 << " column " << j+1 << " reads " << TSV.data[i][j] << endl;
//				cout << "WARNING: this integer (" << isolate[i][j] << ") is large - possible mistake?" << endl;
			}
			if(isolate[i][j]>maxallele[j-1]) maxallele[j-1] = isolate[i][j];
		}
		if(!(isolate[i][nloc+1]>=0 && isolate[i][nloc+1]<=ng)) {
			cout << "Sequence " << i+1 << " final column reads " << TSV.data[i][nloc+1] << endl;
			error("Groups must be numbered from 0 (target pop) to the number of source pops");
		}
		if(isolate[i][nloc+1]>0) {
			++size[isolate[i][nloc+1]-1];
			++size[ng];
		}
		else ++nhuman;
	}

	// Check that all sequence labels are unique and mutually exclusive
	// - this is a test for potential problems in the data file
	Matrix<int> aprofile(MIN(maxST,isolate.nrows()),nloc+1,-1);
	Vector<int> STwhere(maxST+1,-1);
	int NST = 0;
	for(i=0;i<isolate.nrows();i++) {
		const int lab = isolate[i][0];
		// If this ST has already been observed, make sure it is the same
		if(STwhere[lab]!=-1) {
			bool identical = true;
			for(j=0;j<nloc;j++) {
				if(!(aprofile[STwhere[lab]][j]==isolate[i][j+1])) {
					identical = false;
					break;
				}
			}
			if(!identical) {
				cout << "Sequence labelled " << lab << " has two different genotypes:" << endl;
				for(j=0;j<nloc;j++) cout << "\t" << aprofile[STwhere[lab]][j];
				cout << endl;
				for(j=0;j<nloc;j++) cout << "\t" << isolate[i][j+1];
				cout << endl;
				error("Sequence labels must be unique and mutually exclusive");
			}
		}
		// Otherwise, make sure it is different to all previously observed STs
		else {
			bool exclusive = true;
			int ii;
			for(ii=0;ii<NST;ii++) {
				bool identical = true;
				for(j=0;j<nloc;j++) {
					if(aprofile[ii][j]!=isolate[i][j+1]) {
						identical = false;
						break;
					}
				}
				if(identical) {
					exclusive = false;
					break;
				}
			}
			if(!exclusive) {
				for(j=0;j<maxST;j++) if(STwhere[j]==ii) break;
				if(j==maxST) error("Could not find duplicate sequence with different label");
				cout << "Sequence labelled " << lab << " on line " << i+1 << " has already been observed with label " << j << ":" << endl;
				for(j=0;j<nloc;j++) cout << "\t" << aprofile[ii][j];
				cout << endl;
				error("Sequence labels must be unique and mutually exclusive");
			}
			else {
				for(j=0;j<nloc;j++) aprofile[NST][j] = isolate[i][j+1];
				STwhere[lab] = NST;
				++NST;
			}
		}
	}
	if(NST>isolate.nrows()) error("Counted more sequence labels than sequences");
//////////////////////////////////////////////////////////////////////////////////////////////	
	// Create the matrix of human isolates
	human.resize(nhuman,nloc);	///< human loci for each case
	int ih;
	for(i=0,ih=0;i<isolate.nrows();i++) {
		if(isolate[i][nloc+1]==0) {
			for(j=0;j<nloc;j++) {
				human[ih][j] = isolate[i][j+1];
			}
			++ih;
		}
	}
	if(ih!=nhuman) error("Book-keeping problem identifying target isolates");
//////////////////////////////////////////////////////////////////////////////////////////////
	// Count the number of non-human isolates in each group
	/* Calculate the size and number of unique STs in each group */
	Matrix<int> niso(NST,ng+1,0);	///< number of each ST in each group, niso[,ng] is number of each ST in all groups
	size = Vector<int>(ng+1,0);	///< number of STs in each group, size[ng] is number of STs
	nST = Vector<int>(ng+1,0);	///< number of unique STs in each group, nST[ng] is number of unique STs
	for(i=0;i<isolate.nrows();i++) {
		const int lab = isolate[i][0];
		const int gp = isolate[i][nloc+1]-1;
		const int wh = STwhere[lab];
		if(gp>=0) {
			if(niso[wh][gp]==0) ++nST[gp];
			if(niso[wh][ng]==0) ++nST[ng];
			++niso[wh][gp];
			++niso[wh][ng];
			++size[gp];
			++size[ng];
		}
	}
	SIZE = Vector<mydouble>(ng+1);	///< number of isolates in each group, SIZE[ng] is number of STs
	for(i=0;i<=ng;i++) SIZE[i] = mydouble((double)size[i]);

	/* Allocate memory for allele counts */
	int k;
	acount.resize(ng+1,nloc);	///< counts of alleles in each group at each loci
	for(i=0;i<ng+1;i++) {
		for(j=0;j<nloc;j++) {
			acount[i][j].resize(maxallele[j]+1);
			for(k=0;k<=maxallele[j];k++) acount[i][j][k] = 0;
		}
	}
	/* Record the allelic profile for each unique ST in each group */
	MLST.resize(ng);		///< MLST[group,unique_st,loci] -> MLST[group] is profile of each unique ST
	freq.resize(ng);		///< freq[group,unique_st] proportion of unique_st in group
	abun.resize(ng);		///< abun[group,unique_st] number of unique_st in group
	FREQ.resize(ng);		///< FREQ[group,unique_st] same as freq but double
	ABUN.resize(ng);		///< ABUN[group,unique_st] same as abun but double
	for(i=0;i<ng;i++) {
		MLST[i].resize(nST[i],nloc);
		freq[i].resize(nST[i]);
		abun[i].resize(nST[i]);
		FREQ[i].resize(nST[i]);
		ABUN[i].resize(nST[i]);
	}
	Vector<int> ix(ng,0);		///< counter for each group - the allelic profile is copied in separately for each group
	for(i=0;i<NST;i++) { // for each unique ST
		int sc;
		for(sc=0;sc<ng;sc++) { // for each group
			const int ct = niso[i][sc]; // how many of this ST is in this group?
			if(ct>0) {
				for(j=0;j<nloc;j++) { // copy across the MLST profile
					MLST[sc][ix[sc]][j] = aprofile[i][j];
				}
				freq[sc][ix[sc]] = mydouble((double)ct/(double)size[sc]);
				abun[sc][ix[sc]] = mydouble((double)ct);
				FREQ[sc][ix[sc]] = (double)ct/(double)size[sc];
				ABUN[sc][ix[sc]] = (double)ct;
				for(j=0;j<nloc;j++) { // for each loci
					const int allele = aprofile[i][j];			// allele
					acount[sc][j][allele] += (double)ct/(double)size[sc];	// allele frequency
					//acount[ng][j][allele] += ct * (weight = size[ng]/ng/size[sc]) /size[ng]
					/* weighted */// acount[ng][j][allele] += (double)ct/(double)ng/(double)size[sc];
					/* unweighted */ acount[ng][j][allele] += (double)ct/(double)size[ng];
				}
				++ix[sc];
			}
		}
	}

	nalleles = Matrix<int>(ng+1,nloc,0);	///< number of alleles for each group and each loci. nalleles[ng] is total
	for(i=0;i<=ng;i++) {
		for(j=0;j<nloc;j++) {
			for(k=0;k<=maxallele[j];k++) {
				if(acount[i][j][k]>0) {
					++nalleles[i][j];
				}
			}
		}
	}
}

