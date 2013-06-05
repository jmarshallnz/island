#include "cluster.h"

double Cluster::sigma(const double x) {
	if(x<0.05) return 0.608*x;
	else if(x>0.95) return 0.608*(1-x);
	return 0.03;
}

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

// Assumes F is correctly sized
void Cluster::calc_F(Vector<mydouble> &f, Vector<mydouble> &F) {
	int i;
	const int n = f.size()-1;
	f[n] = 0.0;
	for(i=0;i<n;i++) f[n] += f[i];
	for(i=0;i<n;i++) F[i] = (f[i]/f[n]);
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

