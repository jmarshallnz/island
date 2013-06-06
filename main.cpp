/*	Simulation study for "source". 

*/
#include "myutils.h"
#include "cluster.h"

using namespace myutils;

Vector<double> simulate_data(Random &ran)
{
	/* simulate some data */
	double at = 5.0;
	double bt = 0.5;
	double st = 1;
	int     T = 200;
	Vector<double> y(T+1,0);
	y[0] = at / (1-bt);
	for (int t = 1; t <= T; t++)
		y[t] = ran.normal(at + bt*y[t-1],st);
	return y;
}

void gibbs_ar1(Vector<double> &y, double &a, double &b, double &s2, Random &ran)
{
	/*
	 Fits y[t] = a + by[t-1] + e where e=N(0,s2)

	 Note that we don't actually use a,b, or s2 - we estimate them directly
	 from the y's. i.e. this is equivalent to drawing from the posterior p(a,b|y[0..T])
	 then p(s2|a,b,y[0..T]), then p(y[0]|a,b,s2,y[1]).

	 */
	const int T = y.size() - 1;

	/* priors */
	const double a0 = 0;  ///< constant
	const double va = 1000;
	const double b0 = 0;  ///< autocorrelation
	const double vb = 1000;
	const double v0 = 0.01; ///< gamma a
	const double vs = 0.01; ///< gamma b

	/* posterior hyperparams */

	double sx = 0;
	double sxx = 0;
	double sxy = 0;
	for (int t = 0; t < T; t++) {
		sx += y[t];
		sxx += y[t]*y[t];
		sxy += y[t] * y[t+1];
	}
	double sy = sx + y[T] - y[0];
	double syy = sxx + y[T]*y[T] - y[0]*y[0];

	/*
	   posterior for alpha,beta is
		N(theta1, s2*C1)
	   where theta1 = C1(inv(C0)theta + X'y)
	   and C1 = inv(inv(C0) + X'X)
	 */

	/* C1 */
	double m11 = (1/va + T);
	double m00 = (1/vb + sxx);
	double m01 = -sx;
	double det = m00*m11 - m01*m01;
	m00 /= det;
	m11 /= det;
	m01 /= det;

	/* theta1 */
	double t1 = a0 / va + sy;
	double t2 = b0 / vb + sxy;
	double af = (m00*t1 + m01*t2);
	double bf = (m01*t1 + m11*t2);

	/* Cholesky decomp of C1 */
	double l00 = sqrt(m00);
	double l10 = m01/l00;
	double l11 = sqrt(m11 - l10*l10);

	/* estimate s2f */
	double bhat  = (T*sxy - sx*sy) / (T*sxx - sx*sx);
	double s2f = (T*syy - sy*sy - bhat*bhat*(T*sxx - sx*sx))/T/(T-2);

	/* 1: sample from p(a,b|y[0..T]) */
	double z0 = ran.Z(), z1 = ran.Z();
	a = af + sqrt(s2f)*(l00*z0);
	b = bf + sqrt(s2f)*(l10*z0 + l11*z1);

	/* 2: sample from p(s2|a,b,y[0..T]) */

	/* from sum((y - (a + bx)))^2 */
	double rss = syy - 2*a*sy - 2*b*sxy + T*a*a + 2*a*b*sx + b*b*sxx;
	s2 = ran.gamma(2/(vs+rss), (v0 + T)/2);

	/* 3: sample from p(y[0]|a,b,s2,y[1])  */
	y[0] = a + b*y[1] + sqrt(s2)*ran.Z();
}

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

	/*
	Vector<double> y = simulate_data(ran);
	int T = y.size() - 1;

	int burnin = 1000;
	int iters = 10000;

	ofstream out(out_file);
	for (int i = 0; i < burnin + iters; i++) {

		double a, b, s2;
		gibbs_ar1(y, a, b, s2, ran);

		if (i >= burnin) {
			out << a << '\t';
			out << b << '\t';
			out << s2 << '\t';
			out << y[0] << '\n';
		}
	}
	*/

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

