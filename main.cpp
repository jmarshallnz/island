/*	Simulation study for "source". 

*/
#include "myutils.h"
#include "cluster.h"

using namespace myutils;

Vector<double> simulate_data(Random &ran)
{
	/* simulate some data */
	double bt = 0.5;
	double st = 1;
	int     T = 200;
	Vector<double> y(T+1,0);
	y[0] = 0;
	for (int t = 1; t <= T; t++)
		y[t] = ran.normal(bt*y[t-1],st);
	return y;
}

void gibbs_ar1(Vector<double> &y, double &b, double &s2, Random &ran)
{
	/*
	 Fits y[t] = by[t-1] + e where e=N(0,s2)

	 Then b ~ N(bhat, shat2/sxx) and
              s2 ~ G(syy-bhat*sxx/(T-1))

	 where bhat = sxy/sxx and shat2 = (syy - bhat*sxx)/(T-1)

	 Note that we don't actually use b, or s2 - we estimate them directly
	 from the y's. i.e. this is equivalent to drawing from the posterior p(b|y[0..T])
	 then p(s2|b,y[0..T]), then p(y[0]|b,s2,y[1]).

	 */
	const int T = y.size() - 1;

	/* priors */
	const double b0 = 0;  	///< autocorrelation mean
	const double pb = 0.001;///< autocorrelation precision

	const double v0 = 0.01; ///< gamma a
	const double vs = 0.01; ///< gamma b

	/* posterior hyperparams */

	double sxx = 0;
	double sxy = 0;
	for (int t = 0; t < T; t++) {
		sxx += y[t]*y[t];
		sxy += y[t] * y[t+1];
	}
	double syy = sxx + y[T]*y[T] - y[0]*y[0];

	/* estimates from OLS regression */
	double bhat = sxy / sxx;
	double s2f = (syy - bhat*sxy)/(T-1);

	/*
	   posterior for beta is
		N(sxy/(pb + sxx), s^2/(pb + sxx))

	   posterior for sigma is
		IG(v0 + T/2, vs + rss/2)
	 */

	double bf = sxy / (pb + sxx);
	double bs = sqrt(s2f / (pb + sxx));

	/* 1: sample from p(b|y[0..T]) */
	b = bf + bs*ran.Z();

	/* 2: sample from p(s2|b,y[0..T]) */

	/* from sum((y - bx)^2) */
	double rss = syy - 2*b*sxy + b*b*sxx;
	s2 = ran.gamma(1/(vs+rss/2), v0 + T/2);

	/* 3: sample from p(y[0]|b,s2,y[1])  */
	y[0] = y[1]/b + sqrt(s2)*ran.Z();
}

int main(const int argc, const char* argv[]) {
	if(argc < 4) error("SYNTAX: animals_file humans_file design_file out_file [niter thinning alpha seed]");
	const char* animals_file = argv[1];
	const char* humans_file = argv[2];
	const char* design_file = argv[3];
	const char* out_file = argv[4];
	const int niter = argc > 5 ? atoi(argv[5]) : 1000;
	const int thinning = argc > 6 ? atoi(argv[6]) : 10;
	const double alpha = argc > 7 ? atof(argv[7]) : 1.0;
	Random ran;
	ran.setseed(argc > 8 ? atoi(argv[8]) : -3);

	cout << argc << " arguments provided" << endl;
	cout << "Iterations set to " << niter << endl;
	cout << "Thinning set to " << thinning << endl;
	cout << "Alpha set to " << alpha << endl;
	cout << "Seed set to " << ran.getseed() << endl;

	/*
	Vector<double> y = simulate_data(ran);
	int T = y.size() - 1;

	int burnin = 1000;
	int iters = 10000;

	ofstream out(out_file);
	for (int i = 0; i < burnin + iters; i++) {

		double b, s2;
		gibbs_ar1(y, b, s2, ran);

		if (i >= burnin) {
			out << b << '\t';
			out << s2 << '\t';
			out << y[0] << '\n';
		}
	}
	*/

	Cluster clust;
	clust.open_all(animals_file, humans_file, design_file);
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

