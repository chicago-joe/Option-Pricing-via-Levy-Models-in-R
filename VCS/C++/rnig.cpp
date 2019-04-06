#include <iostream>
#include <math.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <string>
#include <sstream>

using namespace std;

long double *Box_Muller(int n, long double u[]);
long double *Beasley_Springer_Moro(int n, long double u[]);
long double Anderson_Darling(int n, long double X[]);
long double phi(long double x);


int main() {
	int n = 2000;

	// Mersenne Twister
	long double Mersenne[n];
	random_device rd;
	mt19937 e2(1);
	uniform_real_distribution<long double> dist(0, 1);
	for (int i = 0; i < n; i++) {
		Mersenne[i] = dist(e2);
	}

	// Poor LCG
	long double x[n];
	x[0] = 1.0;
	for (int i = 1; i < n; i++) {
		x[i] = fmod(1229.0*x[i - 1] + 1.0, 2048.0);
	}
	long double u[n];
	for (int i = 0; i < n; i++) {
		u[i] = x[i] / 2048.0;
		//cout << "u: " <<  u[i] << endl;
	}



	// Print Anderson Statistic for Mersenne 6a
	//long double *result = new long double[n];
	//result = Box_Muller(n,Mersenne);
	//Anderson_Darling(n,result);

	// Print Anderson Statistic for poor LCG 6b
	//long double *result1 = new long double[n];
	//result1 = Box_Muller(n,u);
	//Anderson_Darling(n,result1);

	// Print Anderson Statistic for Mersenne 7
	//long double *result2 = new long double[n];
	//result2 = Beasley_Springer_Moro(n, Mersenne);
	//Anderson_Darling(n,result2);

	// Print Anderson Statistic for poor LCG 7
	long double *result3 = new long double[n];
	result3 = Beasley_Springer_Moro(n, u);
	Anderson_Darling(n, result3);




	return 0;
}

long double *Box_Muller(int n, long double u[]) {
	long double *X = new long double[n];
	long double Y[n];
	long double R_2[n];
	long double theta[n];
	for (int i = 0; i < n; i++) {
		R_2[i] = -2.0*log(u[i]);
		theta[i] = 2.0*M_PI*u[i + 1];
	}
	for (int i = 0; i < n; i++) {
		X[i] = sqrt(-2.0*log(u[i]))*cos(2.0*M_PI*u[i + 1]);
		Y[i] = sqrt(-2.0*log(u[i]))*sin(2.0*M_PI*u[i + 1]);
	}
	return X;
}

long double *Beasley_Springer_Moro(int n, long double u[]) {
	long double y[n];
	long double r[n];
	long double *x = new long double[n];
	// Constants needed for algo
	long double a_0 = 2.50662823884;        long double b_0 = -8.47351093090;
	long double a_1 = -18.61500062529;   long double b_1 = 23.08336743743;
	long double a_2 = 41.39119773534;    long double b_2 = -21.06224101826;
	long double a_3 = -25.44106049637;   long double b_3 = 3.13082909833;

	long double c_0 = 0.3374754822726147; long double c_5 = 0.0003951896511919;
	long double c_1 = 0.9761690190917186; long double c_6 = 0.0000321767881768;
	long double c_2 = 0.1607979714918209; long double c_7 = 0.0000002888167364;
	long double c_3 = 0.0276438810333863; long double c_8 = 0.0000003960315187;
	long double c_4 = 0.0038405729373609;

	for (int i = 0; i < n; i++) {
		y[i] = u[i] - 0.5;
		if (fabs(y[i]) < 0.42) {
			r[i] = y[i] * y[i];
			x[i] = y[i] * (((a_3 * r[i] + a_2)*r[i] + a_1)*r[i] + a_0) / ((((b_3 * r[i] + b_2)*r[i] + b_1)*r[i] + b_0)*r[i] + 1.0);
		}
		else {
			r[i] = u[i];
			if (y[i] > 0) {
				r[i] = 1 - u[i];
			}
			r[i] = log(-log(r[i]));
			x[i] = c_0 + r[i] * (c_1 + r[i] * (c_2 + r[i] * (c_3 + r[i] * (c_4 + r[i] * (c_5 + r[i] * (c_6 + r[i] * (c_7 + r[i] * c_8)))))));
			if (y[i] < 0) {
				x[i] = -x[i];
			}
		}
	}
	ofstream myfile("bb.txt");
	for (int j = 0; j < n; j++) {
		myfile << x[j] << ";" << endl;
	}
	return x;
}

long double phi(long double x) {
	return 0.5 * erfc(-x * M_SQRT1_2);
}

long double Anderson_Darling(int n, long double X[]) {
	sort(X, X + n);
	// Find the mean of X
	long double X_avg = 0.0;
	long double sum = 0.0;
	for (int i = 0; i < n; i++) {
		sum += X[i];
	}
	X_avg = ((long double)sum) / n;


	// Find the variance of X
	long double X_sig = 0.0;
	for (int i = 0; i < n; i++) {
		X_sig += (X[i] - X_avg)*(X[i] - X_avg);
	}
	X_sig /= (n - 1);


	// The values X_i are standardized to create new values Y_i
	long double Y[n];
	for (int i = 0; i < n; i++) {
		Y[i] = (X[i] - X_avg) / (sqrt(X_sig));
		//cout << Y[i] << endl;
	}

	// With a standard normal CDF, we calculate the Anderson_Darling Statistic
	long double A = -n;
	for (int i = 0; i < n; i++) {
		A += -1.0 / (long double)n *(2 * (i + 1) - 1)*(log(phi(Y[i])) + log(1 - phi(Y[n - 1 - i])));
	}
	cout << A << endl;
}