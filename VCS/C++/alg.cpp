
#include <cmath>
#include<vector>
#include <chrono>
#include <random>
#include <iostream>
using namespace std;

uniform_real_distribution<double> undistribution(0, 1);
normal_distribution<double> distribution(0, 1);
random_device rdn{};
mt19937 gen{ rdn() };
default_random_engine generator;


int main(int argc, char* argv[])
{
	alg1(1,15,-5,100,0.5);
	alg2(1,2,3,alg1(1, 15, -5, 1, 0.5),alg1(1, 15, -5, 3, 0.5),0.5, alg1(1, 15, -5, 1, 0.5), alg1(1, 15, -5, 3, 0.5),1);
	system("pause");
}
double alg1(double mu,double alpha,double beta,double t,double sigma) {

	double rand;
	
	double  G1, gama, Z, zeta, zt, G2, Xt;
	gama = pow(pow(alpha, 2) - pow(beta, 2), 0.5);
	for (size_t i = 1; i <= t; i++)
	{
		G1 = distribution(gen);
		Z = pow(G1, 2) / gama;
		zeta = 1.0 / gama * (sigma*i+0.5*Z-pow(sigma*i*Z+pow(Z,2)/4,0.5));
		rand = undistribution(gen);
		if (rand<sigma*i/(sigma*i+gama* zeta))
		{
			zt = zeta;
		}
		else
		{
			zt = pow(sigma, 2) *pow(i, 2) / (pow(gama, 2) * zeta);
		}
		G2 = distribution(gen);
		Xt = mu * i + beta * zt + pow(zt, 0.5)*G2;
		//mu = Xt;
		//cout << Xt << endl;
	}
	
	return zt;
}

double alg2(double ti,double tj,double tk,double zti,double ztk,double sigma,double Wzti,double Wztk,double mu) {
	double ztj;
	double lambd,phi,Q,G1j,s1,s,Uj;
	lambd = pow(sigma, 2)*pow(tk - tj, 2) / (ztk - zti);
	phi = (tk - tj) / (tj - ti);
	G1j = distribution(gen);
	Q = pow(G1j, 2);
	s1 = phi + pow(phi, 2)*Q / (2.0*lambd) - phi / (2.0*lambd)*pow(4 * phi*lambd*Q + pow(phi, 2)*pow(Q, 2), 0.5);
	Uj = undistribution(gen);
	if (Uj<phi*(1+s1)/(1+phi)/(phi+s1))
	{
		s = s1;
	}
	else {
		s = pow(phi, 2) / s1;
	}
	ztj = zti + (ztk - zti) / (1 + s);
	double m = ((ztk - ztj)*Wzti + (ztj - zti)*Wztk) / (ztk - zti);
	double sigma2 = (ztj - zti)*(ztk - ztj) / (ztk - zti);
	double G2j = distribution(gen);
	double Wztj = m + pow(sigma2, 0.5)*G2j;
	double Xtj = mu * tj + Wztj;
	return Xtj;
}