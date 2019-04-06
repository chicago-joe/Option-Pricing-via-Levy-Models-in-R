#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include <math.h>
//depend on http://www.quantopia.net/inverse-normal-cdf/
double get_uniform()
{
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<double> distribution(0.0, 1.0);
	return(distribution(generator));
}
double a1 = 2.50662823884;
double a2 = -18.61500062529;
double a3 = 41.39119773534;
double a4 = -25.44106049637;

double b1 = -8.47351093090;
double b2 = 23.08336743743;
double b3 = -21.06224101826;
double b4 = 3.13082909833;

double c1 = 0.3374754822726147;
double c2 = 0.9761690190917186;
double c3 = 0.1607979714918209;
double c4 = 0.0276438810333863;
double c5 = 0.0038405729373609;
double c6 = 0.0003951896511919;
double c7 = 0.0000321767881768;
double c8 = 0.0000002888167364;
double c9 = 0.0000003960315187;

double get_NIG()
{
	double x = get_uniform();
	double y = x - 0.5;
	double output = 0;
	double z = 0;
	if (0.08 < x&&x < 0.92)
	{
		z = y * y;
		output = y * (a1 + a2 * z + a3 * z*z + a4 * z*z*z) / (1 + b1 * z + b2 * z*z + b3 * z*z*z + b4 * z*z*z*z);
	}
	else
	{
		double k = 0;
		if (y <= 0)
		{
			z = x;
			k = log(-log(z));
			output = -(c1 + c2 * pow(k, 1) + c3 * pow(k, 2) + c4 * pow(k, 3) + c5 * pow(k, 4) + c6 * pow(k, 5) + c7 * pow(k, 6) + c8 * pow(k, 7) + c9 * pow(k, 8));
		}
		else
		{
			z = 1 - x;
			k = log(-log(z));
			output = c1 + c2 * pow(k, 1) + c3 * pow(k, 2) + c4 * pow(k, 3) + c5 * pow(k, 4) + c6 * pow(k, 5) + c7 * pow(k, 6) + c8 * pow(k, 7) + c9 * pow(k, 8);
		}
	}
	return output;
}
int main()
{
	double temp = get_NIG();
	std::cout << temp << std::endl;
	system("pause");
	return 0;
}







