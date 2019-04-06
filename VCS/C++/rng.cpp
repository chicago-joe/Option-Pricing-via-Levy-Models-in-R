#include "statistics.h"
#include <iostream>
#include <vector>
using namespace std;

int main(int argc, char **argv) {

	// Create the Standard Normal Distribution and random draw vectors
	StandardNormalDistribution snd;
	std::vector<double> uniform_draws(20, 0.0);
	std::vector<double> normal_draws(20, 0.0);

	// Simple random number generation method based on RAND
	for (int i = 0; i < uniform_draws.size(); i++) {
		uniform_draws[i] = rand() / static_cast<double>(RAND_MAX);
	}

	// Create standard normal random draws
	// Notice that the uniform draws are unaffected. We have separated
	// out the uniform creation from the normal draw creation, which
	// will allow us to create sophisticated random number generators
	// without interfering with the statistical classes
	snd.random_draws(uniform_draws, normal_draws);

	// Output the values of the standard normal random draws
	for (int i = 0; i < normal_draws.size(); i++) {
		std::cout << normal_draws[i] << std::endl;
	}

	return 0;
}