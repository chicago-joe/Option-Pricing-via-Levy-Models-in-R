#include <ql/quantlib.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif
#include <ql/qldefines.hpp>
using namespace QuantLib;
using namespace boost;
#include <iostream>
#include "main.h"
using namespace std;

namespace {
	using namespace QuantLib;

	BOOST_AUTO_TEST_CASE(testGeometricBrownianMotion) {

		Real startingPrice = 20.16; //closing price for INTC on 12/7/2012
		Real mu = .2312; //one year historical annual return
		Volatility sigma = 0.2116; //one year historical volatility
		Size timeSteps = 255; //trading days in a year
		Time length = 1; //one year

		//instantiate Geometric Brownian Motion (GBM) stochastic process
		const boost::shared_ptr<StochasticProcess>& gbm =
			boost::shared_ptr<StochasticProcess > (new GeometricBrownianMotionProcess(startingPrice, mu, sigma));

		//generate sequence of normally distributed random numbers from uniform distribution using Box-Muller transformation
		BigInteger seed = SeedGenerator::instance().get();
		typedef BoxMullerGaussianRng<MersenneTwisterUniformRng> MersenneBoxMuller;
		MersenneTwisterUniformRng mersenneRng(seed);
		MersenneBoxMuller boxMullerRng(mersenneRng);
		RandomSequenceGenerator<MersenneBoxMuller> gsg(timeSteps, boxMullerRng);
		
        //generate simulated path of stock price
        PathGenerator<RandomSequenceGenerator<MersenneBoxMuller> > gbmPathGenerator(gbm, length, timeSteps, gsg, false);
		const Path& samplePath = gbmPathGenerator.next().value;
		
		//calculate simulated sample returns
		boost::function<Real (Real, Real)> calcLogReturns = [](Real x, Real y) {return std::log(y/x);};
		std::vector<Real> logReturns;
		Path::iterator samplePathBegin = samplePath.begin();
		Path::iterator samplePathEnd = samplePath.end();
		Path::iterator endMinusOne = std::prev(samplePathEnd);
		Path::iterator beginPlusOne = std::next(samplePathBegin);
		
		std::transform(samplePathBegin, endMinusOne, beginPlusOne,
				std::back_inserter(logReturns), calcLogReturns);		
		
		//calculate some general statistics
		GeneralStatistics statistics;

		//returns statistics
		statistics.addSequence(logReturns.begin(), logReturns.end());
		std::cout << boost::format("Standard deviation of simulated returns (Normal): %.4f") % 
				(statistics.standardDeviation() * std::sqrt(255)) << std::endl;

		//price statistics
		statistics.reset();
		statistics.addSequence(samplePath.begin(), samplePath.end());
		std::cout << boost::format("Price statistics: mean=%.2f, min=%.2f, max=%.2f") %
			statistics.mean() % statistics.min() % statistics.max() << std::endl;  

		//write simulated path to a file for charting
        std::ofstream gbmFile;
		gbmFile.open("/home/mick/Documents/blog/geometric-brownian-motion/gbm.dat", std::ios::out);
		for (Size i = 0; i < timeSteps; ++i) {
			gbmFile << boost::format("%d %.4f") % i % samplePath.at(i) << std::endl;
		}
		
		gbmFile.close();

	}
	#include <boost/math/distributions/students_t.hpp>
	
	Real studentTInverse(boost::math::students_t_distribution<> d, const Real& p) {
		return quantile(d,p);
	}

	BOOST_AUTO_TEST_CASE(testGeometricBrownianMotionStudentT) {

		Real startingPrice = 20.16; //closing price for INTC on 12/7/2012
		Real mu = .2312; //one year historical annual return
		Volatility sigma = 0.2116;
		Volatility scaledSigma = std::sqrt(sigma * sigma * 3/5); //one year historical volatility scaled by reciprocal of student-t variance (v/(v - 2)) 
		Size timeSteps = 255; //trading days in a year
		Time length = 1; //one year
		
		//instantiate Geometric Brownian Motion (GBM) stochastic process
		const boost::shared_ptr<StochasticProcess>& gbm =
			boost::shared_ptr<StochasticProcess > (new GeometricBrownianMotionProcess(startingPrice, mu, scaledSigma));

		//random sequence generator will use the Mersenne Twister algorithm
		BigInteger seed = SeedGenerator::instance().get();
		MersenneTwisterUniformRng mersenneRng(seed);
		RandomSequenceGenerator<MersenneTwisterUniformRng> rsg(timeSteps, mersenneRng);

		//instantiate Student T distribution from Boost math library
		boost::math::students_t_distribution<> studentT(5); //5 degrees of freedom - want fat tails!
		boost::function<Real (Real)> icd = boost::bind(studentTInverse, studentT, _1); 

		//samples random numbers from the Student T distribution		
		InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, 
            boost::function<Real (Real)> > invCumRsg(rsg, icd);

		//generates a single path
		PathGenerator<InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, 
            boost::function<Real (Real)> > > gbmPathGenerator(gbm, length, timeSteps, invCumRsg, false);

		const Path& samplePath = gbmPathGenerator.next().value;

		//calculate simulated sample returns
		boost::function<Real (Real, Real)> calcLogReturns = [](Real x, Real y) {return std::log(y/x);};
		std::vector<Real> logReturns;
		Path::iterator samplePathBegin = samplePath.begin();
		Path::iterator samplePathEnd = samplePath.end();
		Path::iterator endMinusOne = std::prev(samplePathEnd);
		Path::iterator beginPlusOne = std::next(samplePathBegin);
		
		std::transform(samplePathBegin, endMinusOne, beginPlusOne,
		    std::back_inserter(logReturns), calcLogReturns);		
		
		//calculate some general statistics
		GeneralStatistics statistics;

		//returns statistics
		statistics.addSequence(logReturns.begin(), logReturns.end());
		std::cout << boost::format("Standard deviation of simulated returns (Student-T): %.4f") % 
	        (statistics.standardDeviation() * std::sqrt(255)) << std::endl;

		//price statistics
		statistics.reset();
		statistics.addSequence(samplePath.begin(), samplePath.end());
		std::cout << boost::format("Price statistics: mean=%.2f, min=%.2f, max=%.2f") %
			statistics.mean() % statistics.min() % statistics.max() << std::endl;  

		//write results to a file 
		std::ofstream gbmFile;
		gbmFile.open("/home/mick/Documents/blog/geometric-brownian-motion/gbm-student.dat", std::ios::out);
		for (Size i = 0; i < timeSteps; ++i) {
		    gbmFile << boost::format("%d %.4f") % i % samplePath.at(i) << std::endl;
		}
		
		gbmFile.close();

	}

    BOOST_AUTO_TEST_CASE(testGeometricBrownianMotionStudentTLowDiscrepancyRSG) {

		Real startingPrice = 20.16; //closing price for INTC on 12/7/2012
		Real mu = .2312; //one year historical annual return
		Volatility sigma = 0.2116;
		Volatility scaledSigma = std::sqrt(sigma * sigma * 3/5); //one year historical volatility scaled by reciprocal of student-t variance (v/(v - 2)) 
		Size timeSteps = 255; //trading days in a year
		Time length = 1; //one year
		
		//instantiate Geometric Brownian Motion (GBM) stochastic process
		const boost::shared_ptr<StochasticProcess>& gbm =
			boost::shared_ptr<StochasticProcess > (new GeometricBrownianMotionProcess(startingPrice, mu, scaledSigma));

		//random sequence generator will use the Mersenne Twister algorithm
		BigInteger seed = SeedGenerator::instance().get();
        HaltonRsg rsg(timeSteps);
        
		//instantiate Student T distribution from Boost math library
		boost::math::students_t_distribution<> studentT(5); //5 degrees of freedom - want fat tails!
		boost::function<Real (Real)> icd = boost::bind(studentTInverse, studentT, _1); 

		//samples random numbers from the Student T distribution		
		InverseCumulativeRsg<HaltonRsg, boost::function<Real (Real)> > invCumRsg(rsg, icd);

		//generates a single path
		PathGenerator<InverseCumulativeRsg<HaltonRsg, boost::function<Real (Real)> > > gbmPathGenerator(gbm, length, timeSteps, invCumRsg, false);

		const Path& samplePath = gbmPathGenerator.next().value;
		
		//calculate simulated sample returns
		boost::function<Real (Real, Real)> calcLogReturns = [](Real x, Real y) {return std::log(y/x);};
		std::vector<Real> logReturns;
		Path::iterator samplePathBegin = samplePath.begin();
		Path::iterator samplePathEnd = samplePath.end();
		Path::iterator endMinusOne = std::prev(samplePathEnd);
		Path::iterator beginPlusOne = std::next(samplePathBegin);
		
		std::transform(samplePathBegin, endMinusOne, beginPlusOne,
				std::back_inserter(logReturns), calcLogReturns);		
		
		//calculate some general statistics
		GeneralStatistics statistics;

		//returns statistics
		statistics.addSequence(logReturns.begin(), logReturns.end());
		std::cout << boost::format("Standard deviation of simulated returns(Student-T / Halton): %.4f") % 
				(statistics.standardDeviation() * std::sqrt(255)) << std::endl;

		//price statistics
		statistics.reset();
		statistics.addSequence(samplePath.begin(), samplePath.end());
		std::cout << boost::format("Price statistics: mean=%.2f, min=%.2f, max=%.2f") %
			statistics.mean() % statistics.min() % statistics.max() << std::endl;  
		
		//write results to a file 
		std::ofstream gbmFile;
		gbmFile.open("/home/mick/Documents/blog/geometric-brownian-motion/gbm-student-halton.dat", std::ios::out);
		for (Size i = 0; i < timeSteps; ++i) {
			gbmFile << boost::format("%d %.4f") % i % samplePath.at(i) << std::endl;
		}
		
		gbmFile.close();
	}
}





int main()
{
	BigInteger seed = 12324;
	MersenneTwisterUniformRng	generator(seed);
	InverseCumulativeRng<MersenneTwisterUniformRng, MoroInverseCumulativeNormal> InvGauss(generator);
	// cout << InvGauss.next().value << endl;
	//for (int i = 0; i < 100; ++i) 
	//{
	//	cout << InvGauss.next().value << endl;
	//}
	//modifiedBesselFunction_i()
	

	Date refDate = Date(04, Jan, 2019);
	Rate riskFreeRate = 0.0321;
	Rate dividendRate = 0.0128;
	Real spot = 52.0;
	Rate vol = 0.2144;
	Calendar cal = TARGET();
	DayCounter dc = ActualActual();

	boost::shared_ptr<YieldTermStructure> rdStruct(new FlatForward(refDate, riskFreeRate, dc));
	boost::shared_ptr<YieldTermStructure> rqStruct(new FlatForward(refDate, dividendRate, dc));
	Handle<YieldTermStructure> rdHandle(rdStruct);
	Handle<YieldTermStructure> rqHandle(rqStruct);

	boost::shared_ptr<SimpleQuote> spotQuote(new SimpleQuote(spot));
	Handle<Quote> spotHandle(spotQuote);

	boost::shared_ptr<BlackVolTermStructure> volQuote(new BlackConstantVol(refDate, cal, vol, dc));
	Handle<BlackVolTermStructure> volHandle(volQuote);

	Real v0 = 0.12, kappa = 1.2, theta = 0.08, sigma = 0.05, rho = -0.6;
	Real lambda = 0.25, nu = 0.0, delta = 0.30;

	boost::shared_ptr<BatesProcess> batesProcess(new BatesProcess(rdHandle, rqHandle, spotHandle, v0,
		kappa, theta, sigma, rho, lambda, nu, delta, HestonProcess::PartialTruncation));

	Time dt = 0.10, t = 0.0;
	Array dw(4), x(2);	
	x[0] = spotQuote->value();		// x is the 2-dimensional process
	x[1] = v0;
	Size numVals = 10;
	for (Size j = 1; j <= numVals; ++j) {
		dw[0] = InvGauss.next().value;
		dw[1] = InvGauss.next().value;
		dw[2] = InvGauss.next().value;
		dw[3] = InvGauss.next().value;

		x = batesProcess->evolve(t, x, dt, dw);
		std::cout << "Time: " << t + dt << ", S(t): " << x[0] << ", V(t): " << x[1] << std::endl;
		t += dt;
	}
}
