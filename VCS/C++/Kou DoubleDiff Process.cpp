#include <ql/qldefines.hpp>
#include <ql/quantlib.hpp>
#include <test-suite/utilities.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif

using namespace QuantLib;
#include <ios>
#include<ostream>
#include <boost/timer.hpp>
using namespace boost;
#include <iostream>
#include <iomanip>
using namespace std;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {
	Integer sessionID() { return 7; }
}
#endif

namespace {
	Real getCalibrationError(
		std::vector<ext::shared_ptr<BlackCalibrationHelper> > & options) {
		Real sse = 0;
		for (Size i = 0; i < options.size(); ++i) {
			const Real diff = options[i]->calibrationError()*100.0;
			sse += diff * diff;
		}
		return sse;
};
/*
	struct CalibrationMarketData {
		Handle<Quote> s0;
		Handle<YieldTermStructure> riskFreeTS, dividendYield;
		std::vector<ext::shared_ptr<BlackCalibrationHelper> > options;
};
*/

void BatesAnalyticVsBlack() 
{
	cout << "Testing analytic Bates engine against Black formula..." << endl;
	SavedSettings backup;
	Date settlementDate = Date::todaysDate();
	Settings::instance().evaluationDate() = settlementDate;
	DayCounter dayCounter = ActualActual();
	Date exerciseDate = settlementDate + 6 * Months;

	ext::shared_ptr<StrikedTypePayoff> payoff(
		new PlainVanillaPayoff(Option::Put, 30));
	ext::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));

	Handle<YieldTermStructure> riskFreeTS(flatRate(0.1, dayCounter));
	Handle<YieldTermStructure> dividendTS(flatRate(0.04, dayCounter));
	Handle<Quote> s0(ext::shared_ptr<Quote>(new SimpleQuote(32.0)));

	Real yearFraction = dayCounter.yearFraction(settlementDate, exerciseDate);
	Real forwardPrice = s0->value()*std::exp((0.1 - 0.04)*yearFraction);
	Real expected = blackFormula(payoff->optionType(), payoff->strike(),
		forwardPrice, std::sqrt(0.05*yearFraction)) *
		std::exp(-0.1*yearFraction);

	const Real v0 = 0.05;
	const Real kappa = 5.0;
	const Real theta = 0.05;
	const Real sigma = 1.0e-4;
	const Real rho = 0.0;
	const Real lambda = 0.0001;
	const Real nu = 0.0;
	const Real delta = 0.0001;

	VanillaOption option(payoff, exercise);

	ext::shared_ptr<BatesProcess> process(
		new BatesProcess(riskFreeTS, dividendTS, s0, v0,
			kappa, theta, sigma, rho, lambda, nu, delta));

	ext::shared_ptr<PricingEngine> engine(new BatesEngine(
		ext::make_shared<BatesModel>(process), 64));

	option.setPricingEngine(engine);
	Real calculated = option.NPV();

	Real tolerance = 2.0e-7;
	Real error = std::fabs(calculated - expected);
	if (error > tolerance) {
		cout << "Failed to reproduce Black price with BatesEngine"
			<< std::fixed
			<< "\n    calculated: " << calculated
			<< "\n    expected:   " << expected
			<< std::scientific
			<< "\n    error:      " << error << endl;
	}

	engine = ext::shared_ptr<PricingEngine>(new BatesDetJumpEngine(
		ext::make_shared<BatesDetJumpModel>(
			process, 1.0, 0.0001), 64));

	option.setPricingEngine(engine);
	calculated = option.NPV();

	error = std::fabs(calculated - expected);
	if (error > tolerance) {
		cout << "failed to reproduce Black price with " \
			"BatesDetJumpEngine"
			<< std::fixed
			<< "\n    calculated: " << calculated
			<< "\n    expected:   " << expected
			<< std::scientific
			<< "\n    error:      " << error << endl;
	}

	engine = ext::shared_ptr<PricingEngine>(new BatesDoubleExpEngine(
		ext::make_shared<BatesDoubleExpModel>(
			process, 0.0001, 0.0001, 0.0001), 64));

	option.setPricingEngine(engine);
	calculated = option.NPV();

	error = std::fabs(calculated - expected);
	if (error > tolerance) {
		cout << "failed to reproduce Black price with BatesDoubleExpEngine"
			<< std::fixed
			<< "\n    calculated: " << calculated
			<< "\n    expected:   " << expected
			<< std::scientific
			<< "\n    error:      " << error << endl;
	}

	engine = ext::shared_ptr<PricingEngine>(new BatesDoubleExpDetJumpEngine(
		ext::make_shared<BatesDoubleExpDetJumpModel>(

			process, 0.0001, 0.0001, 0.0001, 0.5, 1.0, 0.0001), 64));

	option.setPricingEngine(engine);
	calculated = option.NPV();

	error = std::fabs(calculated - expected);
	if (error > tolerance) {
		cout << "failed to reproduce Black price with " \
			"BatesDoubleExpDetJumpEngine"
			<< std::fixed
			<< "\n    calculated: " << calculated
			<< "\n    expected:   " << expected
			<< std::scientific
			<< "\n    error:      " << error << endl;
	}
}

void Bates_vs_Merton76JumpDiffusion() {

	cout << "Testing analytic Bates engine against Merton-76 engine..." << endl;

	SavedSettings backup;

	Date settlementDate = Date::todaysDate();
	Settings::instance().evaluationDate() = settlementDate;

	DayCounter dayCounter = ActualActual();

	ext::shared_ptr<StrikedTypePayoff> payoff(
		new PlainVanillaPayoff(Option::Put, 95));

	Handle<YieldTermStructure> riskFreeTS(flatRate(0.1, dayCounter));
	Handle<YieldTermStructure> dividendTS(flatRate(0.04, dayCounter));
	Handle<Quote> s0(ext::shared_ptr<Quote>(new SimpleQuote(100)));

	//Handle<YieldTermStructure> BlackVolTermStructure((settlementDate, vol, dayCounter);
	Real v0 = 0.0433;
	// FLOATING_POINT_EXCEPTION
	ext::shared_ptr<SimpleQuote> vol(new SimpleQuote(std::sqrt(v0)));
	ext::shared_ptr<BlackVolTermStructure> volTS = flatVol(settlementDate, vol, dayCounter);

	const Real kappa = 0.5;
	const Real theta = v0;
	const Real sigma = 1.0e-4;
	const Real rho = 0.0;

	ext::shared_ptr<SimpleQuote> jumpIntensity(new SimpleQuote(2));
	ext::shared_ptr<SimpleQuote> meanLogJump(new SimpleQuote(-0.2));
	ext::shared_ptr<SimpleQuote> jumpVol(new SimpleQuote(0.2));

	ext::shared_ptr<BatesProcess> batesProcess(new BatesProcess(
		riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho,
		jumpIntensity->value(), meanLogJump->value(), jumpVol->value()));

	ext::shared_ptr<Merton76Process> mertonProcess(
		new Merton76Process(s0, dividendTS, riskFreeTS,
			Handle<BlackVolTermStructure>(volTS),
			Handle<Quote>(jumpIntensity),
			Handle<Quote>(meanLogJump),
			Handle<Quote>(jumpVol)));

	ext::shared_ptr<PricingEngine> batesEngine(new BatesEngine(
		ext::make_shared<BatesModel>(batesProcess), 160));

	const Real mcTol = 0.1;
	ext::shared_ptr<PricingEngine> mcBatesEngine =
		MakeMCEuropeanHestonEngine<PseudoRandom>(batesProcess)
		.withStepsPerYear(2)
		.withAntitheticVariate()
		.withAbsoluteTolerance(mcTol)
		.withSeed(1234);

	ext::shared_ptr<PricingEngine> mertonEngine(
		new JumpDiffusionEngine(mertonProcess, 1e-10, 1000));

	for (int i = 1; i <= 5; i += 2) {
		Date exerciseDate = settlementDate + i * Years;
		ext::shared_ptr<Exercise> exercise(
			new EuropeanExercise(exerciseDate));

		VanillaOption batesOption(payoff, exercise);

		batesOption.setPricingEngine(batesEngine);
		Real calculated = batesOption.NPV();

		batesOption.setPricingEngine(mcBatesEngine);
		Real mcCalculated = batesOption.NPV();

		EuropeanOption mertonOption(payoff, exercise);
		mertonOption.setPricingEngine(mertonEngine);
		Real expected = mertonOption.NPV();

		Real tolerance = 2e-8;
		Real relError = std::fabs(calculated - expected) / expected;
		if (relError > tolerance) {
			cout << "failed to reproduce Merton76 price with semi "
				"analytic BatesEngine"
				<< std::fixed << std::setprecision(8)
				<< "\n    calculated: " << calculated
				<< "\n    expected:   " << expected
				<< "\n    rel. error: " << relError
				<< "\n    tolerance:  " << tolerance<<endl;
		}

		Real mcError = std::fabs(expected - mcCalculated);
		if (mcError > 3 * mcTol) {
			cout << "failed to reproduce Merton76 price with Monte-Carlo "
				"BatesEngine"
				<< std::fixed << std::setprecision(8)
				<< "\n    calculated: " << mcCalculated
				<< "\n    expected:   " << expected
				<< "\n    error: " << mcError
				<< "\n    tolerance:  " << mcTol << endl;
		}
	}
}

namespace {
	struct HestonModelData {
		const char* const name;
		Real v0;
		Real kappa;
		Real theta;
		Real sigma;
		Real rho;
		Real r;
		Real q;
	};

	HestonModelData hestonModels[] = {
		// ADI finite difference schemes for option pricing in the 
		// Heston model with correlation, K.J. in t'Hout and S. Foulon,
		{"'t Hout case 1", 0.04, 1.5, 0.04, 0.3, -0.9, 0.025, 0.0},
		// Efficient numerical methods for pricing American options under 
		// stochastic volatility, Samuli Ikonen and Jari Toivanen,
		{"Ikonen-Toivanen", 0.0625, 5, 0.16, 0.9, 0.1, 0.1, 0.0},
		// Not-so-complex logarithms in the Heston model, 
		// Christian Kahl and Peter Jäckel
		{"Kahl-Jaeckel", 0.16, 1.0, 0.16, 2.0, -0.8, 0.0, 0.0},
		// self defined test cases
		{"Equity case", 0.07, 2.0, 0.04, 0.55, -0.8, 0.03, 0.035 },
	};
}

void BatesAnalytic_vs_MonteCarlo() {
	cout << "Testing analytic Bates engine against Monte-Carlo "
		"engine..." << endl;

	SavedSettings backup;

	Date settlementDate(30, March, 2007);
	Settings::instance().evaluationDate() = settlementDate;

	DayCounter dayCounter = ActualActual();
	Date exerciseDate(30, March, 2012);

	ext::shared_ptr<StrikedTypePayoff> payoff(
		new PlainVanillaPayoff(Option::Put, 100));
	ext::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));


	for (Size i = 0; i < sizeof(hestonModels); ++i) {
		Handle<YieldTermStructure> riskFreeTS(flatRate(hestonModels[i].r,
			dayCounter));
		Handle<YieldTermStructure> dividendTS(flatRate(hestonModels[i].q,
			dayCounter));
		Handle<Quote> s0(ext::shared_ptr<Quote>(new SimpleQuote(100)));

		ext::shared_ptr<BatesProcess> batesProcess(new BatesProcess(
			riskFreeTS, dividendTS, s0,
			hestonModels[i].v0,
			hestonModels[i].kappa,
			hestonModels[i].theta,
			hestonModels[i].sigma,
			hestonModels[i].rho, 2.0, -0.2, 0.1));

		const Real mcTolerance = 0.5;
		ext::shared_ptr<PricingEngine> mcEngine =
			MakeMCEuropeanHestonEngine<PseudoRandom>(batesProcess)
			.withStepsPerYear(20)
			.withAntitheticVariate()
			.withAbsoluteTolerance(mcTolerance)
			.withSeed(1234);

		ext::shared_ptr<BatesModel> batesModel(new BatesModel(batesProcess));

		ext::shared_ptr<PricingEngine> fdEngine(
			new FdBatesVanillaEngine(batesModel, 50, 100, 30));

		ext::shared_ptr<PricingEngine> analyticEngine(
			new BatesEngine(batesModel, 160));

		VanillaOption option(payoff, exercise);

		option.setPricingEngine(mcEngine);
		const Real calculated = option.NPV();

		option.setPricingEngine(analyticEngine);
		const Real expected = option.NPV();

		option.setPricingEngine(fdEngine);
		const Real fdCalculated = option.NPV();

		const Real mcError = std::fabs(calculated - expected);
		if (mcError > 3 * mcTolerance) {
			cout << "failed to reproduce Monte-Carlo price for BatesEngine"
				<< "\n    parameter:  " << hestonModels[i].name
				<< std::fixed << std::setprecision(8)
				<< "\n    calculated: " << calculated
				<< "\n    expected:   " << expected
				<< "\n    error: " << mcError
				<< "\n    tolerance:  " << mcTolerance << endl;
		}
		const Real fdTolerance = 0.2;
		const Real fdError = std::fabs(fdCalculated - expected);
		if (fdError > fdTolerance) {
			cout << "failed to reproduce PIDE price for BatesEngine"
				<< "\n    parameter:  " << hestonModels[i].name
				<< std::fixed << std::setprecision(8)
				<< "\n    calculated: " << fdCalculated
				<< "\n    expected:   " << expected
				<< "\n    error: " << fdError
				<< "\n    tolerance:  " << fdTolerance << endl;
		}
	}
}

void BatesDAXCalibration() {
	/* this example is taken from A. Sepp
	   Pricing European-Style Options under Jump Diffusion Processes
	   with Stochstic Volatility: Applications of Fourier Transform
	   http://math.ut.ee/~spartak/papers/stochjumpvols.pdf
	*/

	cout << "Testing Bates model calibration using DAX volatility data..." << endl;

	SavedSettings backup;

	Date settlementDate(5, July, 2002);
	Settings::instance().evaluationDate() = settlementDate;

	DayCounter dayCounter = Actual365Fixed();
	Calendar calendar = TARGET();

	int t[] = { 13, 41, 75, 165, 256, 345, 524, 703 };
	Rate r[] = { 0.0357,0.0349,0.0341,0.0355,0.0359,0.0368,0.0386,0.0401 };

	std::vector<Date> dates;
	std::vector<Rate> rates;
	dates.push_back(settlementDate);
	rates.push_back(0.0357);
	for (Size i = 0; i < 8; ++i) {
		dates.push_back(settlementDate + t[i]);
		rates.push_back(r[i]);
	}
	// FLOATING_POINT_EXCEPTION
	Handle<YieldTermStructure> riskFreeTS(
		ext::shared_ptr<YieldTermStructure>(
			new ZeroCurve(dates, rates, dayCounter)));

	Handle<YieldTermStructure> dividendTS(
		flatRate(settlementDate, 0.0, dayCounter));

	Volatility v[] =
	{ 0.6625,0.4875,0.4204,0.3667,0.3431,0.3267,0.3121,0.3121,
	  0.6007,0.4543,0.3967,0.3511,0.3279,0.3154,0.2984,0.2921,
	  0.5084,0.4221,0.3718,0.3327,0.3155,0.3027,0.2919,0.2889,
	  0.4541,0.3869,0.3492,0.3149,0.2963,0.2926,0.2819,0.2800,
	  0.4060,0.3607,0.3330,0.2999,0.2887,0.2811,0.2751,0.2775,
	  0.3726,0.3396,0.3108,0.2781,0.2788,0.2722,0.2661,0.2686,
	  0.3550,0.3277,0.3012,0.2781,0.2781,0.2661,0.2661,0.2681,
	  0.3428,0.3209,0.2958,0.2740,0.2688,0.2627,0.2580,0.2620,
	  0.3302,0.3062,0.2799,0.2631,0.2573,0.2533,0.2504,0.2544,
	  0.3343,0.2959,0.2705,0.2540,0.2504,0.2464,0.2448,0.2462,
	  0.3460,0.2845,0.2624,0.2463,0.2425,0.2385,0.2373,0.2422,
	  0.3857,0.2860,0.2578,0.2399,0.2357,0.2327,0.2312,0.2351,
	  0.3976,0.2860,0.2607,0.2356,0.2297,0.2268,0.2241,0.2320 };

	Handle<Quote> s0(ext::shared_ptr<Quote>(new SimpleQuote(4468.17)));
	Real strike[] = { 3400,3600,3800,4000,4200,4400,
					  4500,4600,4800,5000,5200,5400,5600 };


	Real v0 = 0.0433;
	ext::shared_ptr<SimpleQuote> vol(new SimpleQuote(std::sqrt(v0)));
	ext::shared_ptr<BlackVolTermStructure> volTS =
		flatVol(settlementDate, vol, dayCounter);

	const Real kappa = 1.0;
	const Real theta = v0;
	const Real sigma = 1.0;
	const Real rho = 0.0;
	const Real lambda = 1.1098;
	const Real nu = -0.1285;
	const Real delta = 0.1702;

	ext::shared_ptr<BatesProcess> process(
		new BatesProcess(riskFreeTS, dividendTS, s0, v0,
			kappa, theta, sigma, rho, lambda, nu, delta));

	ext::shared_ptr<BatesModel> batesModel(new BatesModel(process));

	ext::shared_ptr<PricingEngine> batesEngine(
		new BatesEngine(batesModel, 64));

	std::vector<ext::shared_ptr<BlackCalibrationHelper> > options;

	for (Size s = 0; s < 13; ++s) {
		for (Size m = 0; m < 8; ++m) {
			Handle<Quote> vol(ext::shared_ptr<Quote>(
				new SimpleQuote(v[s * 8 + m])));

			Period maturity((int)((t[m] + 3) / 7.), Weeks); // round to weeks

			// this is the calibration helper for the bates models
			// FLOATING_POINT_EXCEPTION
			options.push_back(ext::shared_ptr<BlackCalibrationHelper>(
				new HestonModelHelper(maturity, calendar,
					s0->value(), strike[s], vol,
					riskFreeTS, dividendTS,
					BlackCalibrationHelper::ImpliedVolError)));
			options.back()->setPricingEngine(batesEngine);
		}
	}

	// check calibration engine
	LevenbergMarquardt om;
	batesModel->calibrate(options, om,
		EndCriteria(400, 40, 1.0e-8, 1.0e-8, 1.0e-8));

	Real expected = 36.6;
	Real calculated = getCalibrationError(options);
	if (std::fabs(calculated - expected) > 2.5)
		cout << "failed to calibrate the bates model"
			<< "\n    calculated: " << calculated
			<< "\n    expected:   " << expected << endl;

	//check pricing of derived engines
	std::vector<ext::shared_ptr<PricingEngine> > pricingEngines;

	process = ext::shared_ptr<BatesProcess>(
		new BatesProcess(riskFreeTS, dividendTS, s0, v0,
			kappa, theta, sigma, rho, 1.0, -0.1, 0.1));

	pricingEngines.push_back(ext::shared_ptr<PricingEngine>(
		new BatesDetJumpEngine(
			ext::make_shared<BatesDetJumpModel>(
				process), 64)));

	ext::shared_ptr<HestonProcess> hestonProcess(new HestonProcess(
		riskFreeTS, dividendTS, s0, v0, kappa, theta, sigma, rho));

	pricingEngines.push_back(ext::shared_ptr<PricingEngine>(
		new BatesDoubleExpEngine(
			ext::make_shared<BatesDoubleExpModel>(
				hestonProcess, 1.0), 64)));

	pricingEngines.push_back(ext::shared_ptr<PricingEngine>(
		new BatesDoubleExpDetJumpEngine(
			ext::make_shared<BatesDoubleExpDetJumpModel>(
				hestonProcess, 1.0), 64)));

	Real expectedValues[] = { 5896.37,
							  5499.29,
							  6497.89 };

	Real tolerance = 0.1;
	for (Size i = 0; i < pricingEngines.size(); ++i) {
		for (Size j = 0; j < options.size(); ++j) {
			options[j]->setPricingEngine(pricingEngines[i]);
		}

		Real calculated = std::fabs(getCalibrationError(options));
		if (std::fabs(calculated - expectedValues[i]) > tolerance)
			cout << "failed to calculated prices for derived Bates models"
				<< "\n    calculated: " << calculated
				<< "\n    expected:   " << expectedValues[i] << endl;
	}
}

/*test_suite* BatesModelTest::suite() {
	test_suite* suite = BOOST_TEST_SUITE("Bates model tests");
	suite->add(QUANTLIB_TEST_CASE(&BatesModelTest::testAnalyticVsBlack));
	suite->add(QUANTLIB_TEST_CASE(
		&BatesModelTest::testAnalyticAndMcVsJumpDiffusion));
	suite->add(QUANTLIB_TEST_CASE(&BatesModelTest::testAnalyticVsMCPricing));
	// FLOATING_POINT_EXCEPTION
	suite->add(QUANTLIB_TEST_CASE(&BatesModelTest::testDAXCalibration));
	return suite;
}
*/

// main program (working on Kou Double-Exponential Jump Diffusion Process
int main(int, char*[]) {
	try {

		boost::timer timer;
		cout << endl;

		Calendar calendar = TARGET();
		Date todaysDate(10, Jan, 2019);
		Date settlementDate(11, Jan, 2019);
		Settings::instance().evaluationDate() = todaysDate;

		// our options
		Option::Type type(Option::Put);
		Real underlying = 36;
		Real strike = 40;
		Spread dividendYield = 0.00;
		Rate riskFreeRate = 0.06;
		Volatility volatility = 0.20;
		Date maturity(11, Jan, 2020);
		DayCounter dayCounter = Actual365Fixed();

		cout << "Option type = " << type << endl;
		cout << "Maturity = " << maturity << endl;
		cout << "Underlying price = " << underlying << endl;
		cout << "Strike = " << strike << endl;
		cout << "Risk-free interest rate = " << QuantLib::io::rate(riskFreeRate)
			<< endl;
		cout << "Dividend yield = " << QuantLib::io::rate(dividendYield)
			<< endl;
		cout << "Volatility = " << QuantLib::io::volatility(volatility)
			<< endl;
		cout << endl;
		string method;
		cout << endl;

		// write column headings
		Size widths[] = { 35, 14 };
		cout << setw(widths[0]) << left << "Method"
			<< setw(widths[1]) << left << "European"
			<< endl;

		vector<Date> exerciseDates;
		for (int i = 1; i <= 2; i++)
			exerciseDates.push_back(settlementDate + 3 * i*Months);

		ext::shared_ptr<Exercise> europeanExercise(
			new EuropeanExercise(maturity));

		Handle<Quote> underlyingH(
			ext::shared_ptr<Quote>(new SimpleQuote(underlying)));

		// bootstrap the yield/dividend/vol curves
		Handle<YieldTermStructure> flatTermStructure(
			ext::shared_ptr<YieldTermStructure>(
				new FlatForward(settlementDate, riskFreeRate, dayCounter)));
		Handle<YieldTermStructure> flatDividendTS(
			ext::shared_ptr<YieldTermStructure>(
				new FlatForward(settlementDate, dividendYield, dayCounter)));
		Handle<BlackVolTermStructure> flatVolTS(
			ext::shared_ptr<BlackVolTermStructure>(
				new BlackConstantVol(settlementDate, calendar, volatility,
					dayCounter)));
		ext::shared_ptr<StrikedTypePayoff> payoff(
			new PlainVanillaPayoff(type, strike));
		ext::shared_ptr<BlackScholesMertonProcess> bsmProcess(
			new BlackScholesMertonProcess(underlyingH, flatDividendTS,
				flatTermStructure, flatVolTS));

		// options
		VanillaOption europeanOption(payoff, europeanExercise);
		
		// Black-Scholes for European
		method = "Black-Scholes";
		europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
			new AnalyticEuropeanEngine(bsmProcess)));
		cout << setw(widths[0]) << left << method
			<< fixed
			<< setw(widths[1]) << left << europeanOption.NPV()
			<< endl;

		// semi-analytic Heston for European
		method = "Heston semi-analytic";
		ext::shared_ptr<HestonProcess> hestonProcess(
			new HestonProcess(flatTermStructure, flatDividendTS,
				underlyingH, volatility*volatility,
				1.0, volatility*volatility, 0.001, 0.0));
		ext::shared_ptr<HestonModel> hestonModel(
			new HestonModel(hestonProcess));
		europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(
			new AnalyticHestonEngine(hestonModel)));
		cout << setw(widths[0]) << left << method
			<< fixed
			<< setw(widths[1]) << left << europeanOption.NPV()
			<< endl;



		/* THIS IS THE MODEL I AM WORKING ON. IT USES BATES PROCESS AND BATES DOUBLE EXP ENGINE
		IT IS A MEMBER OF HESTON PROCESS, MAYBE THAT IS WHY I HAVE AN ERROR ON BatesProcess(flatTermStructure,....
		please help! */
		method = "Kou Double-Jump Diffusion";
		ext::shared_ptr<BatesProcess> batesProcess(
			new BatesProcess(flatTermStructure, flatDividendTS, underlyingH, volatility*volatility,
				1.0, volatility*volatility, 0.001, 0.0, 1e-14, 1e-14, 1e-14));
		ext::shared_ptr<BatesDoubleExpModel> batesDoubleExpModel(new BatesDoubleExpModel(batesProcess));
		europeanOption.setPricingEngine(ext::shared_ptr<PricingEngine>(new BatesDoubleExpEngine(batesDoubleExpModel)));
		cout << setw(widths[0]) << left << method << fixed << setw(widths[1]) << left << europeanOption.NPV()
			<< endl;

		


		// random QuantLib Timer & Debug Functions below!
		// just ignore them 
		// or comment them out for now
		double seconds = timer.elapsed();
		int hours = int(seconds / 3600);
		seconds -= hours * 3600;
		int minutes = int(seconds / 60);
		seconds -= minutes * 60;
		cout << " \nRun completed in ";
		if (hours > 0)
			cout << hours << " h ";
		if (hours > 0 || minutes > 0)
			cout << minutes << " m ";
		cout << fixed << setprecision(0)
			<< seconds << " s\n" << endl;
		return 0;

		} catch (std::exception& e) {
			cerr << e.what() << endl;
			return 1;
		} catch (...) {
			cerr << "unknown error" << endl;
			return 1;
		}

	system("pause");

};