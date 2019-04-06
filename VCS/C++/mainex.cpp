// QlTesting.cpp : Defines the entry point for the console application.
//#include "stdafx.h"

#include <ql/quantlib.hpp>
using namespace QuantLib;


// Author: Dimitri Reiswich



// #include "TestingMacros1.h"


//#include "Matrix3.h"

// #include "DateTesting1.h"
// #include "DateTesting2.h"
// #include "DateTesting3.h"
// #include "DateTesting4.h"

 // #include "Singleton2.h"
// #include "Singleton3.h"
// #include "LazyObject3.h"



#include "CalendarTesting1.h"
#include "CalendarTesting2.h"
#include "CalendarTesting3.h"

#include "DayCounter1.h"

//#include "Schedule2.h"
//#include "Schedule3.h"
//#include "Schedule4.h"


//#include "Integration1.h"
//#include "Integration2.h"
//#include "Integration3.h"
//#include "Integration4.h"
/*
#include "Solver2.h"
#include "Solver3.h"
#include "Solver4.h"
#include "Solver5.h"

#include "Copulas1.h"

#include "Interpolations1.h"
#include "Interpolations2.h"
#include "Interpolations4.h"
#include "Interpolations6.h"

#include "Matrix1.h"
#include "Matrix2.h"

#include "Optimizer1.h"
#include "Optimizer3.h"
#include "Optimizer4.h"
#include "Optimizer5.h"

#include "RandomNumbers1.h"
#include "RandomNumbers2.h"
#include "RandomNumbers3.h"
#include "RandomNumbers4.h"
#include "RandomNumbers5.h"

#include "Null1.h"
#include "Design2.h"
#include "Design3.h"
#include "Design4.h"
#include "Design4a.h"
#include "Design4b.h"
#include "Design5.h"
#include "Design6.h"
#include "Design7.h"
//#include "Design8.h"
//#include "Design10.h"
//#include "Design9.h"
//#include "Design12.h"

#include "Factory1.h"
#include "Factory2.h"
#include "Factory3.h"
#include "Factory4.h"

#include "Handle1.h"
#include "Handle2.h"*/

//#include "YieldCurve2.h"
//#include "YieldCurve4.h"
//#include "YieldCurve5.h"
// #include "YieldCurve6.h"
#include "YieldCurve9.h"
#include "YieldCurve10.h"

/*
#include "Indexes1.h"

#include "GetYieldCurve.h"
#include "GetBondYieldCurve.h"
#include "GetZeroYieldCurve.h"

#include "Bonds1.h"
#include "Bonds2.h"
#include "Bonds3.h"

#include "Coupon2.h"
#include "Coupon3.h"
#include "Coupon4.h"

#include "Swaps1.h"
#include "Swaps2.h"

#include "Test1.h"
#include "Test2.h"

#include "Visitor1.h"
#include "Visitor2.h"
#include "Visitor3.h"

#include "BSPricingEngines1.h"
#include "BlackScholesCalculator.h"

#include "VolatilityObjects1.h"
#include "VolatilityObjects2.h"

#include "StochasticProcesses1.h"
#include "StochasticProcesses2.h"
#include "StochasticProcesses3.h"
#include "StochasticProcesses4.h"*/



// #include "Integration4.h"

// #include "DateTesting1.h"


// #include "YieldCurve9.h"

int _tmain(int __argc __argv[])
{
	try{

//		testingSingleton1();
// testIntegration4();


		// DateTesting1();

// testingLazyObject1();

//	testingMacros1();

		
//	DateTesting1();
//	DateTesting2();
//	DateTesting3();
//	DateTesting4();


//  	testingSingleton1();
//		testingSingleton2();

//		testingMatrix3();

//		testingStochasticProcesses1();
//		testingStochasticProcesses2();
//		testingStochasticProcesses3();

//		testingStochasticProcesses4();

	//	testingVolatilityObjects1();
	//	testingVolatilityObjects2();


	//	testingBsPricingEngines1();
	//	testingBlackScholesCalculator();

	//	testingVisitor1();

	//	testingNull();

	//	testingBonds1();
	//	testingBonds2();
	//	testingBonds3();

	//	testingSwaps1();
	//	testingSwaps2();


	//	testingCoupons1();
	//	testingCoupons2();
	//	testingCoupons3();

	//	testingFactory1();
	//	testingFactory2();

	//	testingNull1();

	//	testingYields1();
	//	testingYields2();
	//	testingYields3();
		testingYields4();
	//	testingYields5();

	//	testingIndexes1();
	//	testHandle1();

	//	testingDesignPatterns1();
	//	testingDesignPatterns2();
	//	testingDesignPatterns2a();
	//	testingDesignPatterns3();
	//	testingDesignPatterns4();


//	testingRandomNumbers1();
//	testingRandomNumbers2();
//	testingRandomNumbers3();
//	testingRandomNumbers4();
//	testingRandomNumbers5();

//	testOptimizer2();

//	testingInterpolations1();
//	testingInterpolations2();
//	testingInterpolations3();
//	testingInterpolations4();

//	testingMatrix1();
//	testingMatrix2();

//	testCopulas1();

//	testIntegration1();
//	testIntegration2();
//	testIntegration3();
//	testIntegration4();

//	testSolver1();
//	testSolver2();

//	CalendarTesting1();
//	CalendarTesting2();
//	CalendarTesting3();

//	dayCounterTesting1();
//	testingMacros1();
//	testingMacros2();

//	testingSchedule1();
//	testingSchedule3();
//	testingSchedule4();


	} catch (std::exception& e) {
        std::cout << e.what() << std::endl;
		    return 1;
    } catch (...) {
        std::cout << "unknown error" << std::endl;
	       return 1;
    }


	return 0;
}




