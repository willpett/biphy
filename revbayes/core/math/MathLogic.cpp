/**
 * @file MathLogic
 * This file contains some logic function on numbers.
 *
 * The following functions check whether two double-precision real
 * numbers are equal. They are described in:
 *
 * Knuth, D. E. 1981. The art of computer programming: Seminumerical
 *    algorithms, Volume 2. Addison-Wesley.
 *
 * Note that approximately equal to is more stringent than essentially
 * equal to.
 *
 * @brief Implementation of simple RlBoolean algebra on numbers.
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes core development team
 * @license GPL version 3
 * @version 1.0
 * @since 2011-03-17, version 1.0
 *
 * $Id$
 */

#include "MathLogic.h"
#include <cassert>
#include <cmath>
#include "Constants.h"
#include "Exception.h"
#include "Settings.h"

/** Compares two doubles for equality */
bool Math::compApproximatelyEqual(double a, double b) {
    
    double epsilon = Settings::userSettings().getTolerance();
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Compares two doubles for equality */
bool Math::compApproximatelyEqual(double a, double b, double epsilon) {
    
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Compares two doubles for equality */
bool Math::compEssentiallyEqual(double a, double b) {
    
    double epsilon = Settings::userSettings().getTolerance();
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Compares two doubles for equality */
bool Math::compEssentiallyEqual(double a, double b, double epsilon) {
    
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Tests whether one number is greater than another */
bool Math::compDefinitelyGreaterThan(double a, double b) {
    
    double epsilon = Settings::userSettings().getTolerance();
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Tests whether one number is greater than another */
bool Math::compDefinitelyGreaterThan(double a, double b, double epsilon) {

    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Tests whether one number is less than another */
bool Math::compDefinitelyLessThan(double a, double b) {
    
    double epsilon = Settings::userSettings().getTolerance();
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


/** Tests whether one number is less than another */
bool Math::compDefinitelyLessThan(double a, double b, double epsilon) {
	//return fabs(b - a) >= epsilon;
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}


bool Math::isAComputableNumber(double x)
{
    if ( x != x )
        return false;
    if ( x == Constants::Double::neginf )
        return false;
    if ( x == Constants::Double::nan )
        return false;
    if ( x == Constants::Double::inf )
        return false;
    
    return true;
}

/** Tests whether a double is finite */
bool Math::isFinite(double x) {
    
    return x == x;
}


/** Tests whether a double is actually an interger */
bool Math::isInt(double x){

    double int_part;
    return ( modf ( x, &int_part ) == 0.0 );
}


/** Tests whether a double is NAN (not a number) */
bool Math::isNan(double x) {

    return x != x;
}


/** Returns the maximum of two real numbers */
double Math::max(double x, double y){

	return (x < y) ? y : x;
}


/** Returns the minimum of two real numbers */
double Math::min(double x, double y) {

    return (x < y) ? x : y;
}

double Math::sign(double a) {

	return a < 0 ? -1 : (a == 0 ? 0 : 1);
}

double Math::sign(double a, double b	) {

	return fabs(a)*sign(b);
}


