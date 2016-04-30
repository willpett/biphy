//
//  MathHelper.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 11/27/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "MathHelper.h"


double Math::Helper::fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}

double Math::Helper::fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}
