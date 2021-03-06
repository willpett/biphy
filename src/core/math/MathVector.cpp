/**
 * @file MathMatrix
 * This file contains the math on matrices.
 *
 * @brief Implementation of matrix algebra.
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


#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "MathVector.h"

#include "Exception.h"
#include "Settings.h"

// Vector Functions


/*!
 * This function normalizes a vector such that its sum is some value.
 *
 * \brief Vector normalization function.
 * \param x is a reference to the vector to be normalized .
 * \param sum is the desired sum for the normalized vector.
 * \return Does not return a value. 
 * \throws Does not throw an error.
 */
void Math::normalize(std::vector<double>& x, double sum) {
    
    double s = 0.0;
    for (size_t i=0; i<x.size(); i++)
        s += x[i];
    double factor = sum / s;
    for (size_t i=0; i<x.size(); i++)
        x[i] *= factor;
}




