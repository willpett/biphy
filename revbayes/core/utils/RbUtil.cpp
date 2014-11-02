//
//  RbUtil.cpp
//  RevBayesCore
//
//  Created by Sebastian Hoehna on 7/6/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "RbException.h"
#include "RbUtil.h"
#include "Clade.h"


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<bool>& x) {
    o << "( ";
    for (std::vector<bool>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        if ( *it ) {
            o << "TRUE";
        } else {
            o << "FALSE";
        }
    }
    o << " )";
    
    return o;
}

std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<double>& x) {
    o << "( ";
    for (std::vector<double>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<int>& x) {
    o << "( ";
    for (std::vector<int>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<unsigned int>& x) {
    o << "( ";
    for (std::vector<unsigned int>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<std::string>& x) {
    o << "( ";
    for (std::vector<std::string>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<Clade>& x) {
    o << "( ";
    for (std::vector<Clade>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<Taxon>& x) {
    o << "( ";
    for (std::vector<Taxon>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}


std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<TimeTree>& x) {
    o << "( ";
    for (std::vector<TimeTree>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}

/*
std::ostream& RevBayesCore::operator<<(std::ostream& o, const std::vector<Trace>& x) {
    o << "( ";
    for (std::vector<Trace>::const_iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            o << " , ";
        }
        o << *it;
    }
    o << " )";
    
    return o;
}
*/

std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<bool>& x) {
	std::string tmp;
	is >> tmp;
    for (std::vector<bool>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> tmp;
        if ( tmp == "TRUE" ) {
            *it = true;
        } else {
        	*it = false;
        }
    }
    is >> tmp;

    return is;
}

std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<double>& x) {
    for (std::vector<double>::iterator it = x.begin(); it != x.end(); ++it)
    	is >> *it;

    return is;
}


std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<int>& x) {
	std::string tmp;
	is >> tmp;
    for (std::vector<int>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> *it;
    }
    is >> tmp;

    return is;
}


std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<unsigned int>& x) {
	std::string tmp;
	is >> tmp;
    for (std::vector<unsigned int>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> *it;
    }
    is >> tmp;

    return is;
}


std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<std::string>& x) {
	std::string tmp;
	is >> tmp;
    for (std::vector<std::string>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> *it;
    }
    is >> tmp;

    return is;
}


std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<TimeTree>& x) {
	std::string tmp;
	is >> tmp;
    for (std::vector<TimeTree>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> *it;
    }
    is >> tmp;

    return is;
}

/*
std::istream& RevBayesCore::operator>>(std::istream& is, std::vector<Trace>& x) {
    std::string tmp;
	is >> tmp;
    for (std::vector<Trace>::iterator it = x.begin(); it != x.end(); ++it) {
        if ( it != x.begin() ) {
            is >> tmp;
        }
        is >> *it;
    }
    is >> tmp;

    return is;
}
*/

std::vector<int> RevBayesCore::operator+(const std::vector<int>& x, const std::vector<int>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only add vectors of same size!");
    }
    
    std::vector<int> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] + y[i];
    }
    
    return z;
}


std::vector<double> RevBayesCore::operator+(const std::vector<double>& x, const std::vector<double>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only add vectors of same size!");
    }
    
    std::vector<double> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] + y[i];
    }
    
    return z;
}


std::vector<int> RevBayesCore::operator-(const std::vector<int>& x, const std::vector<int>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only subtract vectors of same size!");
    }
    
    std::vector<int> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] - y[i];
    }
    
    return z;
}


std::vector<double> RevBayesCore::operator-(const std::vector<double>& x, const std::vector<double>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only subtract vectors of same size!");
    }
    
    std::vector<double> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] - y[i];
    }
    
    return z;
}


std::vector<int> RevBayesCore::operator*(const std::vector<int>& x, const std::vector<int>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only multiply vectors of same size!");
    }
    
    std::vector<int> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] * y[i];
    }
    
    return z;
}


std::vector<double> RevBayesCore::operator*(const std::vector<double>& x, const std::vector<double>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only multiply vectors of same size!");
    }
    
    std::vector<double> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] * y[i];
    }
    
    return z;
}


std::vector<double> RevBayesCore::operator/(const std::vector<int>& x, const std::vector<int>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only divide vectors of same size!");
    }
    
    std::vector<double> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] / double(y[i]);
    }
    
    return z;
}


std::vector<double> RevBayesCore::operator/(const std::vector<double>& x, const std::vector<double>& y) 
{
    
    size_t n = x.size();
    
    if ( n != y.size() )
    {
        throw RbException("Can only divide vectors of same size!");
    }
    
    std::vector<double> z(n,0);
    
    for (size_t i = 0; i < n; ++i) 
    {
        z[i] = x[i] / y[i];
    }
    
    return z;
}
