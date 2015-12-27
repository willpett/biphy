/**
 * @file
 * This file contains the implementation of Exception, which
 * is used to handle eceptions in RevBayes.
 *
 * @brief Implementation of Exception
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */

#include "Exception.h"

#include <string>
#include <iostream>


/** Static string with names of exception types for printing */
std::string Exception::exceptionName[] = { "Default", "Quit", "Missing Variable" };


/** Default constructor */
Exception::Exception(void) : exceptionType(DEFAULT), message() {
}


/** Message constructor */
Exception::Exception(const std::string& msg) :exceptionType(DEFAULT), message(msg) {
}


/** Message constructor from stringstream */
Exception::Exception(const std::ostringstream& msg) : exceptionType(DEFAULT), message(msg.str()) {
}


/** General constructor */
Exception::Exception(exceptionT type, const std::string& msg) : exceptionType(type), message(msg) {
}


std::string Exception::getMessage(void) const {

    return message;
}


void Exception::print(std::ostream &o) const {
    
    std::string errorType;
    switch (exceptionType) {
        case DEFAULT:
            errorType = "Error";
            break;
        case QUIT:
            errorType = "Quit";
            break;
        case MISSING_VARIABLE:
            errorType = "Missing Variable";
            break;
            
        default:
            errorType = "Error";
    }
    o << errorType << ":\t" << message;
}

void Exception::setMessage(std::string msg) {

    message = msg;
}


