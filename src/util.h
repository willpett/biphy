/*
 * util.h
 *
 *  Created on: Oct 29, 2014
 *      Author: walker
 */

#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_

using namespace std;

bool fexists(const std::string& filename) {
  ifstream ifile(filename.c_str());
  return ifile;
}

const char digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

int IsInt(string s)	{
	int returnValue = 1;
	unsigned int i = 0;
	if ((s[0] == '+') || (s[0] == '-')) i++;
	if (i == s.length())	returnValue = 0;

	while ( returnValue && (i < s.length()) )	{
		int j = 0;
		while ((j<10) &&  (digit[j] != s[i]))	{
			j++;
		}
		if (j == 10)	{
			returnValue = 0;
		}
		i++;
	}
	return returnValue;
}

int IsFloat(string s)	{
	int returnValue = 1;
	unsigned int i = 0;

	while ( returnValue && (i < s.length()) )	{
		int j = 0;
		while ( (j<10) && (digit[j] != s[i]))	{
			j++;
		}
		if (j == 10)	{
			if ( ! ( (s[i] == '.') || (s[i] == '-') || (s[i] == '+') || (s[i] == 'e') ) )	{
				returnValue = 0;
			}
		}
		i++;
	}
	return returnValue;
}

void our_terminate (void);

namespace {
    static const bool SET_TERMINATE = set_terminate(our_terminate);
}

void our_terminate (void) { // try 1
	static bool tried_throw = false;
    try {
    	if(!tried_throw++) throw;
    }
    catch (RbException& e) {
    	cerr << "RbException:\t" << e.getMessage() << '\n';
    }
}


#endif /* SRC_UTIL_H_ */
