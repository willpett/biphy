/**
 * @file
 * This file helper functions for manipulating RlStrings in RevBayes.
 *
 * @brief Namespace containing helper functions for RlStrings
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date$
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 *
 * $Id$
 */
#ifndef RlStringUtilities_H
#define RlStringUtilities_H

#include <string>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

namespace StringUtilities {
    
    std::string                 formatStringForScreen(  const std::string   s,
                                                        std::string         firstLinePad,
                                                        std::string         hangingPad,
                                                        size_t              screenWidth);                           //!< Format RlString for output to screen
    std::string                 getFileContentsAsString(std::string& s);                                            //!< Convert the file contents to a RlString
    std::string                 getLastPathComponent(const std::string& s);                                         //!< Find the last component of a file path
    void                        stringSplit(std::string str, std::string delim, std::vector<std::string>& results); //!< Split a RlString into pieces
    void                        toLower(std::string& str);                                                          //!< Convert RlString's characters to lower case
    bool                        isNumber(std::string& s);                                                           //!< Checks if a RlString is a number
}

const char digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    
inline int Int(std::string s)   {
    return std::atoi( s.c_str() );
}

inline double Double(std::string s) {
    return std::atof( s.c_str() );
}

inline int GoPastNext (std::istream& is, const char inChar)  {
    unsigned char c;
    if (! is.eof()) {
        do  {
            is >> c;
        }
        while( c != inChar && ! is.eof());
    }
    return ! is.eof();
}

inline std::string ReadLine(std::istream& is) {
    std::string str = "";
    char c;
    do  {
        is.get(c);
        if (c != '\n')  {
            str += c;
        }
    }   while (c != '\n' && ! is.eof());
    return str;
}

inline void GoPastNextWord(std::istream& is, const std::string inWord)    {

    unsigned int k = 0;
    char c;
    while ((!is.eof()) && (k<inWord.length()))  {
        is.get(c);
        if ((c >=65) && (c <= 90))  {
            c += 32;
        }
        char ca = inWord[k];
        if ((ca >=65) && (ca <= 90))    {
            ca += 32;
        }
        if (c == ca)    {
            k++;
        }
        else    {
            k=0;
        }
    }
        
}

inline int EquivalentStrings(std::string a, std::string b)    {

    if (a.length() != b.length())   {
        return 0;
    }
    unsigned int k = 0;
    int cont = 1;
    while ((k < a.length()) && (cont))  {
        char ca = a[k];
        char cb = b[k];
        if ((ca >=65) && (ca <= 90))    {
            ca += 32;
        }
        if ((cb >=65) && (cb <= 90))    {
            cb += 32;
        }
        if (ca != cb)   {
            cont = 0;
        }
        k++;
    }
    return cont;    
}

inline void GoPastNextLine(std::istream& is, const std::string inLine)    {
    std::string theLine;
    do  {
        theLine = ReadLine(is);
        std::cerr << theLine << "\n";
    }
    while(! EquivalentStrings(theLine,inLine));
}



inline std::string StringReplace(char c, std::string by, std::string s)    {
    std::string tmp;
    for (unsigned int i=0; i<s.length(); i++)   {
        if (s[i] == c)  {
            tmp += by;
        }
        else    {
            tmp += s[i];
        }
    }
    return tmp;
}   
    
inline int EmptyLine(std::string s)  {

    int unsigned n = 0;
    while ( (n<s.length()) && ((s[n] == ' ') || (s[n] == '\t') || (s[n] == '\n')) ) {
        n++;
    }
    return (n == s.length()) ;
}

inline std::string Filter(std::string input, char c)  {

    std::string temp = "";
    for (int unsigned i=0; i<input.length(); i++)   {
        if (input[i] != c)  {
            temp += input[i];
        }
    }
    return temp;
}

inline int IsInt(std::string s)  {
    int returnValue = 1;
    unsigned int i = 0;
    if ((s[0] == '+') || (s[0] == '-')) i++;
    if (i == s.length())    returnValue = 0;

    while ( returnValue && (i < s.length()) )   {
        int j = 0;
        while ( (j<10) && (digit[j] != s[i]))   {
            j++;
        }
        if (j == 10)    {
            returnValue = 0;
        }
        i++;
    }
    return returnValue;
}

inline int IsFloat(std::string s)    {
    int returnValue = 1;
    unsigned int i = 0;
    if ((s[0] == '+') || (s[0] == '-')) i++;
    if (i == s.length())    returnValue = 0;
    
    while ( returnValue && (i < s.length()) )   {
        int j = 0;
        while ((j<10) && (digit[j] != s[i]))    {
            j++;
        }
        if (j == 10)    {
            if ( ! ( (s[i] == '.') || (s[i] == '-') || (s[i] == '+') || (s[i] == 'e') ) )   {
                returnValue = 0;
            }
        }
        i++;
    }
    return returnValue;
}

inline int IsDigit(char c)  {
    int returnValue = 0;
    int i=0;
    while (! returnValue && i<10)   {
        returnValue = (c == digit[i]);
        i++;
    }
    return returnValue;
}

inline double Decimal(double d, int ndigit) {

    double precision = 1;
    for (int k=0; k<ndigit; k++)    {
        precision *= 10;
    }
    return  ((double) ((int) (precision * d + 0.1/precision))) / precision;
}

std::string                     operator+(const std::string& A, double B);
std::string                     operator+(const std::string& A, int B);
std::string                     operator+(const std::string& A, size_t B);
std::string                     operator+(double A, const std::string& B);
std::string                     operator+(int A, const std::string& B);
std::string                     operator+(size_t A, const std::string& B);


#endif
