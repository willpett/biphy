/**
 * @file
 * This file contains the declaration of Matrix,
 * a container type used to hold value matrices for the inference machinery.
 *
 * @brief Declaration of Matrix
 *
 * (c) Copyright 2009- under GPL version 3
 * @date Last modified: $Date: 2012-03-10 12:55:11 +0100 (Sat, 10 Mar 2012) $
 * @author The RevBayes Development Core Team
 * @license GPL version 3
 * @version 1.0
 * @since 2012-05-08, version 1.0
 *
 * $Id: Matrix.h 1327 2012-03-10 11:55:11Z hoehna $
 */

#ifndef MatrixReal_H
#define MatrixReal_H

#include <cstddef>
#include <iostream>
#include <vector>

namespace RevBayesCore {
    
    class MatrixReal {
        
    public:
        MatrixReal(size_t n, size_t k);
        MatrixReal(size_t n, size_t k, double v);
        MatrixReal(const MatrixReal& m);
        virtual                                ~MatrixReal(void);
        
        
        // overloaded operators
        MatrixReal&                             operator=(const MatrixReal& m);
        std::vector<double>&                    operator[](size_t index);
        const std::vector<double>&              operator[](size_t index) const;
        
        // global operators
        MatrixReal&                             operator+=(double b);                                               //!< operator += for scalar 
        MatrixReal&                             operator-=(double b);                                               //!< operator -= for scalar 
        MatrixReal&                             operator*=(double b);                                               //!< operator *= for scalar 
        MatrixReal&                             operator+=(const MatrixReal& B);                                    //!< operator += 
        MatrixReal&                             operator-=(const MatrixReal& B);                                    //!< operator -= 
        MatrixReal&                             operator*=(const MatrixReal& B);                                    //!< operator *= (matrix multiplication)
        MatrixReal                              operator+(double b) const;                                          //!< operator + for matrix + scalar 
        MatrixReal                              operator-(double b) const;                                          //!< operator - for scalar 
        MatrixReal                              operator*(double b) const;                                          //!< operator * for scalar 
        MatrixReal                              operator+(const MatrixReal& B) const;                               //!< operator + 
        MatrixReal                              operator-(const MatrixReal& B) const;                               //!< operator - 
        MatrixReal                              operator*(const MatrixReal& B) const;                               //!< operator * (matrix multiplication) 
        std::vector<double>                     operator*(const std::vector<double> &b) const;                                          //!< operator * for scalar 
        

        std::vector<std::vector<double> >::const_iterator       begin(void) const;
        std::vector<std::vector<double> >::iterator             begin(void);
        std::vector<std::vector<double> >::const_iterator       end(void) const;
        std::vector<std::vector<double> >::iterator             end(void);
        
        // utility funcions
        void                                    clear(void);
        size_t                                  getNumberOfColumns(void) const;
        size_t                                  getNumberOfRows(void) const;
        size_t                                  size(void) const;
        void                                    resize(size_t r, size_t c);
        
    protected:
        std::vector<std::vector<double> >       elements;
        size_t                                  nRows;
        size_t                                  nCols;
    };
    
    // Global functions using the class
    std::ostream&                       operator<<(std::ostream& o, const MatrixReal& x);                                           //!< Overloaded output operator

    
//    MatrixReal                            operator+(const MatrixReal& A);                                             //!< Unary operator + 
//    MatrixReal                            operator-(const MatrixReal& A);                                             //!< Unary operator - 
//    MatrixReal                            operator/(const MatrixReal& A, const MatrixReal& B);                        //!< operator / for matrix / matrix 
//    MatrixReal                            operator+(double a, const MatrixReal& B);                            //!< operator + for scalar + matrix 
//    MatrixReal                            operator-(double a, const MatrixReal& B);                            //!< operator - for scalar - matrix 
//    MatrixReal                            operator*(double a, const MatrixReal& B);                            //!< operator * for scalar * matrix 
//    MatrixReal                            operator/(double a, const MatrixReal& B);                            //!< operator / for scalar / matrix 
//    MatrixReal                            operator*(const MatrixReal& A, double b);                            //!< operator * for matrix * scalar 
//    MatrixReal                            operator/(const MatrixReal& A, double b);                            //!< operator / for matrix / scalar 
//    MatrixReal&                           operator/=(MatrixReal& A, double b);                                 //!< operator /= for scalar 
    
}

#endif

