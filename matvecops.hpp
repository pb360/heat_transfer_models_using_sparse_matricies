#ifndef MATVECOPS_HPP
#define MATVECOPS_HPP

# include <cmath>
# include <fstream>
# include <iostream>
# include <vector>

/* add() 
 * Adds two vectors of equal size.
 */
std::vector<double> add(std::vector<double> &v1, std::vector<double> &v2); 

/* minus() 
 * subtracts two vectos of equal size 
 */
std::vector<double> minus(std::vector<double> &v1, std::vector<double> &v2); 

/* scalar_mult() 
 * multiplies a vector by a scalar
 */
std::vector<double>  scalar_mult(double scalar, std::vector<double> &v); 

/* dot() 
 * finds the dot product of two vectors
 */
double dot(std::vector<double> &v1, std::vector<double> &v2); 

/* CSR_VEC_MULT() 
 * multiplies a matrix in CSR form given with vectors (value, row_ptr, 
 * & col_idx) by another vector, x. Returns the new vector 
 */
std::vector<double> CSR_VEC_MULT( std::vector<double> &val,
                   				  std::vector<int>    &row_ptr,
				                  std::vector<int>    &col_idx,
				                  std::vector<double> &x); 

#endif 