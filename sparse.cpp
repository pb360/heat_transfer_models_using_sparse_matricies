//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
#include <array>
#include <cmath>
#include <fstream>
#include <iostream> 
#include <string> 

#include "sparse.hpp"
#include "COO2CSR.hpp"
#include "CGSolver.hpp"
#include "matvecops.hpp"

SparseMatrix::SparseMatrix(std::string inputfile, std::string soln_prefix) {
    // read input file and get parameters describing matrix
    double length, width, h, T_c, T_h;

    std::ifstream infile;
    infile.open(inputfile);

    infile >> length >> width >> h;
    infile >> T_c    >> T_h;

    nrows = int (width/h - 1);
    ncols = int (length/h); 

    // Fill Solution vector (u) and RHS vector (b)
    std::vector<double> b_temp(( (unsigned long) (nrows*ncols)), 0);  
    b = b_temp; 
    std::vector<double> u_temp((nrows*ncols), 1);  
    u = u_temp; 
    
    // Populate COO matrix & RHS vectors:  vals,  rows,  cols,  b 
    for (int j = 0; j < ncols; j++) {
        for (int i = 0; i < nrows; i++ ) {
        // initializing indicy values to add to matrix 
        // visual of what each point indicy suffix maps too
        //        p1
        //  p4  (i,j)   p2
        //        p3
        int row = nrows*j+i; // current row of matrix were on

            // top row... no p1 as --> T_h value to RHS 
            if (i == 0) {
                int i2, i3, i4; 
                int j2, j3, j4; 

                i2 = i;     j2 = j+1;  
                i3 = i+1;   j3 = j;  
                i4 = i;     j4 = j-1;  

                // top left entry of array... wrap p4 to far right column 
                if (j == 0) { 
                    j4 = ncols-1;
                }

                // top right of "array"... wrap p2 to far left, column 0
                if (j == ncols-1) { 
                    j2 = 0;
                }

                // Move hot boundary value to RHS 
                b[j*nrows+i] += T_h;

                // (i, j) --> (I, J)... map grid to heat matrix coordinates
                int J2 = nrows*j2 + i2; 
                int J3 = nrows*j3 + i3; 
                int J4 = nrows*j4 + i4; 

                // append these values to vectors 
                rows.push_back(row); rows.push_back(row); rows.push_back(row);
                cols.push_back(J2); cols.push_back(J3); cols.push_back(J4);
                vals.push_back(-1);  vals.push_back(-1);  vals.push_back(-1);
            }

            // bottom row... no p3 --> T_c value to RHS 
            else if (i == nrows-1) {
                int i1, i2, i4; 
                int j1, j2, j4; 

                i1 = i-1;   j1 = j; 
                i2 = i;     j2 = j+1;
                i4 = i;     j4 = j-1;

                // bottom left of array... wrap p4 to far right column 
                if (j == 0) { 
                    j4 = ncols-1;
                }
                // bottom right of "array"... wrap p2 to far left, column 0
                if (j == ncols-1) { 
                    j2 = 0;
                }

                // Move cold boundary value to RHS
                double T_c_j = get_cold_boundary_value(j, ncols, length, T_c);
                b[j*nrows+i] += T_c_j;

                // (i, j) --> (I, J)... map grid to heat matrix coordinates
                int J1 = nrows*j1 + i1; 
                int J2 = nrows*j2 + i2; 
                int J4 = nrows*j4 + i4; 

                // append these values to vectors 
                rows.push_back(row); rows.push_back(row); rows.push_back(row);
                cols.push_back(J1); cols.push_back(J2); cols.push_back(J4);
                vals.push_back(-1);  vals.push_back(-1);  vals.push_back(-1);
            }

            // if not top or bottom row 
            else { 
                int i1, i2, i3, i4; 
                int j1, j2, j3, j4; 
                
                i1 = i-1;   j1 = j; 
                i2 = i;     j2 = j+1;
                i3 = i+1;   j3 = j; 
                i4 = i;     j4 = j-1;
                
                // Left column of "array"... wrap p4 to far right column 
                if (j == 0) { 
                    j4 = ncols-1;
                }

                // Right column of "array"... wrap p2 to far left, column 0
                if (j == ncols-1) { 
                    j2 = 0;
                }

                // (i, j) --> (I, J)... map grid to heat matrix coordinates
                int J1 = nrows*j1 + i1;
                int J2 = nrows*j2 + i2; 
                int J3 = nrows*j3 + i3; 
                int J4 = nrows*j4 + i4; 
                // append values to vectors 
                rows.push_back(row); rows.push_back(row); 
                rows.push_back(row); rows.push_back(row);
                cols.push_back(J1);  cols.push_back(J2); 
                cols.push_back(J3);  cols.push_back(J4);
                vals.push_back(-1);  vals.push_back(-1);  
                vals.push_back(-1);  vals.push_back(-1);
            }
            // append i-th, j-th cofficient (all cases)
            int J = nrows*j + i;
            rows.push_back(row);
            cols.push_back(J);
            vals.push_back(4);
        } 
    }  

    // COO -->> CSR 
    COO2CSR(vals, rows, cols);

    // Solver 
    double tol = 0.00005;
	int iter_count = CGSolver(vals, rows, cols, b, u, tol, soln_prefix);
}


double SparseMatrix::get_cold_boundary_value(int j, 
                                            int col_count, 
                                            double length, 
                                            double T_c){ 
    double x;
    double minus_ten = -10;
    double J = j;
    double NCOLS = col_count;
    x = J/NCOLS*length;

    double exp_term = minus_ten*(x - length/2)*(x - length/2);
    double cold_boundary_temp = -T_c*(exp(exp_term) - 2);
    return cold_boundary_temp;
}


/* Method to convert COO matrix to CSR format using provided function */
void SparseMatrix::ConvertToCSR(){
	// COO2CSR(vals, rows, cols);
}
