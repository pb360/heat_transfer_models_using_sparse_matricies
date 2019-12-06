//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
#ifndef SPARSE_HPP
#define SPARSE_HPP

#include <vector>
#include <tuple>

class SparseMatrix
{
  private:
    double length, width, h, T_c, T_h;
    std::vector<int> i_idx;
    std::vector<int> j_idx;
    std::vector<double> a;
    int ncols;
    int nrows;

    std::vector<int> rows;     
    std::vector<int> cols;     
    std::vector<double> vals; 
    std::vector<double> b;  // RHS of A*u = b
    std::vector<double> u;  // solution vector of A*u = b


    /* TODO: Add any additional private data attributes and/or methods  you need */


  public:
    // constructor
    SparseMatrix(std::string inputfile, std::string soln_prefix);

    /* get the value of the boundary for */
    double get_cold_boundary_value(int j, 
                                    int col_count,
                                    double length, 
                                    double T_c);

    /*  */
    void map(std::vector<double> &point);

    /* Method to convert COO matrix to CSR format using provided function */
    void ConvertToCSR();

    
};

#endif /* SPARSE_HPP */
