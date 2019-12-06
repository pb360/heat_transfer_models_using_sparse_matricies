//3456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789
# include <cmath>
# include <fstream>
# include <iostream>
# include <vector>

# include "matvecops.hpp"

unsigned int CGSolver(std::vector<double> &val,
                   			 std::vector<int>    &row_ptr,
				             std::vector<int>    &col_idx,
				             std::vector<double> &b,
				             std::vector<double> &u0,
				             double              tol,
							 std::string         soln_prefix)  {

	// initialize variables that will be needed for algorithm 
	std::vector<double> r0((row_ptr.size()-1), 0);   // first residual
	std::vector<double> r1((row_ptr.size()-1), 0);   // second residual
	std::vector<double> p0((row_ptr.size()-1), 0);   // first p
	std::vector<double> A_u0((row_ptr.size()-1), 0); // matrix vector product 
	std::vector<double> A_p0((row_ptr.size()-1), 0); // matrix vector product 
	std::vector<double> beta_p0((row_ptr.size()-1), 0);     // scalar_vector
	std::vector<double> alpha_p0((row_ptr.size()-1), 0);    // scalar_vector
	std::vector<double> alpha_A_p0((row_ptr.size()-1), 0);  //  alpha*p0
	
	// start cg algorithm    note _ indicates multiplication (ie: A*u0 = A_u0)
	A_u0 = CSR_VEC_MULT(val, row_ptr, col_idx, u0); 
	r0   = minus(b, A_u0); // first residual
	double norm_r0 = sqrt(dot(r0,r0)); // get norm_r0
	p0 = r0; // set p0... need initial residual  through iterations

	unsigned int niter = 0; 
	while (niter < b.size() ) { 
		niter++;

		A_p0 = CSR_VEC_MULT(val, row_ptr, col_idx, p0);
		double alpha = dot(r0, r0)/dot(p0, A_p0);
		alpha_p0 = scalar_mult(alpha, p0);
		u0 = add(u0, alpha_p0);
		alpha_A_p0 = scalar_mult(alpha, A_p0);
		r1 = minus(r0, alpha_A_p0);
		double norm_r = sqrt(dot(r1,r1));

		// when tolerance is met return niter and write solution 
		if (norm_r/norm_r0 < tol) { 
			// Write final Solution File 
			std::cout << "writing file \n";
			std::string filename = soln_prefix + std::to_string(niter);
			std::cout << "filename: " << filename << std::endl;

			std::ofstream outfile;
			outfile.open(filename);
			outfile.precision(3);
			for (unsigned int n = 0 ; n < u0.size(); n++) { 
				outfile << u0[n] << " ";
			}
			return niter; 
		}

		double beta = dot(r1, r1)/dot(r0, r0);
		r0 = r1; 
		beta_p0 = scalar_mult(beta, p0);
		
		// set p0 to new value for next iteration 
		p0 = add(r1, beta_p0);

		if (niter % 10 == 0) { 
			// Write final Solution File 
			std::cout << "writing file \n";
			std::string filename = soln_prefix + std::to_string(niter);
			std::cout << "filename: " << filename << std::endl;

			std::ofstream outfile;
			outfile.open(filename);
			outfile.precision(3);
			for (unsigned int n = 0 ; n < u0.size(); n++) { 
				outfile << u0[n] << " ";
			}
		}
	}
	
	return niter;	
}















