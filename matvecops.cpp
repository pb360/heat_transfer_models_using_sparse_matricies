# include <cmath>
# include <fstream>
# include <iostream>
# include <vector>

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
std::vector<double> add(std::vector<double> &v1, std::vector<double> &v2) {
	std::vector<double> v_new(v1.size(), 0);
	for (unsigned int n = 0; n < v1.size(); n++){ 
		v_new[n] = v1[n] + v2[n]; 
	}
	return v_new;
}

std::vector<double> minus(std::vector<double> &v1, std::vector<double> &v2) {
	std::vector<double> v_new(v1.size(), 0);
	for (unsigned int n = 0; n < v1.size(); n++){ 
		v_new[n] = v1[n] - v2[n]; 
	}
	return v_new;
}

std::vector<double>  scalar_mult(double scalar, std::vector<double> &v){
	std::vector<double> v_new(v.size(), 0);
	for (unsigned int n = 0 ; n < v.size(); n++){ 
		v_new[n] = scalar*v[n]; 
	}
	return v_new;
}

double dot(std::vector<double> &v1, std::vector<double> &v2) {
	double dot = 0;
	for (unsigned int n = 0 ; n < v1.size(); n++){ 
		dot += v1[n]*v2[n]; 
	}
	return dot;
}


std::vector<double> CSR_VEC_MULT( std::vector<double> &val,
                   				  std::vector<int>    &row_ptr,
				                  std::vector<int>    &col_idx,
				                  std::vector<double> &x) {
	// val_index walks through each value in the val vector
	unsigned int val_index = 0;

	// X_new is the vector that is the product... filled in over loops 
	std::vector<double> x_new(x.size(), 0);

	// loop through a rows worth of cols entries or the transformation of 
	// one element of b
	for (unsigned int n = 0; n < x.size() ; n++) { 
		int c_id_0 =  row_ptr[n];
		int c_id_1 =  row_ptr[n+1];

		// new entry for the x to be created 
		double new_entry = 0;

		// loop through number of values in a row corresponding to entry of x
		for (int i =0; i < (c_id_1 - c_id_0); i++) {
			int col = col_idx[val_index];
			double v   = val[val_index];
			double x_n = x[col];

			new_entry +=  v*x_n;
			val_index++; //increment the walk through whole value vector
		}
		
		x_new[n] = new_entry;
	}

	return x_new;	
}



