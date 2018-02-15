#include <vector>
#include <stdafx.h>
#include <ap.h>
#include "data_types.h"

namespace bandwidth_steering {

#ifndef OPTIMIZATION_PARAMS
#define OPTIMIZATION_PARAMS 
struct optimization_parameters{
	// minimize ||Ax - b||
	// subject to Cx <= C[n+1]
	alglib::real_2d_array A;
	alglib::real_1d_array b;
	alglib::real_2d_array C;
	alglib::integer_1d_array ct; // basically tells the optimizer the i-th constrain has what inequality sign
};
#endif


// Reads the matrix of the extracted traffic matrix
bool read_etm(std::string filename, std::vector<std::vector<double_t>>& traffic_matrix);

void row_normalize_matrix(std::vector<std::vector<double_t>>traffic_matrix);

void configure_constraints(int_t num_groups, 
							double_t link_per_group,
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix);

void rounding_solution(std::vector<std::vector<double_t>>& qp_solution, 
						std::vector<std::vector<int_t>>& traffic_matrix, 
						double_t link_per_group);

void bandwidth_steering(std::string filename, 
						std::vector<std::vector<double_t>>& traffix_matrix, 
						double_t link_per_group,
						double_t threshold);
int main(int argc, char** argv);
}
