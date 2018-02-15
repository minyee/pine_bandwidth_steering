#include "bandwidth_steering.h"
#include "matrix_op.h"
#include <iostream>
#include <fstream>

namespace bandwidth_steering {

/**
 * Given a filename, reads in the matrix file and then performs the reading process into
 * a traffic matrix
 **/
bool read_etm(std::string filename, 
				std::vector<std::vector<double_t>>& traffic_matrix, 
				int_t& matrix_size) {

	ifstream file_stream (filename);
	// file stream has to be opened
	if (!file_stream.is_open()) {
		return false;
	}
	int_t matrix_size;
	file_stream >> matrix_size; 
	int_t entry;
	int_t index = 0;
	int_t row = 0;
	int_t col = 0;
	
	// Initialize the traffic matrix sizes so that we don't get
	// access faults later
	for (int_t i = 0; i < matrix_size; i++) {
		std::vector<double> row_entries;
		for (int_t j = 0; j < matrix_size; j++) {
			row_entries.push_back(0);
		}
		traffic_matrix.push_back(row_entries);
	}

	int_t limit = matrix_size * matrix_size;
	while (!file_stream.eof() && (index < limit)) {
		file_stream >> entry;
		row = index / matrix_size;
		col = index % matrix_size;
		auto row_entries = traffic_matrix[row];
		row_entries[col] = entry;
		index++;
	}
	file_stream.close();
	return true;
};

/**
 * Given a matrix, normalizes each row to the value given by the argument normalize_to
 * For example, if normalize_to = 3, then the sum of each row will equal 3
 **/
void row_normalize_matrix(std::vector<std::vector<double_t>>& matrix, double_t normalize_to) {
	for (auto row : matrix) {
		double_t row_sum = 0;
		for (auto entry : row) {
			row_sum += entry;
		}
		for (auto entry : row) {
			entry = entry / row_sum;
		}
	}
	return;
};

/**
 * Configures the objective function for this quadratic program
 * NOTE: Note that this is the translation objective function
 *		where the solution vector is not b, but b - e, where e is
 *		an element in the ETM
 * MORE IMPORTANTLY: the b_vector is not the bandwidth matrix vector 
 * 					but rather the b vector in the objective function of QP
 **/
void configure_objective_function(int_t matrix_size,
									std::string& A_matrix,
									std::string& b_vector,
								  	optimization_parameters* opt_params) {
	A_matrix += "[";
	b_vector += "[";
	for (int i = 0; i < matrix_size ; i++) {
		A_matrix += "[";
		for (int j = 0; j < matrix_size ; j++) {
			if (j < matrix_size - 1) {
				if (i == j) A_matrix += "2, ";
				else A_matrix += "0, ";
			} else {
				if (i == j) A_matrix += "2";
				else A_matrix += "0";
			}
		}
		if (i < matrix_size - 1) {
			A_matrix += "], ";
			b_vector += "0,";	
		} else {
			A_matrix += "]";
			b_vector += "0";
		}
	}
	b_vector += "]";
	A_matrix += "]"; // close the final loop

};

/**
 * Configures the individual constraints for each x variable
 **/
void configure_constraints_individual(int_t matrix_size, 
							double_t link_per_group,
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	int_t iterations = matrix_size * matrix_size;
	for (int_t i = 0; i < iterations; i++) {
		C_str += "[";
		for (int_t j = 0; j < iterations; j++) {
			if (j < iterations - 1) {
				if (i == j) C_str += "1, ";
				else C_str += "0, ";
			} else {
				if (i == j) C_str += "1";
				else C_str += "0";
			}
		}
		int_t row = i / matrix_size;
		int_t col = i % matrix_size;
		double_t remainder = link_per_group - traffic_matrix[row][col];
		C_str += (std::to_string(remainder) + "], ");
		d_str += "-1, ";
	}
};

void configure_constraints_rows(int_t matrix_size, 
							double_t link_per_group;
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	int_t iterations = matrix_size * matrix_size;
	// NOTE: write one row for each element in the x vector
	for (int_t row = 0; row < iterations; row++) {
		C_str += "[";
		for (int_t i = 0; i < iterations; i++) {
			if ((i / matrix_size) == row)
				C_str += "1, ";
			else 
				C_str += "0, ";
		}
		C_str += "0], ";
		d_str += "-1, ";
	}
};

void configure_constraints_cols(int_t matrix_size, 
							double_t link_per_group;
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	int_t iterations = matrix_size * matrix_size;
	for (int_t col = 0; col < iterations; col++) {
		C_str += "[";
		double_t col_sum = 0;
		for (int_t i = 0; i < iterations; i++) {
			if ((i % matrix_size) == col) {
				int_t row = i / matrix_size;
				int_t col = i % matrix_size;
				col_sum += traffic_matrix[row][col];
				C_str += "1, ";
			} else {
				C_str += "0, ";
			}
		}
		if (col < iterations - 1) {
			C_str += (std::to_string(link_per_group - col_sum) + "], ");
			d_str += "-1, ";
		} else {
			C_str += (std::to_string(link_per_group - col_sum) + "]");
			d_str += "-1";
		}
	}
};

bool begin_optimization(int_t matrix_size, 
						optimization_parameters* opt_params,
						std::vector<std::vector<double_t>>& solution) {
	alglib::minqpstate state;
    alglib::minqpreport rep;
    alglib::real_1d_array x;
    // create solver, set quadratic/linear terms
    alglib::minqpcreate(num_groups*num_groups, state);
    alglib::minqpsetquadraticterm(state, opt_params->A);
    alglib::minqpsetlinearterm(state, opt_params->b);
    alglib::minqpsetlc(state, opt_params->C, opt_params->ct);

    // Set scale of the parameters.
    // It is strongly recommended that you set scale of your variables.
    // Knowing their scales is essential for evaluation of stopping criteria
    // and for preconditioning of the algorithm steps.
    // You can find more information on scaling at http://www.alglib.net/optimization/scaling.php
    //alglib::minqpsetscale(state, s);

    //
    // Solve problem with BLEIC-based QP solver.
    //
    // This solver is int64_tended for problems with moderate (up to 50) number
    // of general linear constraint64_ts and unlimited number of box constraint64_ts.
    //
    // Default stopping criteria are used.
    //
    alglib::minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
    alglib::minqpoptimize(state);
    alglib::minqpresults(state, x, rep);
    if (rep.terminationtype < 0) {
    	return false;
    }
    std::vector<double> traffic_vector;
    const double* psol = x.getcontent();
    int_t iterations = matrix_size * matrix_size;
    for (int_t i = 0; i < iterations; i++) {
    	int_t row = i / matrix_size;
    	int_t col = i % matrix_size;
    	solution[row][col] = (double_t) psol[i];
    }
    return true;
};

void configure_constraints(int_t matrix_size, 
							double_t link_per_group;
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	C_str += "[";
	d_str += "[";
	configure_constraints_individual(matrix_size, link_per_group, opt_params, C_str, d_str, traffic_matrix);
	configure_constraints_rows(matrix_size, link_per_group, opt_params, C_str, d_str);
	configure_constraints_cols(matrix_size, link_per_group, opt_params, C_str, d_str, traffic_matrix);
	C_str += "]";
	d_str += "]";
	return;
};

/**
 * The main bandwidth steering algorithm.
 * This is the highest level of function call. 
 * NOTE: the argument link_per_group defines the 
 * 			constraint that is going to be used for the 
 * 			number of incoming and outgoing links each group has
 **/
void bandwidth_steering(std::string filename, 
						std::vector<std::vector<double_t>>& traffix_matrix, 
						double_t link_per_group) {

	std::vector<std::vector<double_t>> traffic_matrix;
	int_t matrix_size;
	if (!read_etm(filename, traffic_matrix, matrix_size)) {
		std::cerr << "file cannot be opened" << std::endl;
		std::exit(-1);
	}
	if (link_per_group < 0) {
		link_per_group = (double_t) matrix_size;
	}

	// Step 0: gotta do some preprocessing by row normalizing each row
	row_normalize_matrix(traffic_matrix, link_per_group);
	// if we are here, that means that the traffic matrix has been read in
	optimization_parameters opt_params;
	std::string A_str;
	std::string b_str;
	std::string C_str;
	std::string d_str;

	// Step 1: configure the constraints and objective
	configure_objective_function(matrix_size, A_str, b_str, opt_params);
	configure_constraints(matrix_size, link_per_group, C_str, d_str, opt_params, traffic_matrix);
	opt_params->A = opt_params->C = alglib::real_2d_array(A_str.c_str());
	opt_params->b = alglib::real_1d_array(b_str.c_str());
	opt_params->C = alglib::real_2d_array(C_str.c_str());
	opt_params->ct = alglib::integer_1d_array(d_str.c_str());

	// Step 2: begin the quadratic program
	std::vector<std::vector<double_t>> solution;
	if (!begin_optimization(matrix_size, opt_params, solution)) {
		std::cerr << "Solution did not converge" << std::endl;
		std::exit(-1);
	}
	// Step 3: translate the results into bandwidth matrix
	std::vector<std::vector<double_t>> actual_solution;
	translate_results(solution, traffic_matrix, actual_solution);
	rounding_solution(actual_solution, traffic_matrix);
	return;
};


int main(int argc, std::string[] argv) {
	int default_mode;
	if (argc > 2) {
		default_mode = 1; 
	} else {
		default_mode = 0;
	}
	std::string filename = argv[1];
	std::vector<std::vector<double_t>> traffic_matrix;
	double_t link_per_group;
	if (default_mode == 0) {
		link_per_group = -1;
	} else {
		link_per_group = std::stoi(argv[2]);
	}
	bandwidth_steering(filename, traffic_matrix, link_per_group);
	std::cout << "Exited Cleanly" << std::endl;
	return 0;
};

}