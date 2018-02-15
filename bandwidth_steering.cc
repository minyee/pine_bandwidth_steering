#include "bandwidth_steering.h"
#include "matrix_op.h"
#include <iostream>
#include <fstream>
#include <optimization.h>
#include <algorithm>
namespace bandwidth_steering {

/**
 * Given a filename, reads in the matrix file and then performs the reading process into
 * a traffic matrix
 **/
bool read_etm(std::string filename, 
				std::vector<std::vector<double_t>>& traffic_matrix, 
				int_t& matrix_size) {

	std::ifstream file_stream (filename);
	// file stream has to be opened
	if (!file_stream.is_open()) {
		return false;
	}
	
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
		std::vector<double_t>& row_entries = traffic_matrix[row];
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
	int_t iter = matrix_size * matrix_size;
	for (int i = 0; i < iter ; i++) {
		A_matrix += "[";
		for (int j = 0; j < iter ; j++) {
			if (j < iter - 1) {
				if (i == j) A_matrix += "2, ";
				else A_matrix += "0, ";
			} else {
				if (i == j) A_matrix += "2";
				else A_matrix += "0";
			}
		}
		if (i < iter - 1) {
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
			if (i == j) C_str += "1, ";
			else C_str += "0, ";
		}
		int_t row = i / matrix_size;
		int_t col = i % matrix_size;
		double_t remainder = -traffic_matrix[row][col];
		C_str += (std::to_string(remainder) + "], ");
		d_str += "1, ";
	}
};

void configure_constraints_rows(int_t matrix_size, 
							double_t link_per_group,
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	int_t iterations = matrix_size * matrix_size;
	// NOTE: write one row for each element in the x vector
	for (int_t row = 0; row < matrix_size; row++) {
		C_str += "[";
		for (int_t i = 0; i < iterations; i++) {
			if ((i / matrix_size) == row) {
				C_str += "1, ";
			} else  {
				C_str += "0, ";
			}
		}
		C_str += "0], ";
		d_str += "-1, ";
	}
	//std::cout << C_str << std::endl;
};

/**
 * Configures the column constraints of the 
 **/
void configure_constraints_cols(int_t matrix_size, 
							double_t link_per_group,
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	int_t iterations = matrix_size * matrix_size;
	for (int_t col = 0; col < matrix_size; col++) {
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
		if (col < matrix_size - 1) {
			C_str += (std::to_string(link_per_group - col_sum) + "], ");
			d_str += "-1, ";
		} else {
			C_str += (std::to_string(link_per_group - col_sum) + "]");
			d_str += "-1";
		}
	}
	
};

/**
 * Initializes the optimization process
 **/
bool begin_optimization(int_t matrix_size, 
						optimization_parameters* opt_params,
						std::vector<std::vector<double_t>>& solution) {
	alglib::minqpstate state;
    alglib::minqpreport rep;
    alglib::real_1d_array x;
    // create solver, set quadratic/linear terms
    alglib::minqpcreate(matrix_size*matrix_size, state);
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
    solution.resize(matrix_size);
    for (int_t i = 0; i < matrix_size; i++) 
		solution[i].resize(matrix_size);

    for (int_t i = 0; i < iterations; i++) {
    	int_t row = i / matrix_size;
    	int_t col = i % matrix_size;
    	solution[row][col] = (double_t) psol[i];
    }
    return true;
};

void configure_constraints(int_t matrix_size, 
							double_t link_per_group,
							std::string& C_str, 
							std::string& d_str, 
							optimization_parameters* opt_params, 
							std::vector<std::vector<double_t>>& traffic_matrix) {
	C_str += "[";
	d_str += "[";
	configure_constraints_individual(matrix_size, link_per_group, C_str, d_str, opt_params, traffic_matrix);
	configure_constraints_rows(matrix_size, link_per_group, C_str, d_str, opt_params, traffic_matrix);
	configure_constraints_cols(matrix_size, link_per_group, C_str, d_str, opt_params, traffic_matrix);
	C_str += "]";
	d_str += "]";
	std::cout << C_str << std::endl;
	return;
};

/**
 * Translates the result into the actual bandwidth allocation matrix form which result is stored
 * in the actual_solution matrix
 **/
void translate_results(int_t matrix_size, 
						std::vector<std::vector<double_t>>& solution,
						std::vector<std::vector<double_t>>& traffic_matrix,
						std::vector<std::vector<double_t>>& actual_solution) {
	actual_solution.resize(matrix_size);
	matrix_op::square_matrix_sum(matrix_size, solution, traffic_matrix, actual_solution);
	return;
}

/**
 * Called during the final stages of the algorithm to run down results
 * This rounding algorithm is adjustable to be somewhat conservative in that it will
 * set what links will go to 0 if the traffic doesn't exceed the threshold
 **/
void rounding_solution(int_t matrix_size, 
						std::vector<std::vector<double_t>>& traffic_matrix,
						std::vector<std::vector<double_t>>& actual_solution,
						std::vector<std::vector<int_t>>& actual_solution_integer,
						double_t threshold,
						int_t link_per_group) {
	// Step 1: Round everything down
	actual_solution_integer.resize(matrix_size);
	for (int_t i = 0; i < matrix_size; i++) {
		actual_solution_integer[i].resize(matrix_size);
		for (int_t j = 0; j < matrix_size; j++) {
			actual_solution_integer[i][j] = (int_t) traffic_matrix[i][j];
			actual_solution_integer[i][j] = 
				(traffic_matrix[i][j] >= threshold) ? std::max((int_t) actual_solution[i][j], (int_t) 1) : 0;
		}
	}
	std::vector<int_t> row_links_left(matrix_size);
	std::vector<int_t> col_links_left(matrix_size);

	// Step 2: Now scan across to reround actual_solution_integer
	for (int_t i = 0; i < matrix_size; i++) {
		int_t col_sum = 0;
		int_t row_sum = 0;
		for (int_t j = 0; j < matrix_size; j++) {
			row_sum += actual_solution_integer[i][j];
			col_sum += actual_solution_integer[j][i];
		}
		row_links_left[i] = link_per_group - row_sum;
		col_links_left[i] = link_per_group - col_sum;
	}
	for (int_t row = 0; row < matrix_size; row++) {
		while (row_links_left[row] != 0) {
			double_t curr_min = 10000;
			int_t col_entry = 0;
			// NOTE: in this case we want to take out links from the row, so find columns that are over occupied
			if (row_links_left[row] < 0) {
				for (int_t col = 0; col < matrix_size; col++) {
					if (actual_solution_integer[row][col] <= 1 || 
							col_links_left[col] >= 0)
						continue; // skip over entries with very few links
									// also skip over column entries which are already occupied fully or underoccupied
					if (curr_min > pow((double_t) actual_solution_integer[row][col] - traffic_matrix[row][col], 2)) {
						col_entry = col;
						curr_min = pow((double_t) actual_solution_integer[row][col] - traffic_matrix[row][col], 2);
					}
				}
				actual_solution_integer[row][col_entry]--;
				col_links_left[col_entry]++;
				row_links_left[row]++;
			} 
			// NOTE: in this case we want to put in links into the row, so find columns that are under occupied
			else {
				for (int_t col = 0; col < matrix_size; col++) {
					if (actual_solution_integer[row][col] <= 1 || 
							col_links_left[col] <= 0)
						continue; // skip over entries with very few links
									// also skip over column entries which are already occupied fully or overoccupied
					if (curr_min > pow((double_t) actual_solution_integer[row][col] - traffic_matrix[row][col], 2)) {
						col_entry = col;
						curr_min = pow((double_t) actual_solution_integer[row][col] - traffic_matrix[row][col], 2);
					}
				}
				actual_solution_integer[row][col_entry]++;
				col_links_left[col_entry]--;
				row_links_left[row]--;
			}
		}
	}
}

void print_matrix(std::vector<std::vector<double_t>>& mat) {
	int_t size = mat.size();
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			std::cout << std::to_string(mat[i][j]) << " ";
		}
		std::cout << std::endl;
	}
}

/**
 * The main bandwidth steering algorithm.
 * This is the highest level of function call. 
 * NOTE: the argument link_per_group defines the 
 * 			constraint that is going to be used for the 
 * 			number of incoming and outgoing links each group has
 **/
void bandwidth_steering(std::string filename, 
						std::vector<std::vector<double_t>>& traffix_matrix, 
						double_t link_per_group, 
						double_t threshold) {
	std::cout << "cp1" << std::endl;
	std::vector<std::vector<double_t>> traffic_matrix;
	int_t matrix_size;
	if (!read_etm(filename, traffic_matrix, matrix_size)) {
		std::cerr << "file cannot be opened" << std::endl;
		std::exit(-1);
	}
	if (link_per_group < 0) {
		link_per_group = (double_t) matrix_size - 1;
	}

	// Step 0: gotta do some preprocessing by row normalizing each row
	row_normalize_matrix(traffic_matrix, link_per_group);
	print_matrix(traffic_matrix);
	std::cout << "cp2" << std::endl;
	// if we are here, that means that the traffic matrix has been read in
	optimization_parameters opt_params;
	std::string A_str;
	std::string b_str;
	std::string C_str;
	std::string d_str;

	// Step 1: configure the constraints and objective
	configure_objective_function(matrix_size, A_str, b_str, &opt_params);
	std::cout << "cp3" << std::endl;
	std::cout << A_str << std::endl;
	std::cout << b_str << std::endl;
	configure_constraints(matrix_size, link_per_group, C_str, d_str, &opt_params, traffic_matrix);
	//std::cout << C_str << std::endl;
	std::cout << d_str << std::endl;
	opt_params.A = alglib::real_2d_array(A_str.c_str());
	opt_params.b = alglib::real_1d_array(b_str.c_str());
	opt_params.C = alglib::real_2d_array(C_str.c_str());
	opt_params.ct = alglib::integer_1d_array(d_str.c_str());
	std::cout << "cp4" << std::endl;
	// Step 2: begin the quadratic program
	std::vector<std::vector<double_t>> solution;
	if (!begin_optimization(matrix_size, &opt_params, solution)) {
		std::cerr << "Solution did not converge" << std::endl;
		std::exit(-1);
	}
	std::cout << "cp5" << std::endl;
	// Step 3: translate the results into bandwidth matrix
	std::vector<std::vector<double_t>> actual_solution;
	translate_results(matrix_size, solution, traffic_matrix, actual_solution);
	std::cout << "cp6" << std::endl;
	std::vector<std::vector<int_t>> final_integer_soln;
	rounding_solution(matrix_size, 
						actual_solution, 
						traffic_matrix, 
						final_integer_soln,
						threshold, 
						(int_t) link_per_group);
	return;
};
}

using namespace std;

int main(int argc, char** argv) {
	int default_mode;
	if (argc > 2) {
		default_mode = 1; 
	} else {
		default_mode = 0;
	}
	double_t threshold = 0.1;
	std::string filename = argv[1];
	std::vector<std::vector<double_t>> traffic_matrix;
	double_t link_per_group;
	if (default_mode == 0) {
		link_per_group = -1;
	} else {
		link_per_group = std::stoi(argv[2]);
	}

	bandwidth_steering::bandwidth_steering(filename, traffic_matrix, link_per_group, threshold);
	std::cout << "Exited Cleanly" << std::endl;
	return 0;
};