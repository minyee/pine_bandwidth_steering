#include <vector>
#include "bandwidth_steering.h"
#include "data_types.h"
namespace matrix_op {

//typedef matrix_double_t std::vector<std::vector<double_t>>;

void vector_sum(std::vector<double_t>& a, std::vector<double_t>& b, int_t size, std::vector<double_t>& sol){
	sol.resize(size);
	for (int_t i = 0; i < size; i++) {
		sol[i] = (a[i] + b[i]);
	}
};

// calculates a transposed multiplied by b, a.k.a transpose(a) * b.
int_t inner_product(std::vector<int_t>& a, std::vector<int_t>& b, int_t size) {
	int_t sol = 0;
	for (int_t i = 0; i < size; i++) {
		sol += a[i] * b[i];
	}
	return sol;
};

void matrix_vector_multiply(std::vector< std::vector<int_t> >& matrix, 
											std::vector<int_t>& a, 
											int_t num_rows, 
											int_t num_cols, 
											std::vector<int_t>& sol) {
	sol.resize(num_rows);
	for (int_t i = 0; i < num_rows; i++) {
		sol[i] = inner_product(matrix[i],a,num_cols);
	}
};

}