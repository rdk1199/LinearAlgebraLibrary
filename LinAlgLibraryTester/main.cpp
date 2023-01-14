
#include <iostream>
#include "LinearSystemSolver.h"

using std::cout;
using std::vector;
using std::endl;

int main()
{

	//Matrix<float> test(5, 5);

	
	Matrix<float> test({{4, 1, .25},
						{2 ,3 ,.5},
						{1, 1.5, 1} });
	
	//vector<vector<float>> kernel = test.kernel();

	vector<float> test_b = { 1, 2, 3};

	vector<float> soln = LinearSystemSolver<float>::solve_jacobi_direct(test, test_b, { 100000,100000121.0001,1000000 }, 1000000);

	//cout << "kernel: " << kernel << endl;

	cout << "soln: " << soln << endl;

	cout << "product: " << test * soln << endl;
}
