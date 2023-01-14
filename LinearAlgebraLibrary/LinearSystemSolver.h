#pragma once

#include "Matrix.h"
#include "VectorFunctions.h"

template <class T>
class LinearSystemSolver
{
private:

	static std::vector<T> jacobi_iterate(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start);
	static std::vector<T> gauss_seidel_iterate(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start);
	static std::vector<T> gauss_SOR_iterate(const Matrix<T>& A, const std::vector<T>& b, float relax_coeff, const std::vector<T>& start);

public:
	static std::vector<T> solve_row_reduction(const Matrix<T>& A, const std::vector<T>& b); //returns a particular solution to Ax = b, or {} if none exists
	static std::vector<T> solve_jacobi(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float err, int max_iter = 100); //iteratively solves Ax = b, starting at start and continuing until error bound is < err
	static std::vector<T> solve_jacobi_direct(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, int max_iter = 100);
	static std::vector<T> solve_gauss_seidel(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float err, int max_iter = 100);
	static std::vector<T> solve_gauss_SOR(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float relax_coeff, int max_iter = 100);

};

class InvalidLinearSystem : public std::exception
{
	virtual const char* what() const throw()
	{
		return "invalid system passed";
	}
};

template<class T>
std::vector<T> LinearSystemSolver<T>::solve_row_reduction(const Matrix<T>& A, const std::vector<T>& b)
{
	if (A.n_rows() != b.size())
	{
		throw InvalidLinearSystem();
	}

	Matrix<T> aug(A.n_rows(), A.n_columns() + 1);

	for (int i = 0; i < A.n_rows(); i++)
	{
		for (int j = 0; j < A.n_columns(); j++)
		{
			aug[i][j] = A.at(i, j);
		}
		aug[i][A.n_columns()] = b[i];
	}

	Matrix<T> aug_rre = aug.reduced_row_echelon();

	std::vector<T> soln_vector(A.n_columns());

	int current_col = 0;
	int current_row = 0;

	while (current_col < A.n_columns())
	{
		if (current_row == aug_rre.n_rows() || aug_rre.at(current_row, current_col) == 0) // free variable found
		{
			soln_vector[current_col] = 0;
			current_col++;
		}

		else if (aug_rre.at(current_row, current_col) == 1) //pivot variable found
		{
			soln_vector[current_col] = aug_rre[current_row][aug_rre.n_columns() - 1];

			current_row++;
			current_col++;
		}

		else //invalid rref
		{
			throw InvalidRRef();
		}
	}

	while (current_row < aug_rre.n_rows())
	{
		if (aug_rre[current_row][aug_rre.n_columns() - 1] != 0)
		{
			return {};
		}
		current_row++; 
	}
	return soln_vector;
}

template<class T>
std::vector<T> LinearSystemSolver<T>::solve_jacobi(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float err, int max_iter)
{
	if (!A.square() || b.size() != A.n_rows() || start.size() != b.size() || err <= 0 || max_iter <= 0)
	{
		throw InvalidLinearSystem();
	}

	T q = A.jacobi_infinity_norm(); //contraction constant

	if (q >= 1)
	{
		std::cout << "Warning, Jacobi iteration may not converge for matrix " << std::endl << A << std::endl;
	}

	std::vector<T> out = jacobi_iterate(A, b, start);

	float init_jump = p_norm(start - out);

	float current_error = q * init_jump / (1 - q);

	int iterations = 1;

	while (current_error > err && iterations < max_iter)
	{
		out = jacobi_iterate(A, b, out);
		current_error *= q;
		iterations++;
	}

	std::cout << "current error: " << current_error << ", number of iterations: " << iterations << std::endl;

	return out;
}

template<class T>
std::vector<T> LinearSystemSolver<T>::solve_jacobi_direct(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, int max_iter)
{
	if (!A.square() || b.size() != A.n_rows() || start.size() != b.size() || max_iter <= 0)
	{
		throw InvalidLinearSystem();
	}

	std::vector<T> out = start;

	for (int i = 0; i < max_iter; i++)
	{
		out = jacobi_iterate(A, b, out);
	}

	return out;
}


template<class T>
std::vector<T> LinearSystemSolver<T>::jacobi_iterate(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start)
{
	std::vector<T> out(start.size());

	for (int j = 0; j < start.size(); j++)
	{
		if (A.at(j, j) == 0)
		{
			throw MatrixZeroEntry();
		}

		out[j] = b[j];

		for (int k = 0; k < A.n_rows(); k++)
		{
			if (j != k)
			{
				out[j] -= A.at(j, k) * start[k];
			}
		}

		out[j] = out[j]/A.at(j, j);
	}

	return out;
}


template<class T>
std::vector<T> LinearSystemSolver<T>::solve_gauss_seidel(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float err, int max_iter)
{
	if (!A.square() || b.size() != A.n_rows() || start.size() != b.size() || err <= 0 || max_iter <= 0)
	{
		throw InvalidLinearSystem();
	}

	T p = A.sassenfeld_const();

	std::cout << "Sassenfeld constant: " << p << std::endl;

	if (p >= 1)
	{
		std::cout << "Warning; Gauss-Seidel method may not converge for " << std::endl << A << std::endl;
	}

	std::vector<T> out = gauss_seidel_iterate(A, b, start);
	
	//cout << "first: " << out << endl;

	float init_jump = p_norm(start - out);

	float current_error = p * init_jump / (1 - p);

	int iterations = 1;

	while (current_error > err && iterations < max_iter)
	{
		out = gauss_seidel_iterate(A, b, out);
		current_error *= p;
		iterations++;
	}

	std::cout << "current error: " << current_error << ", number of iterations: " << iterations << std::endl;

	return out;
}

template<class T>
std::vector<T> LinearSystemSolver<T>::gauss_seidel_iterate(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start)
{
	std::vector<T> out(start.size());

	for (int j = 0; j < A.n_rows(); j++)
	{
		if (A.at(j, j) == 0)
		{
			throw MatrixZeroEntry();
		}

		out[j] = b[j];

		for (int k = 0; k < j; k++)
		{
			out[j] -= A.at(j, k) * out[k];
		}

		for (int k = j + 1; k < A.n_rows(); k++)
		{
			out[j] -= A.at(j, k) * start[k];
		}

		out[j] /= A.at(j, j);

	}

	return out;
}

template <class T>
std::vector<T> LinearSystemSolver<T>::solve_gauss_SOR(const Matrix<T>& A, const std::vector<T>& b, const std::vector<T>& start, float relax_coeff, int max_iter)
{
	if (!A.square() || b.size() != A.n_rows() || start.size() != b.size() || max_iter <= 0)
	{
		throw InvalidLinearSystem();
	}

	if (relax_coeff <= 0 || relax_coeff >= 2)
	{
		std::cout << "warning; relaxation coefficient should be between 0 and 2";
	}

	std::vector<T> out = gauss_SOR_iterate(A, b, relax_coeff, start);
	int iterations = 1;

	while (iterations < max_iter)
	{
		out = gauss_SOR_iterate(A, b, relax_coeff, out);
		iterations++;
	}

	return out; 
}

template<class T>
std::vector<T> LinearSystemSolver<T>::gauss_SOR_iterate(const Matrix<T>& A, const std::vector<T>& b, float relax_coeff, const std::vector<T>& start)
{
	std::vector<T> out(start.size());

	for (int j = 0; j < start.size(); j++)
	{	
		if (A.at(j, j) == 0)
		{
			throw MatrixZeroEntry();
		}

		out[j] = b[j];

		for (int k = 0; k < j; k++)
		{
			out[j] -= A.at(j, k) * out[k];
		}

		for (int k = j; k < start.size(); k++)
		{
			out[j] -= A.at(j, k) * start[k];
		}

		out[j] *= (relax_coeff / A.at(j, j));
		out[j] += start[j];
	}

	return out; 
}