#pragma once

#include <iostream>
#include <vector>
#include <cmath>

/*
using std::vector;
using std::ostream;
using std::cout;
using std::endl;
using std::abs;
*/

template<class T> 
class Matrix
{

private:
	std::vector<std::vector<T>> matrix;

public:
	std::vector<T>& operator[](int i);
	std::vector<T> at(int i) const; 
	T at(int i, int j) const;

	int n_rows() const;
	int n_columns() const;

	bool square() const;

	Matrix<T>(int rows, int columns);
	Matrix<T>(std::vector<std::vector<T>> entries);

	//elementary row operations
	void swap_rows(int first, int second);
	void multiply_row(int row_index, T coeff);
	void add_multiple_of_row(int sum_index, T coeff, int summand_index);

	static Matrix<T> identity(int n);

	Matrix<T> transpose() const;
	Matrix<T> inverse() const; //will return 0 matrix if not invertible

	Matrix<T> row_echelon() const;
	Matrix<T> reduced_row_echelon() const;

	std::vector<std::vector<T>> kernel() const; //returns a basis of the kernel

	T determinant() const;
	T trace() const;
	int rank() const;
	int nullity() const;

	T jacobi_infinity_norm() const;
	T sassenfeld_const() const; 
};

class MatrixConstructionError : public std::exception
{
	virtual	const char* what() const throw()
	{
		return "matrix construction error";
	}
};

class MatrixArithmeticError : public std::exception
{
	virtual	const char* what() const throw()
	{
		return "matrix arithmetic error";
	}
};

class MatrixRowOperationError : public std::exception
{
	virtual	const char* what() const throw()
	{
		return "matrix row operation error";
	}
};

class InvalidRRef : public std::exception
{
	virtual const char* what() const throw()
	{
		return "invalid reduced row echelon";
	}
};

class MatrixZeroEntry : public std::exception
{
	virtual const char* what() const throw()
	{
		return "zero entry where one shouldnt be";
	}
};

template<class T>
Matrix<T>::Matrix(int rows, int columns) :
	matrix({})
{
	if (rows <= 0 || columns <= 0)
	{
		throw MatrixConstructionError();
	}

	std::vector<T> default_row = {};

	for (int i = 0; i < columns; i++)
	{
		default_row.push_back(T());
	}

	for (int i = 0; i < rows; i++)
	{
		matrix.push_back(default_row);
	}
}

template<class T>
Matrix<T>::Matrix(std::vector<std::vector<T>> entries) :
	matrix(entries)
{
	if (entries.size() == 0 || entries[0].size() == 0)
	{
		throw MatrixConstructionError();
	}

	for (int i = 1; i < matrix.size(); i++)
	{
		if (matrix[i].size() != matrix[0].size())
		{
			throw MatrixConstructionError();
		}
	}
}

template<class T>
std::vector<T> Matrix<T>::at(int i) const
{
	return matrix[i];
}

template<class T>
T Matrix<T>::at(int i, int j) const
{
	return matrix[i][j];
}

template<class T>
std::vector<T>& Matrix<T>::operator[](int i)
{
	return matrix[i];
}

template<class T>
int Matrix<T>::n_rows() const
{
	return matrix.size();
}

template<class T>
int Matrix<T>::n_columns() const
{
	return matrix[0].size();
}

template<class T>
bool Matrix<T>::square() const
{
	return matrix.size() == matrix[0].size();
}

template<class T>
void Matrix<T>::swap_rows(int first, int second)
{
	if (!(first >= 0 && first < matrix.size() && second >= 0 && second < matrix.size()))
	{
		throw MatrixRowOperationError();
	}

	std::vector<T> temp = matrix[first];

	matrix[first] = matrix[second];
	matrix[second] = temp;
}

template<class T>
void Matrix<T>::multiply_row(int row_index, T coeff)
{
	if (row_index < 0 || row_index >= matrix.size())
	{
		throw MatrixRowOperationError();
	}

	for (int j = 0; j < matrix[row_index].size(); j++)
	{
		matrix[row_index][j] = coeff * matrix[row_index][j];
	}
}

template<class T>
void Matrix<T>::add_multiple_of_row(int sum_index, T coeff, int summand_index)
{
	if (!(sum_index >= 0 && sum_index < matrix.size() && summand_index >= 0 && summand_index < matrix.size()))
	{
		throw MatrixRowOperationError();
	}

	for (int j = 0; j < matrix[sum_index].size(); j++)
	{
		matrix[sum_index][j] = matrix[sum_index][j] + coeff * matrix[summand_index][j];
	}
}

template<class T>
Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B)
{
	if (A.n_rows() != B.n_rows() || A.n_columns() != B.n_columns())
	{
		throw MatrixArithmeticError();
	}

	Matrix<T> sum(A.n_rows(), A.n_columns());

	for (int i = 0; i < A.n_rows(); i++)
	{
		for (int j = 0; j < A.n_columns(); j++)
		{
			sum[i][j] = A.at(i,j) + B.at(i,j);
		}
	}

	return sum;
}

template<class T>
Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B)
{
	if (A.n_columns() != B.n_rows())
	{
		throw MatrixArithmeticError();
	}

	Matrix<T> product(A.n_rows(), B.n_columns());

	for (int i = 0; i < product.n_rows(); i++)
	{
		for (int j = 0; j < product.n_columns(); j++)
		{
			for (int k = 0; k < A.n_columns(); k++)
			{
				product[i][j] += A.at(i, k) * B.at(k, j);
			}
		}
	}

	return product;
}

template<class T>
Matrix<T> operator*(const T& c,  const Matrix<T>& A)
{
	Matrix<T> product(A.n_rows(), A.n_columns());

	for (int i = 0; i < product.n_rows(); i++)
	{
		for (int j = 0; j < product.n_columns(); j++)
		{
			product[i][j] = c * A.at(i, j);
		}
	}

	return product; 
}

template<class T>
std::vector<T> operator*(const Matrix<T>& A, std::vector<T>& x)
{
	if (x.size() != A.n_columns())
	{
		throw MatrixArithmeticError();
	}

	std::vector<T> product(A.n_rows());

	for (int i = 0; i < product.size(); i++)
	{
		for (int j = 0; j < x.size(); j++)
		{
			product[i] += A.at(i,j) * x[j];
		}
	}

	return product;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& matrix) 
{
	for (int i = 0; i < matrix.n_rows(); i++)
	{
		for (int j = 0; j < matrix.n_columns(); j++)
		{
			os << matrix.at(i,j) << " ";
		}

		os << '\n';
	}

	return os;
}

template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
	os << "{";
	for (int i = 0; i < vec.size(); i++)
	{
		os << vec[i] << ", ";
	}
	os << "}";

	return os;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const
{
	Matrix<T> transpose(matrix[0].size(), matrix.size());

	for (int i = 0; i < matrix[0].size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			transpose[i][j] = matrix.at(j,i);
		}
	}

	return transpose;
}

template<class T>
Matrix<T> Matrix<T>::inverse() const
{
	Matrix<T> inv(n_rows(), n_columns());

	if (n_rows() != n_columns()) //not square, so not invertible
	{
		return inv; 
	}

	Matrix<T> aug(n_rows(), 2 * n_columns());

	for (int i = 0; i < n_rows(); i++)
	{
		for (int j = 0; j < n_columns(); j++)
		{
			aug[i][j] = at(i, j);
		}

		aug[i][i + n_rows()] = 1;
	}

	//cout << "augmented matrix : " << endl << aug << endl;

	Matrix<T> aug_reduced = aug.reduced_row_echelon();

	for (int i = 0; i < n_rows(); i++)
	{
		if (aug_reduced[i][i] != 1) //indicates matrix is not invertible
		{
			return inv;
		}
	}

	for (int i = 0; i < n_rows(); i++)
	{
		for (int j = 0; j < n_columns(); j++)
		{
			inv[i][j] = aug_reduced.at(i, j + n_columns());
		}
	}
	
	return inv;

}

template<class T>
Matrix<T> Matrix<T>::row_echelon() const
{
	Matrix<T> row_echelon = matrix;

	int top = 0;
	int current_col = 0;

	while (current_col < n_columns() && top < n_rows())
	{
		int i = top;

		while (current_col < n_columns() && row_echelon[i][current_col] == 0)
		{
			i++;

			if (i == n_rows() && current_col < n_columns())
			{
				current_col++;
				i = top;
			}
		}

		if (current_col == n_columns())
		{
			return row_echelon;
		}

		row_echelon.swap_rows(top, i);

		T pivot = row_echelon[top][current_col];

		for (int j = top + 1; j < n_rows(); j++)
		{
			row_echelon.add_multiple_of_row(j, -row_echelon[j][current_col] / pivot, top);
			row_echelon[j][current_col] = 0;
		}

		top++;
		current_col++;
	}

	return row_echelon;
}

template<class T>
Matrix<T> Matrix<T>::reduced_row_echelon() const
{
	Matrix<T> r_row_echelon = row_echelon();

	int current_row = 0;
	int current_col = 0;

	while (current_row < n_rows() && current_col < n_columns())
	{
		while (current_col < n_columns() && r_row_echelon[current_row][current_col] == 0)
		{
			current_col++;
		}

		if (current_col == n_columns()) //have found a row of all 0's, so finished
		{
			return r_row_echelon;
		}

		r_row_echelon.multiply_row(current_row, static_cast<T>(1) / r_row_echelon[current_row][current_col]);
		r_row_echelon[current_row][current_col] = static_cast<T>(1);

		for (int i = 0; i < n_rows(); i++)
		{
			if (i != current_row)
			{
				r_row_echelon.add_multiple_of_row(i, -r_row_echelon[i][current_col], current_row);
				r_row_echelon[i][current_col] = static_cast<T>(0);
			}
		}

		current_row++;
		current_col++;
	}

	return r_row_echelon;	
}

template<class T>
Matrix<T> Matrix<T>::identity(int n)
{
	Matrix<T> id(n, n);

	for (int i = 0; i < n; i++)
	{
		id[i][i] = 1;
	}

	return id;
}

template<class T>
std::vector<std::vector<T>> Matrix<T>::kernel() const
{
	Matrix<T> rre_matrix = reduced_row_echelon();
	
	Matrix<T> basis_matrix(n_columns(), n_columns()); // the non-zero columns of this matrix correspond to the kernel basis

	int current_row = 0;
	int current_col = 0;

	while (current_col < n_columns())
	{
		if (current_row == n_rows() || rre_matrix.at(current_row, current_col) == 0) // free variable found
		{
			basis_matrix[current_col][current_col] = 1;
			current_col++;
		}

		else if (rre_matrix.at(current_row, current_col) == 1) //pivot variable found
		{
			for (int i = current_col + 1; i < n_columns(); i++)
			{
				basis_matrix[current_col][i] = -rre_matrix.at(current_row, i);
			}

			current_row++;
			current_col++;
		}

		else //invalid rref
		{
			throw InvalidRRef();
		}
	}

	std::vector<std::vector<T>> basis = {};

	for (int j = 0; j < basis_matrix.n_columns(); j++)
	{
		std::vector<T> column_vec = {};
		bool zero_column = true;

		for (int i = 0; i < basis_matrix.n_rows(); i++)
		{
			column_vec.push_back(basis_matrix.at(i, j));

			if (basis_matrix.at(i, j) != 0)
			{
				zero_column = false;
			}
		}

		if (!zero_column)
		{
			basis.push_back(column_vec);
		}
	}

	return basis;
}

template<class T>
T Matrix<T>::trace() const
{
	if (n_rows() != n_columns())
	{
		throw MatrixArithmeticError();
	}

	T trace;

	for (int i = 0; i < n_rows(); i++)
	{
		trace += matrix[i][i];
	}

	return trace; 
}

template<class T>
T Matrix<T>::jacobi_infinity_norm() const
{
	T max = -1;

	for (int j = 0; j < n_rows(); j++)
	{
		if (at(j, j) == 0)
		{
			throw MatrixZeroEntry();
		}

		T row_sum = 0;
		for (int k = 0; k < n_columns(); k++)
		{
			if (k != j)
			{
				row_sum += abs(at(j, k));
			}
		}

		row_sum /= abs(at(j, j));

		if (row_sum > max)
		{
			max = row_sum;
		}
	}

	return max; 
}

template<class T>
T Matrix<T>::sassenfeld_const() const
{
	std::vector<T> p(n_rows());

	if (at(0, 0) == 0)
	{
		throw MatrixZeroEntry();
	}

	for (int k = 1; k < n_rows(); k++)
	{
		p[0] += abs(at(0, k));
	}

	p[0] /= abs(at(0, 0));

	T p_max = p[0];

	for (int j = 1; j < n_rows(); j++)
	{

		if (at(j, j) == 0)
		{
			throw MatrixZeroEntry();
		}

		p[j] = 0;

		for (int k = 0; k < j; k++)
		{
			p[j] += abs(at(j, k)) * p[k];
		}

		for (int k = j + 1; k < n_rows(); k++)
		{
			p[j] += abs(at(j, k));
		}

		p[j] /= abs(at(j, j));

		if (p[j] > p_max)
		{
			p_max = p[j];
		}
	}

	return p_max;
}