#include "Matrix.h"

template<class T>
Matrix<T>::Matrix(int rows, int columns):
	matrix({})
{
	if (rows <= 0 || columns <= 0)
	{
		throw MatrixConstructionError();
	}

	vector<T> default_row = {};
	default_row.resize(columns);

	for (int i = 0; i < columns; i++)
	{
		default_row[i] = T();
	}

	matrix.resize(rows);

	for (int i = 0; i < rows; i++)
	{
		matrix[i] = default_row;
	}
}

template<class T>
Matrix<T>::Matrix(vector<vector<T>> entries):
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
vector<T>& Matrix<T>::operator[](int i)
{
	return matrix[i];
}

template<class T>
int Matrix<T>::n_rows()
{
	return matrix.size();
}

template<class T>
int Matrix<T>::n_columns()
{
	return matrix[0].size();
}

template<class T>
void Matrix<T>::swap_rows(int first, int second)
{
	if (!(first >= 0 && first < matrix.size() && second >= 0 && second < matrix.size()))
	{
		throw MatrixRowOperationError();
	}

	vector<T> temp = matrix[first];

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
		matrix[row_index][j] = T * matrix[row_index][j];
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
Matrix<T>& operator+(const Matrix<T>& A, const Matrix<T>& B)
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
			sum[i][j] = A[i][j] + B[i][j];
		}
	}

	return sum;
}

template<class T>
Matrix<T>& operator*(const Matrix<T>& A, const Matrix<T>& B)
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
				product[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

template<class T>
ostream& operator<<(ostream& os, const Matrix<T>& matrix)
{
	for (int i = 0; i < matrix.n_rows(); i++)
	{
		for (int j = 0; j < matrix.n_columns(); j++)
		{
			os << matrix[i][j] << " ";
		}

		os << '\n';
	}

	return os;
}

template<class T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> transpose(matrix[0].size(), matrix.size());

	for (int i = 0; i < matrix[0].size(); i++)
	{
		for (int j = 0; j < matrix.size(); j++)
		{
			transpose[i][j] = matrix[j][i];
		}
	}

	return transpose; 
}