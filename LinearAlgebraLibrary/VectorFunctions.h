#pragma once

#include<vector>
#include<cmath>

class VectorArithmeticError : public std::exception
{
	virtual const char* what() const throw()
	{
		return "vector arithmetic error";
	}
};

template<class T>
float p_norm(std::vector<T> vec, float p = -1)
{
	if (p == -1)
	{
		float max = -1;
		for (int i = 0; i < vec.size(); i++)
		{
			if (std::abs(vec[i]) > max)
			{
				max = std::abs(vec[i]);
			}
		}

		return max;
	}

	else
	{
		float out = 0;
		for (int i = 0; i < vec.size(); i++)
		{
			out += std::pow(std::abs(vec[i]), p);
		}

		return std::pow(out, 1 / p);
	}
}

template<class T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	if (a.size() != b.size())
	{
		throw VectorArithmeticError();
	}

	std::vector<T> out;

	for (int i = 0; i < a.size(); i++)
	{
		out.push_back(a[i] + b[i]);
	}

	return out; 
}

template<class T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b)
{
	if (a.size() != b.size())
	{
		throw VectorArithmeticError();
	}

	std::vector<T> out;

	for (int i = 0; i < a.size(); i++)
	{
		out.push_back(a[i] - b[i]);
	}

	return out;
}