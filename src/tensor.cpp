#include "Tensor.h"


Tensor::Tensor(int32_t rows, int32_t cols)
	: rows(rows), cols(cols)
{
	resize(rows);
	for (int i = 0; i < rows; i++)
	{
		(*this)[i].resize(cols);
	}
}

Tensor::Tensor(int32_t rows, int32_t cols, const Eigen::MatrixXd& X)
	: Tensor(rows, cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			(*this)(i, j) = X;
		}
	}
}

Eigen::MatrixXd& Tensor::operator()(int32_t i, int32_t j)
{
	return (*this)[i][j];
}

const Eigen::MatrixXd& Tensor::operator()(int32_t i, int32_t j) const
{
	return (*this)[i][j];
}

Eigen::MatrixXd Tensor::operator*(const Eigen::MatrixXd& X) const
{
	Eigen::MatrixXd out(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out(i, j) = (*this)(i, j).cwiseProduct(X).sum();
		}
	}
	return out;
}

Tensor Tensor::operator*(const Tensor& X) const
{
	Tensor out(rows, X.cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < X.cols; j++)
		{
			out(i, j) = (*this)(i, 0) * X(0, j);
			for (int k = 1; k < cols; k++)
			{
				out(i, j) += (*this)(i, k) * X(k, j);
			}
		}
	}
	return out;
}

Tensor Tensor::mult(const Eigen::MatrixXd& X) const
{
	Tensor out(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out(i, j) = X * (*this)(i, j).transpose();
		}
	}
	return out;
}

Tensor Tensor::operator+(const Tensor& X) const
{
	Tensor out(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out(i, j) = (*this)(i, j) + X(i, j);
		}
	}
	return out;
}

Tensor Tensor::operator-(const Tensor& X) const
{
	Tensor out(rows, cols);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out(i, j) = (*this)(i, j) - X(i, j);
		}
	}
	return out;
}

double Tensor::squaredNorm() const
{
	double out = 0.0;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			out += (*this)(i, j).squaredNorm();
		}
	}
	return out;
}

Eigen::MatrixXd Tensor::flatten() const
{
	Eigen::MatrixXd out(rows * cols, rows * cols);
	for (int j = 0; j < cols; j++)
	{
		for (int i = 0; i < rows; i++)
		{
			Eigen::MatrixXd M = (*this)(i, j);
			out.row(i + j * rows) = Eigen::Map<Eigen::RowVectorXd>(M.data(), M.size());
		}
	}
	return out;
}


std::string Tensor::toString() const
{
	std::stringstream ss;
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			ss << "." << std::endl << (*this)[i][j] << std::endl;
		}
	}
	return ss.str();
}

Tensor::operator std::string() const
{
	return toString();
}
