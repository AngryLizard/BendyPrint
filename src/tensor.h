
#include <Eigen/Core>
#include <vector>

struct Tensor : public std::vector<std::vector<Eigen::MatrixXd>>
{
	Tensor(int32_t rows, int32_t cols);
	Tensor(int32_t rows, int32_t cols, const Eigen::MatrixXd& X);
	Eigen::MatrixXd& operator()(int32_t i, int32_t j);
	const Eigen::MatrixXd& operator()(int32_t i, int32_t j) const;
	Eigen::MatrixXd operator*(const Eigen::MatrixXd& X) const;
	Tensor operator*(const Tensor& X) const;
	Tensor mult(const Eigen::MatrixXd& X) const;
	Tensor operator+(const Tensor& X) const;
	Tensor operator-(const Tensor& X) const;
	double squaredNorm() const;

	Eigen::MatrixXd flatten() const;

	std::string toString() const;
	operator std::string() const;

	int32_t rows;
	int32_t cols;
};