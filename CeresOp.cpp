#include "CeresOp.h"
void OptimizationProblem::Solve()
{
	double pose[3] = { initial_pose_.x(),
					   initial_pose_.y(),
					   initial_pose_.z() };
	ceres::Problem problem_;
	ceres::Solver::Options options_;
	options_.max_num_iterations = 100;
	options_.linear_solver_type = ceres::DENSE_QR;
	// 这里的scale_factor_还没有开始使用
	problem_.
		AddResidualBlock(
		OccupiedSpaceError::CreateAutoDiffCostFunction(Laserdata_ptr_,scale_factor_,pgm_ptr_)
		,NULL,pose);
	ceres::Solver::Summary summary;
	ceres::Solve(options_, &problem_, &summary);

}




// BicubicInterpolation 
template<typename T>
void BicubicInerpolator::CubicInterpolation<T>(const T& t1, const T& t2, const T& t3,
	const T& t4,const T & x ,T * estimation,T *estimation_derivation)
{
	const T a = 0.5 * (-p0 + 3.0 * p1 - 3.0 * p2 + p3);
	const T b = 0.5 * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3);
	const T c = 0.5 * (-p0 + p2);
	const T d = p1;

	// f = ax^3 + bx^2 + cx + d
	if (f != NULL) {
		*estimation = d + x * (c + x * (b + x * a));
	}

	// dfdx = 3ax^2 + 2bx + c
	if (dfdx != NULL) {
		*estimation_derivation = c + x * (2.0 * b + 3.0 * a * x);
	}

}