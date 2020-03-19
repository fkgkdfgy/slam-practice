#include "CeresOp.h"
void OptimizationProblem::Solve()
{
	double pose[3] = { initial_pose_.x(),
					   initial_pose_.y(),
					   initial_pose_.z()};
	ceres::Problem problem_;
	ceres::Solver::Options options_;
	options_.max_num_iterations = 100;
	options_.linear_solver_type = ceres::DENSE_QR;
	// �����scale_factor_��û�п�ʼʹ��
	problem_.
		AddResidualBlock(
		OccupiedSpaceError::CreateAutoDiffCostFunction(Laserdata_ptr_,scale_factor_,pgm_ptr_)
		,NULL,pose);
	ceres::Solver::Summary summary;
	ceres::Solve(options_, &problem_, &summary);
}

