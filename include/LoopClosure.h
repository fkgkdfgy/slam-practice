#pragma once
#include "common.h"
#include "PGM.h"
#include "M2MPGM.h"
#include "CSM.h"
#include "HillClimbFiner.h"
#include "ceres/ceres.h"
typedef Submap_Scan::Node_info::status NodeStatus;
ceres::CostFunction * CreateSPAErrorFucntion(const Array3d & relative_pose);

template <typename T>
static std::array<T, 3> ComputeUnscaledError(
	const Array3d & relative_pose, const T * const start,
	const T * const end) {
	// 
	const T cos_theta_i = cos(start[2]);
	const T sin_theta_i = sin(start[2]);
	const T delta_x = end[0] - start[0];
	const T delta_y = end[1] - start[1];

	// 求解相对量 是end到start 的表示
	const T h[3] = { cos_theta_i * delta_x + sin_theta_i * delta_y,
					-sin_theta_i * delta_x + cos_theta_i * delta_y,
					end[2] - start[2] };
	return { {T(relative_pose(0)) - h[0],
			 T(relative_pose(1)) - h[1],
				 T(relative_pose(2)) - h[2]} };
}

template <typename T>
std::array<T, 3> ScaleError(const std::array<T, 3> & error,
	double translation_weight, double rotation_weight) {

	return { {
		error[0] * translation_weight,
		error[1] * translation_weight,
		error[2] * rotation_weight
	} };

}

class SPAErrorFunctoin
{
public:
	SPAErrorFunctoin(const Eigen::Array3d & relative_pose_):
	relative_pose(relative_pose_){}
	template<typename T>
	bool operator()(const T* const start_pose,const T*const end_pose,
	T* residual) const 
	{
		const std::array<T, 3> error =
			ScaleError(ComputeUnscaledError(
				relative_pose,
				start_pose, end_pose),
				translation_weight,
				rotation_weight);
		std::copy(std::begin(error), std::end(error), residual);
		return true;
	}
	static void UpdateWeight(const double translation_weight_,const double rotation_weight_)
	{
		translation_weight = translation_weight_;
		rotation_weight = rotation_weight_;
		initializationFlag = 1;
	}
private:
	Eigen::Array3d relative_pose;
	static bool initializationFlag;
	static double translation_weight;
	static double rotation_weight;
};

class LoopClosure
{
public:
	struct LoopResult
	{
		LoopResult() { Good = 0; }
		LoopResult(float score_, bool Good_,
			Array3d & BBM_match_pose_) :score(score_), Good(Good_),
		BBM_match_pose(BBM_match_pose_){}
		float score;
		bool Good;
		Array3d BBM_match_pose;
	};
	LoopClosure() { submap_num = 0; LoopCount = 0; }
	void Insert(const PGMBuilder* temp);
	LoopResult BBMMatch(CSM & csmer_,M2MPGM & m2mpgm_,const vector<Array2d> & laserdata_
					   ,const Array3d & initial_pose, const bool world_flag);

	// ���������������Node ��֮ǰ����submap֮������λ��
	void ComputeConstraintForNode(CSM& csmer_, const vector<Array2d>& laserdata_
		, const Array3d& initial_pose, map<int, Submap_Scan>& SubmapScanInfoPool_, int ScanIdx);
	bool CheckLoopFlag();
	bool CheckSolveFlag() const  { return Ready; }
	void Solve(map<int, Submap_Scan>& SubmapScanInfoPool_);
	const map<int,Eigen::Array3d> & GetSolveResult()const {return node_pose_result;}
private:
	vector<const PGMBuilder*> PGMPool;
	vector<M2MPGM*> M2MPGMPool;
	int submap_num;
	bool Ready;
	int LoopCount;
	Array3d LastLoopPose;
	bool NotFind;
	float distance_drive;
	Array3d CurrentlyScanPose;
	map<int,Eigen::Array3d> submap_pose_result;
	map<int,Eigen::Array3d> node_pose_result;
};

ceres::CostFunction * CreateSPAErrorFucntion(const Array3d & relative_pose)
{
	return new ceres::AutoDiffCostFunction<SPAErrorFunctoin,
	3, // residual
	3, // start pose
	3  // end pose
	>(new SPAErrorFunctoin(relative_pose));
}