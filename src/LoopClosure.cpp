#include "LoopClosure.h"

double SPAErrorFunctoin::translation_weight = 0;
double SPAErrorFunctoin::rotation_weight = 0;
bool SPAErrorFunctoin::initializationFlag =0;

void LoopClosure::Insert(const PGMBuilder* temp)
{
	PGMPool.push_back(temp);
	M2MPGMPool.push_back(new M2MPGM(*temp));
	submap_num++;
}

LoopClosure::LoopResult LoopClosure::BBMMatch(CSM& csmer_,  M2MPGM& m2mpgm_, const vector<Array2d>& laserdata_
	, const Array3d& initial_pose, const bool world_flag)
{	
	// ֮��Ҫ�޸�BBM �Ĵ�ֻ���
	auto result_ = csmer_.BranchAndBoundMatch(m2mpgm_, laserdata_, initial_pose);
	//CSM��generateaWindow֮����Ҫ��һ��
	Array3d BBM_result = result_.bestpose + Array3d(result_.candidate(0)
		, result_.candidate(1), 0) * m2mpgm_.GetPixelSize();
	if (result_.bestscore == 0)
		return LoopResult();
	else
		return LoopResult(result_.bestscore, bool(1), BBM_result);
}

// ����CheckLoopFlag֮���ٽ����������
void LoopClosure::ComputeConstraintForNode(CSM& csmer_,const vector<Array2d>& laserdata_
	, const Array3d& initial_pose, map<int, Submap_Scan>& SubmapScanInfoPool_,int ScanIdx)
{
	if (SubmapScanInfoPool_.size() == 2)
	{
		return;
	}
	auto distance_lastloop = initial_pose - LastLoopPose;
	if (Vector2d(distance_lastloop(0), distance_lastloop(1)).norm() < 6 )
	{
		return;
	}

	if (NotFind)
	{
		auto distance_lastscan = initial_pose - CurrentlyScanPose;
		distance_drive += Vector2d(distance_lastscan(0), distance_lastscan(1)).norm();
		CurrentlyScanPose = initial_pose;
		if (distance_drive < 3)
			return ;
	}
	// SubmapIdx ��0��ʼ��
	int SubmapIdx = M2MPGMPool.size()-1;
	bool FindFlag = 0;
	for (auto& element : M2MPGMPool)
	{
		// ���о�������
		
		auto distance_submap = initial_pose - SubmapScanInfoPool_[SubmapIdx].mean_pose;
		if (Vector2d(distance_submap(0), distance_submap(1)).norm() < 2)
		{
			// ����Ƿ�֮���Ѿ����ӹ���
			auto iter_ = SubmapScanInfoPool_[SubmapIdx].Node_info_pool.find(ScanIdx);
			if (iter_ == SubmapScanInfoPool_[SubmapIdx].Node_info_pool.end())
			{
				LoopResult BBM_match_result = BBMMatch(csmer_, *element, laserdata_, initial_pose, 1);
				// ����ѵ����ʵ�ƥ�� >0.6 ��score
				if (BBM_match_result.Good)
				{
					HillClimbFiner finermatcher(*PGMPool[SubmapIdx], csmer_);
					float matching_score = -1.f;
					Array3d finer_pose = finermatcher.FinerEstimaiton(laserdata_, BBM_match_result.BBM_match_pose, matching_score);
					SubmapScanInfoPool_[SubmapIdx].Insert(SubmapIdx, ScanIdx, NodeStatus(0), finer_pose, matching_score);
					SubmapIdx--;
					FindFlag = 1;
					break;
				}
				else
				{
					SubmapIdx--;
					continue;
				}
			}
			else
			{
				SubmapIdx--;
				continue;
			}
		}
		else { SubmapIdx--; }
	}
	if (FindFlag)
	{
		LastLoopPose = initial_pose;
		LoopCount = 0;
	}
	else
	{
		NotFind = 1;
		distance_drive = 0;
		CurrentlyScanPose = initial_pose;
		LoopCount++;
	}
		//else
	//{
	//	// �´μ���ؼ�֡�ٽ���һ�μ����ʱ���ٽ��лػ�һ��
	//	LoopCount = 9;
	//}
}

bool LoopClosure::CheckLoopFlag() 
{ 
	LoopCount++;
	if (LoopCount == 10)
	{
		LoopCount = 0;
		return true;
	}
	return false;
}

void LoopClosure::Solve(map<int, Submap_Scan>& SubmapScanInfoPool_)
{
	// 建立问题
	ceres::Problem problem_;
	ceres::Problem::Options options_;
	// 建立数据存贮
	map<int, std::array<double, 3>> C_submaps;
	map<int, std::array<double, 3>> C_nodes;
	
	//给数据存储添加数据
	// submap
	for （auto & element_pair:SubmapScanInfoPool_)
	{
		
	}


}