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
	// 建立数据存贮
	map<int, std::array<double, 3>> C_submaps;
	map<int, std::array<double, 3>> C_nodes;
	map<int, Eigen::Array3d> R_submaps;
	//给数据存储添加数据初始值
	for(auto & element_pair : SubmapScanInfoPool_)
	{
		// 这里的是一个像素坐标
	    Eigen::Vector2i temp_coordinate = PGMPool[element_pair.first]->GetScale().min();
		// 这里得到的是真实的世界坐标
		Eigen::Array3d submap_coordinate_insert = Array3d(temp_coordinate(0),temp_coordinate(1),0) * PGMPool[element_pair.first]->GetPixelSize();
		C_submaps.insert(pair<int,std::array<double,3>>(element_pair.first,{submap_coordinate_insert(0),submap_coordinate_insert(1),submap_coordinate_insert(2)}));
		R_submaps.insert(pair<int,Eigen::Array3d>(element_pair.first,submap_coordinate_insert));
		// 插入Node 初始值
		for (auto & element_node_pair : element_pair.second.Node_info_pool)
		{
			if(element_node_pair.second.status_ == 1)
			C_nodes.insert(pair<int,std::array<double,3>>(element_node_pair.first,{element_node_pair.second.world_pose(0),
																				   element_node_pair.second.world_pose(1),
																				   element_node_pair.second.world_pose(2)}));
		}
	}

	// 给问题添加变量
	for (auto &  element_submap:C_submaps)
	{
		problem_.AddParameterBlock(element_submap.second.data(),3);
		if (element_submap.first == 0)
		problem_.SetParameterBlockConstant(element_submap.second.data());
	}
	for (auto & element_node:C_nodes)
	{
		problem_.AddParameterBlock(element_node.second.data(),3);
	}

	// 给问题添加约束
	for(auto & element_submap:SubmapScanInfoPool_)
	{
		Eigen::Array3d submap_pose_ = R_submaps[element_submap.first];
		double sin_ = std::sin(submap_pose_(2));
		double cos_ = std::cos(submap_pose_(2));
		for (auto & element_node:element_submap.second.Node_info_pool)
		{
			// 计算相对位姿
			Eigen::Array3d node_pose_ = element_node.second.world_pose;
			Eigen::Array3d delta_pose_ = node_pose_-submap_pose_;
			Eigen::Array3d relative_pose_ (cos_*delta_pose_(0) + sin_*delta_pose_(1),
										    -sin_*delta_pose_(0) + cos_ * delta_pose_(1),
											delta_pose_(2));
			problem_.AddResidualBlock(CreateSPAErrorFucntion(relative_pose_),
			(element_node.second.status_ == 1)? (nullptr):(new ceres::HuberLoss(20))
			,C_submaps[element_submap.first].data(),C_nodes[element_node.first].data());
		}
	}
	
	ceres::Solver::Summary summary_;
	ceres::Solver::Options options_;
	// Global Optimization 整体优化的配置
	options_.use_nonmonotonic_steps = false;
	options_.num_threads = 6;
	options_.max_num_iterations = 30;
	
	ceres::Solve(options_,&problem_,&summary_);

	// 存储结果

	for (auto & element_node:C_nodes)
	{
		auto temp = element_node.second;
		node_pose_result.insert(pair<int,Eigen::Array3d>(element_node.first,Array3d(temp[0],temp[1],temp[2])));
	}
	// submap 的结果暂时用不到就先不写
}