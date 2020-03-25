#include <iostream>
#include "slam_2d_v01.h"
#include "PGM.h"
#include "Test.h"
#include "CSM.h"
#include "M2MPGM.h"
#include "HillClimbFiner.h"
#include "LoopClosure.h"
using namespace std;
typedef Submap_Scan::Node_info::status NodeStatus;
// CastRays 时间消耗大户
extern int count_0;
int main()
{

    // Initialization Part
    MapBuilder mapbuilder_;
    // 设置雷达参数
    mapbuilder_.SetLaserParameter();
    // 设置另外的一些参数， 这里的参数比较重要在之后匹配的时候都能用上
    mapbuilder_.SetMapParameter();
    // 获取所有的的Scan 数据
    mapbuilder_.GetAllData("E:\\OpenCV\\opencv_learning\\ConsoleApplication1\\b.txt", mapbuilder_.AllLaserData);
    MapBuilder::Pose CurrentPose = MapBuilder::Pose(0, 0, 0);
    MapBuilder::LocalMap_ localMap;
    MapBuilder::Pose pose_guess;
    MapBuilder::Pose new_pose;
    
    // preparation for CSM
    PGMBuilder::Initialization();
    CSM Roughcsmer_;
    CSM Finercsmer_;
    deque<PGMBuilder*> pgm_set;
    vector<vector<Array2d>*> keyscanpool;
    vector<Array3d> keyposepool;
    map<int, Submap_Scan> SubmapScanInfoPool;
    Array3d motion_filter_pose(0, 0, 0);

    int submap_scale = 10;
    double pixelSize = 0.1;
    vector<vector<Array2d>> point_pool;

    // 回环的准备工作
    LoopClosure loopclosure_;


    // 主体
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // 读取单帧数据
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter 对数据做一次滤波 
        MapBuilder::LaserData &  voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, pixelSize);

        // 对第一帧进行初始化
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            // 为地图插入第一帧
            LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
            pgm_set.push_back(new PGMBuilder());

            pgm_set.back()->InsertFirstScan(temp, pixelSize);
            // Submap 添加信息
            auto result_flag = SubmapScanInfoPool.insert(pair<int, Submap_Scan>(pgm_set.back()->GetSubmapIdx(),
                Submap_Scan(pgm_set.back()->GetSubmapIdx())));
            assert(result_flag.second == true);
            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].Insert(pgm_set.back()->GetSubmapIdx(), ScanIdx,NodeStatus(1) ,CurrentPose,-1.f);
            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.back()->GetScale().min()));
            Roughcsmer_.GenerateSearchWindow(temp,
                0.2 // CSM 平移搜素范围
                , 0.1, // CSM 旋转搜索范围
                pixelSize);
            continue;
        }

        // 以匀速运动模型，对当前帧提供一个位姿的估计值
        pose_guess = mapbuilder_.RoughPoseEstimation(CurrentPose, ScanIdx);

        // CSM匹配,粗匹配
        Array3d new_pose = Roughcsmer_.Match(voxel_filtered_scan, pose_guess, *pgm_set.front());

        // HillClimb 精匹配
        float matching_score = -1.f;
        HillClimbFiner Finermatcher_(*pgm_set.front(), Roughcsmer_);
        new_pose = Finermatcher_.FinerEstimaiton(voxel_filtered_scan, new_pose,matching_score);
        cout << "the pose from the FinerEstimation is " << new_pose(0) << "   " << new_pose(1) << "   " << new_pose(2) << endl;

        // 根据相对运动判断 是否添加关键帧 这些if 函数也就是MotionFilter
        if (abs(new_pose(0) - motion_filter_pose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - motion_filter_pose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(2) - motion_filter_pose(2)) > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            motion_filter_pose = new_pose;
            // 这里跟踪丢失的参数和CSM的匹配参数相等
            if (abs(new_pose(0) - CurrentPose(0)) > 0.5 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.5)
            {
                cout << "Tracking LOST!!!!" << endl;
                // TODO 这里后面会直接跳转到重定位的模块来进行重定位
                break;
            }
            else
            {
                // 以下这些代码最后会被集成在Mapbuilder里
                // 插入
                LaserData temp_node_submap = BasicFunction::Transform(voxel_filtered_scan, new_pose);
                LaserFan temp_insert(temp_node_submap, new_pose);
                
                // 子图插入Node 
                if (pgm_set.front()->laserfan_num < submap_scale)
                {
                    pgm_set.front()->Insert(temp_insert);
                    // new
                    SubmapScanInfoPool[pgm_set.front()->GetSubmapIdx()].Insert(pgm_set.front()->GetSubmapIdx(), ScanIdx, NodeStatus(1), new_pose,matching_score);
                    SubmapScanInfoPool[pgm_set.front()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.front()->GetScale().min()));
                    if (pgm_set.size() != 1)
                    {
                        pgm_set.back()->Insert(temp_insert);
                        //new
                        SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].Insert(pgm_set.back()->GetSubmapIdx(), ScanIdx, NodeStatus(1), new_pose, matching_score);
                        SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.back()->GetScale().min()));
                    }
                    else
                    {
                        if (pgm_set.front()->laserfan_num == (submap_scale / 2) + 1)
                        {
                            pgm_set.push_back(new PGMBuilder());

                            pgm_set.back()->InsertFirstScan(temp_insert, pixelSize);
                            // submap 数据和搜集
                            auto result_flag = SubmapScanInfoPool.insert(pair<int, Submap_Scan>(pgm_set.back()->GetSubmapIdx(), 
                                Submap_Scan(pgm_set.back()->GetSubmapIdx())));
                            assert(result_flag.second == true);
                            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].Insert(pgm_set.back()->GetSubmapIdx(), ScanIdx, NodeStatus(1), new_pose, matching_score);
                            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.back()->GetScale().min()));

                        }
                    }
                }
                else
                {
                    loopclosure_.Insert(pgm_set.front());
                    if (pgm_set.size() != 1)
                        pgm_set.pop_front();
                    pgm_set.push_back(new PGMBuilder());
                    pgm_set.front()->Insert(temp_insert);
                    SubmapScanInfoPool[pgm_set.front()->GetSubmapIdx()].Insert(pgm_set.front()->GetSubmapIdx(), ScanIdx, NodeStatus(1), new_pose, matching_score);
                    SubmapScanInfoPool[pgm_set.front()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.front()->GetScale().min()));
                    pgm_set.back()->InsertFirstScan(temp_insert, pixelSize);
                    //submap 数据收集
                    auto result_flag = SubmapScanInfoPool.insert(pair<int, Submap_Scan>(pgm_set.back()->GetSubmapIdx(),
                        Submap_Scan(pgm_set.back()->GetSubmapIdx())));
                    assert(result_flag.second == true);
                    SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].Insert(pgm_set.back()->GetSubmapIdx(), ScanIdx, NodeStatus(1), new_pose, matching_score);
                    SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.back()->GetScale().min()));
                }

                // voxel_filtered_scan 保存在MapBuilder
                // PGMBuilder M2MPGM 放在loopclosure


                // 写一个回环的逻辑
                // 明确一个PGMBuilder 的逻辑 也就是子图的坐标实际上直接是scale_.min() + 角度为0 
                // 然后子节点都要根据这个来算相对位姿 
                // 因为插入数值的时候 扫描点都已经被转到了世界坐标系下
                // 还是按照Cartographer的逻辑对这个代码进行构建 也就是一共有两种 约束 submap 内部节点的约束和 submap 外部节点的约束
                
                // 1. 给每个KeyScan Node 匹配每一个Submap 一开始都进行全匹配
                // 2. 给每个完成的submap中的Node 做匹配 这一步好像和第一步有重叠
                // 3. outside 的点到达一定数目之后就进行一次回环优化
                // 第一步

                // 设置回环的要求 
                if (1)
                {
                    loopclosure_.ComputeConstraintForNode(Roughcsmer_, voxel_filtered_scan, new_pose, SubmapScanInfoPool, ScanIdx);
                    bool Ready = loopclosure_.CheckSolveFlag();
                    if (Ready)
                        loopclosure_.Solve();
                    else
                        continue;
                }

            }
        }

        // TODO tracking lost and relocation;
        

        // 添加当前位姿到路径
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);
    }

    cout << "MISSION COMPLETED" << endl;
    //mapbuilder_.Drawer(point_pool, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);
    return 0;
}

