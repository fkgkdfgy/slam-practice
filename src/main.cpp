#include <iostream>
#include "slam_2d_v01.h"
#include "PGM.h"
#include "Test.h"
#include "CSM.h"
#include "M2MPGM.h"
using namespace std;
// CastRays 时间消耗大户

int main()
{ 

    // Initialization Part
    MapBuilder mapbuilder_;
    // 设置雷达参数
    mapbuilder_.SetLaserParameter();
    // 设置另外的一些参数， 这里的参数比较重要在之后匹配的时候都能用上
    mapbuilder_.SetMapParameter();
    //
    mapbuilder_.GetAllData("/home/liu/CLionProjects/LidarSLAM2d/b.txt", mapbuilder_.AllLaserData);
    MapBuilder::Pose CurrentPose = MapBuilder::Pose(0, 0, 0);
    MapBuilder::LocalMap_ localMap;
    MapBuilder::Pose pose_guess;
    MapBuilder::Pose new_pose;

    // preparation for CSM
    PGMBuilder::Initialization();
    CSM Roughcsmer_;
    CSM Finercsmer_;
    deque<PGMBuilder*> pgm_set;
    vector<PGMBuilder*> pgmpool;
    Array3d motion_filter_pose(0, 0, 0);

    int submap_scale = 10;
    double pixelSize = 0.05;
    // 主体
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // 读取单帧数据
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter 对数据做一次滤波 
        MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, pixelSize);

        // 对第一帧进行初始化
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            // 为地图插入第一帧
            LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
            pgm_set.push_back(new PGMBuilder());
            pgm_set.back()->InsertFirstScan(temp, pixelSize);
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

        // TODO CERES 精匹配



        // 根据相对运动判断 是否添加关键帧 这些if 函数也就是MotionFilter
        if (abs(new_pose(0) - motion_filter_pose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - motion_filter_pose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(2) - motion_filter_pose(2)) > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            motion_filter_pose = new_pose;
            // 这里跟踪丢失的参数和CSM的匹配参数相等
            if (abs(new_pose(0) - CurrentPose(0)) > 0.3 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.3 ||
                abs(new_pose(2) - CurrentPose(2)) > 0.2)
            {
                cout << "Tracking LOST!!!!" << endl;
                // TODO 这里后面会直接跳转到重定位的模块来进行重定位
                break;
            }
            else
            {
                // 以下这些代码最后会被集成在Mapbuilder里
                LaserData temp_node_submap = BasicFunction::Transform(voxel_filtered_scan, new_pose);
                LaserFan temp_insert(temp_node_submap, new_pose);
                // 子图插入Node 
                if (pgm_set.front()->laserfan_num < submap_scale)
                {
                    pgm_set.front()->Insert(temp_insert);
                    if (pgm_set.size() != 1)
                    {
                        pgm_set.back()->Insert(temp_insert);
                    }
                    else
                    {
                        if (pgm_set.front()->laserfan_num == (submap_scale / 2)+1)
                        {
                            pgm_set.push_back(new PGMBuilder());
                            pgm_set.back()->InsertFirstScan(temp_insert, pixelSize);
                        }
                    }
                }
                else
                {
                    pgmpool.push_back(pgm_set.front());
                    if(pgm_set.size()!=1)
                    pgm_set.pop_front();
                    pgm_set.push_back(new PGMBuilder());
                    pgm_set.front()->Insert(temp_insert);
                    pgm_set.back()->InsertFirstScan(temp_insert,pixelSize );
                }
            }
        }

        // TODO tracking lost and relocation;


        // 添加当前位姿到路径
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);
    }

    cout << "MISSION COMPLETED" << endl;
    mapbuilder_.Drawer(mapbuilder_.AllPoints, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);
    return 0;
}

