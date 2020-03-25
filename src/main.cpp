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
// CastRays ʱ�����Ĵ�
extern int count_0;
int main()
{

    // Initialization Part
    MapBuilder mapbuilder_;
    // �����״����
    mapbuilder_.SetLaserParameter();
    // ���������һЩ������ ����Ĳ����Ƚ���Ҫ��֮��ƥ���ʱ��������
    mapbuilder_.SetMapParameter();
    // ��ȡ���еĵ�Scan ����
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

    // �ػ���׼������
    LoopClosure loopclosure_;


    // ����
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // ��ȡ��֡����
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter ��������һ���˲� 
        MapBuilder::LaserData &  voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, pixelSize);

        // �Ե�һ֡���г�ʼ��
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            // Ϊ��ͼ�����һ֡
            LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
            pgm_set.push_back(new PGMBuilder());

            pgm_set.back()->InsertFirstScan(temp, pixelSize);
            // Submap ������Ϣ
            auto result_flag = SubmapScanInfoPool.insert(pair<int, Submap_Scan>(pgm_set.back()->GetSubmapIdx(),
                Submap_Scan(pgm_set.back()->GetSubmapIdx())));
            assert(result_flag.second == true);
            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].Insert(pgm_set.back()->GetSubmapIdx(), ScanIdx,NodeStatus(1) ,CurrentPose,-1.f);
            SubmapScanInfoPool[pgm_set.back()->GetSubmapIdx()].UpdateScale_Min(Array2i(pgm_set.back()->GetScale().min()));
            Roughcsmer_.GenerateSearchWindow(temp,
                0.2 // CSM ƽ�����ط�Χ
                , 0.1, // CSM ��ת������Χ
                pixelSize);
            continue;
        }

        // �������˶�ģ�ͣ��Ե�ǰ֡�ṩһ��λ�˵Ĺ���ֵ
        pose_guess = mapbuilder_.RoughPoseEstimation(CurrentPose, ScanIdx);

        // CSMƥ��,��ƥ��
        Array3d new_pose = Roughcsmer_.Match(voxel_filtered_scan, pose_guess, *pgm_set.front());

        // HillClimb ��ƥ��
        float matching_score = -1.f;
        HillClimbFiner Finermatcher_(*pgm_set.front(), Roughcsmer_);
        new_pose = Finermatcher_.FinerEstimaiton(voxel_filtered_scan, new_pose,matching_score);
        cout << "the pose from the FinerEstimation is " << new_pose(0) << "   " << new_pose(1) << "   " << new_pose(2) << endl;

        // ��������˶��ж� �Ƿ����ӹؼ�֡ ��Щif ����Ҳ����MotionFilter
        if (abs(new_pose(0) - motion_filter_pose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - motion_filter_pose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(2) - motion_filter_pose(2)) > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            motion_filter_pose = new_pose;
            // ������ٶ�ʧ�Ĳ�����CSM��ƥ��������
            if (abs(new_pose(0) - CurrentPose(0)) > 0.5 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.5)
            {
                cout << "Tracking LOST!!!!" << endl;
                // TODO ��������ֱ����ת���ض�λ��ģ���������ض�λ
                break;
            }
            else
            {
                // ������Щ�������ᱻ������Mapbuilder��
                // ����
                LaserData temp_node_submap = BasicFunction::Transform(voxel_filtered_scan, new_pose);
                LaserFan temp_insert(temp_node_submap, new_pose);
                
                // ��ͼ����Node 
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
                            // submap ���ݺ��Ѽ�
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
                    //submap �����ռ�
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
        

        // ���ӵ�ǰλ�˵�·��
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);
    }

    cout << "MISSION COMPLETED" << endl;
    //mapbuilder_.Drawer(point_pool, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);
    return 0;
}

