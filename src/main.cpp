#include <iostream>
#include "slam_2d_v01.h"
#include "PGM.h"
#include "Test.h"
#include "CSM.h"
#include "M2MPGM.h"
using namespace std;
// CastRays ʱ�����Ĵ�

int main()
{ 

    // Initialization Part
    MapBuilder mapbuilder_;
    // �����״����
    mapbuilder_.SetLaserParameter();
    // ���������һЩ������ ����Ĳ����Ƚ���Ҫ��֮��ƥ���ʱ��������
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
    // ����
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // ��ȡ��֡����
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter ��������һ���˲� 
        MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, pixelSize);

        // �Ե�һ֡���г�ʼ��
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            // Ϊ��ͼ�����һ֡
            LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
            pgm_set.push_back(new PGMBuilder());
            pgm_set.back()->InsertFirstScan(temp, pixelSize);
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

        // TODO CERES ��ƥ��



        // ��������˶��ж� �Ƿ���ӹؼ�֡ ��Щif ����Ҳ����MotionFilter
        if (abs(new_pose(0) - motion_filter_pose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - motion_filter_pose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(2) - motion_filter_pose(2)) > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            motion_filter_pose = new_pose;
            // ������ٶ�ʧ�Ĳ�����CSM��ƥ��������
            if (abs(new_pose(0) - CurrentPose(0)) > 0.3 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.3 ||
                abs(new_pose(2) - CurrentPose(2)) > 0.2)
            {
                cout << "Tracking LOST!!!!" << endl;
                // TODO ��������ֱ����ת���ض�λ��ģ���������ض�λ
                break;
            }
            else
            {
                // ������Щ�������ᱻ������Mapbuilder��
                LaserData temp_node_submap = BasicFunction::Transform(voxel_filtered_scan, new_pose);
                LaserFan temp_insert(temp_node_submap, new_pose);
                // ��ͼ����Node 
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


        // ��ӵ�ǰλ�˵�·��
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);
    }

    cout << "MISSION COMPLETED" << endl;
    mapbuilder_.Drawer(mapbuilder_.AllPoints, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);
    return 0;
}

