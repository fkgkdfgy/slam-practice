#include "Test.h"
using namespace std;

void TEST_PGM_INITIALIZATION()
{
    PGMBuilder::Initialization();
    cout << PGMBuilder::InitializationFlag << endl;
    uint16_t count = 3000;
    uint16_t temp = 0;
    cout << " Print the Value to Probability" << endl;
    while (temp < 32768)
    {
        cout << PGMBuilder::ValueToProbabilityChart[temp] << "    ";
        temp += count;
    }
    cout << endl;
    temp = 0;
    cout << " Print the Value to CorrespondenceCost" << endl;
    while (temp < 32768)
    {
        cout << PGMBuilder::ValueToCorrespondenceCost[temp] << "    ";
        temp += count;
    }
    cout << endl;
    temp = 0;
    cout << " Print the HitAddChartCorrespondenceCost" << endl;
    while (temp < 32768)
    {
        cout << PGMBuilder::ValueToCorrespondenceCost[(PGMBuilder::HitAddChartCorrespondenceCost[temp]) - 32768] << "    ";
        temp += count;
    }
    cout << endl;
    temp = 0;
    cout << " Print the MIssAddChartCorrespondenceCost" << endl;
    while (temp < 32768)
    {
        cout << PGMBuilder::ValueToCorrespondenceCost[PGMBuilder::MissAddChartCorrespondenceCost[temp] - 32768] << "    ";
        temp += count;
    }
    cout << endl;
}
void TEST_PGM_INSERTSCAN() 
{
    // PGMBuilder Test  Part
// Initialization test part
    MapBuilder mapbuilder_;
    mapbuilder_.SetLaserParameter();
    mapbuilder_.SetMapParameter();
    mapbuilder_.GetAllData("E:\\OpenCV\\opencv_learning\\ConsoleApplication1\\b.txt", mapbuilder_.AllLaserData);
    MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(0, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
    MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);
    LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
    PGMBuilder::Initialization();
    cout << "if the initialization of PGMbuilder is completed " << "   ";
    cout << PGMBuilder::InitializationFlag << endl;
    PGMBuilder pgmbuilder_;
    pgmbuilder_.InsertFirstScan(temp, 0.2);
    pgmbuilder_.Finished();
    Vector2i a = pgmbuilder_.scale_.min() - Vector2i(3, 3);
    Vector2i b = pgmbuilder_.scale_.max() + Vector2i(-1, 20);
    vector<Array2d> temp_v;
    temp_v.push_back(Array2d(a(0) / 5, a(1) / 5));
    temp_v.push_back(Array2d(b(0) / 5, b(1) / 5));

    LaserFan temp_1;
    temp_1.end = temp_v;
    temp_1.pose = Array3d(-10, -10, 32);


    cout << endl << endl << endl << "Isert a new scan" << endl;
    pgmbuilder_.Insert(temp_1);
}
void TEST_CSM()
{
    CSM csm_;
    // TEST CSM 思路 CSM 负责粗匹配
    // Initialization Part
    MapBuilder mapbuilder_;
    // 设置雷达参数
    mapbuilder_.SetLaserParameter();
    // 设置另外的一些参数， 这里的参数比较重要在之后匹配的时候都能用上
    mapbuilder_.SetMapParameter();
    // 获取所有的的Scan 数据
    mapbuilder_.GetAllData("E:\\OpenCV\\opencv_learning\\ConsoleApplication1\\b.txt", mapbuilder_.AllLaserData);
    MapBuilder::Pose CurrentPose = MapBuilder::Pose(10, -2, 0);
    MapBuilder::LocalMap_ localMap;
    MapBuilder::Pose pose_guess;
    MapBuilder::Pose new_pose;
    // TEST_CSM
    CSM csmer_;
    vector<PGMBuilder> pgm_set;
    vector<PGMBuilder> pgmpool;
    // 主体

    mapbuilder_.ScanIdx = 1;
    // 读取单帧数据
    MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(0, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
    // Voxel Filter 对数据做一次滤波 
    MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

    // TEST_CSM
    PGMBuilder::Initialization();
    cout << "if the initialization of PGMbuilder is completed " << "   ";
    cout << PGMBuilder::InitializationFlag << endl;

    // 构建插入地图的第一帧
    LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
    pgm_set.push_back(PGMBuilder());

    // 插入地图第一帧

    pgm_set.back().InsertFirstScan(temp, double(0.05));
    csmer_.GenerateSearchWindow(temp, 0.4, 0.5, 0.05);


    // 准备构建地图的第二帧
    LaserData temp_2 = BasicFunction::Transform(voxel_filtered_scan, Array3d(10, -2, 0.4));
    LaserFan temp_insert(temp_2, Array3d(10, -2, 0.4));

    // 插入第二帧
    for (int i = 0; i < 6; i++)
    {
        pgm_set.back().Insert(temp);
    }
    pgm_set.back().Insert(temp_insert);

    // 构建一个相对第二帧的微小运动
    MapBuilder::LaserData temp_1 = BasicFunction::Transform(voxel_filtered_scan, -1 * Array3d(0, 0, 0));
    Array3d pose = csmer_.Match(temp_1, Array3d(0, 0, 0), pgm_set.back());
    pose = csmer_.Match(temp_1, Array3d(10.2, -2.3, 0.2), pgm_set.back());
    cout << "the result is " << pose(0) << "," << pose(1) << "," << pose(2) << endl;
}
// 构建的多级地图是可以使用的
void TEST_M2MPGM()
{
    // 测试地图生成的正确性

    CSM csm_;
    // TEST CSM 思路 CSM 负责粗匹配
    // Initialization Part
    MapBuilder mapbuilder_;
    // 设置雷达参数
    mapbuilder_.SetLaserParameter();
    // 设置另外的一些参数， 这里的参数比较重要在之后匹配的时候都能用上
    mapbuilder_.SetMapParameter();
    // 获取所有的的Scan 数据
    mapbuilder_.GetAllData("E:\\OpenCV\\opencv_learning\\ConsoleApplication1\\b.txt", mapbuilder_.AllLaserData);
    MapBuilder::Pose CurrentPose = MapBuilder::Pose(10, -2, 0);
    MapBuilder::LocalMap_ localMap;
    MapBuilder::Pose pose_guess;
    MapBuilder::Pose new_pose;
    // TEST_CSM
    CSM csmer_;
    vector<PGMBuilder> pgm_set;
    vector<PGMBuilder> pgmpool;
    // 主体

    mapbuilder_.ScanIdx = 1;
    // 读取单帧数据
    MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(0, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);

    // Voxel Filter 对数据做一次滤波 
    MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

    // TEST_CSM
    PGMBuilder::Initialization();
    cout << "if the initialization of PGMbuilder is completed " << "   ";
    cout << PGMBuilder::InitializationFlag << endl;

    // 构建插入地图的第一帧
    LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
    pgm_set.push_back(PGMBuilder());

    // 插入地图第一帧
    pgm_set.back().InsertFirstScan(temp, double(0.1));
    csmer_.GenerateSearchWindow(temp, 0.4, 0.15, 0.1);

    // 让地图的一个点变得明显
    Array2i point(0, 0);
    pgm_set.back().ProbabilityValue[BasicFunction::FlatToIndex(point, pgm_set.back().GetMapLimits())] = 2;
    cout << pgm_set.back().ProbabilityValue[15] << endl;
    int height = 3;
    M2MPGM m2mpgm_(make_shared<PGMBuilder>(pgm_set.back()), height);
    cout << "The m2mpgm_ creating is finished  : " << m2mpgm_.GetFinished() << endl;

    // show the different level
    Array2i start(0, 16);
    MapLimits temp_limits = m2mpgm_.GetMapLimits();
    for (int i = 0; i <= height; i++)
    {
        cout << "Show the " << i << " th" << endl;
        const vector<float>& temp = m2mpgm_.GetMap(i);
        for (int row = 15; row >= 0; row--)
        {
            for (int col = 0; col < 20; col++)
            {
                cout << round(temp[(row) * (m2mpgm_.GetMapLimits().scaleXY(0) + int(1 << i) - 1) + col] * 100) << "   ";

            }
            cout << endl;
        }
        cout << endl;
    }

}
void TEST_BBM_SIMPLE()
{
    // 第一个能不能找到正确的值
 // 第二个时间加速的问题
 // 第三个能不能对submap 做一个对比

    CSM csm_;
    // TEST CSM 思路 CSM 负责粗匹配
    // Initialization Part
    MapBuilder mapbuilder_;
    // 设置雷达参数
    mapbuilder_.SetLaserParameter();
    // 设置另外的一些参数， 这里的参数比较重要在之后匹配的时候都能用上
    mapbuilder_.SetMapParameter();
    // 获取所有的的Scan 数据
    mapbuilder_.GetAllData("E:\\OpenCV\\opencv_learning\\ConsoleApplication1\\b.txt", mapbuilder_.AllLaserData);
    MapBuilder::Pose CurrentPose = MapBuilder::Pose(10, -2, 0);
    MapBuilder::LocalMap_ localMap;
    MapBuilder::Pose pose_guess;
    MapBuilder::Pose new_pose;
    // TEST_CSM
    CSM csmer_;
    vector<PGMBuilder> pgm_set;
    vector<PGMBuilder> pgmpool;
    // 主体

    mapbuilder_.ScanIdx = 1;
    // 读取单帧数据
    MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(0, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);

    // Voxel Filter 对数据做一次滤波 
    MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

    // TEST_CSM
    PGMBuilder::Initialization();
    cout << "if the initialization of PGMbuilder is completed " << "   ";
    cout << PGMBuilder::InitializationFlag << endl;

    // 构建插入的内容
    vector<Array2d> temp_insert = BasicFunction::Transform(voxel_filtered_scan, Array3d(30, -21, 2.7));

    // 构建插入地图的第一帧
    LaserFan temp(temp_insert, Array3d(0, 0, 0));
    LaserFan temp_laser(voxel_filtered_scan, Array3d::Zero());
    pgm_set.push_back(PGMBuilder());

    // 插入地图第一帧
    pgm_set.back().InsertFirstScan(temp, double(0.1));
    csmer_.GenerateSearchWindow(temp_laser, 0.4, 0.15, 0.1);

    // 构建多层地图
    M2MPGM m2mpgm_(make_shared<PGMBuilder>(pgm_set.back()));
    //Array3d estimation = Array3d(10.2, - 7.3, 1.55);
    //Array3d pose = csmer_.Match(voxel_filtered_scan, estimation, pgm_set.back());
    //cout << pose(0) << "   " << pose(1) << "    " << pose(2) << endl;


    // 构建一个大运动
    // 初值赋值4.7时，顺利通过BBM 找到最优也就是说
    // 10 -7 的时候，答案也是错误的就需要看一下了
    CSM::MatchResult_BBM result = csmer_.BranchAndBoundMatch(make_shared<M2MPGM>(m2mpgm_), voxel_filtered_scan, Array3d(0, 0, 0));
    Array3d BBM_result = result.bestpose + Array3d(result.candidate(0), result.candidate(1), 0) * pgm_set.back().GetPixelSize();
    cout << result.bestscore << endl;
    cout << "the result pose is " << BBM_result(0) << "   " << BBM_result(1) << "  " << BBM_result(2) << endl;

}
void TEST_BBM_DIFF()
{}
void TEST_CERES()
{


}
void TEST_KDtree()
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
    // 主体
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // 读取单帧数据
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter 对数据做一次滤波 
        MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

        // 对第一帧进行初始化
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;
            // Kdtree part
            //vector<Vector2d> voxel_vector2d;
            //voxel_vector2d.reserve(voxel_filtered_scan.size());
            //for (auto& element : voxel_filtered_scan)
            //{
            //    voxel_vector2d.push_back(Vector2d(element));
            //}
            //mapbuilder_.kdtree_.Initialization(voxel_vector2d);
            continue;
        }

        // 从全局地图点中抽取 部分点来做局部地图用来匹配
        if (mapbuilder_.MapParameter_.miniUpdated == true)
        {
            //if (ScanIdx == 342)
            //    cout << "check here!" << endl;
            localMap = mapbuilder_.ExtractLocalMap_Kdtree(mapbuilder_.AllPoints, voxel_filtered_scan, CurrentPose, mapbuilder_.MapParameter_.BorderSize);
            // 把局部地图从点云地图变成栅格地图
            mapbuilder_.OccGridRough =
                mapbuilder_.ComputeOccGrid(localMap, mapbuilder_.MapParameter_.RoughpixelSize);
            // 产生更高分辨率的地图
            mapbuilder_.OccGridFiner =
                mapbuilder_.ComputeOccGrid(localMap, mapbuilder_.MapParameter_.RoughpixelSize / 2);
        }

        // 以匀速运动模型，对当前帧提供一个位姿的估计值
        pose_guess = mapbuilder_.RoughPoseEstimation(CurrentPose, ScanIdx);

        // 进行粗匹配
        MapBuilder::MatchResult result_temp;
        if (mapbuilder_.MapParameter_.miniUpdated == true)
        {
            result_temp =
                mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridRough, mapbuilder_.MapParameter_.fastResolution, voxel_filtered_scan, pose_guess, mapbuilder_.MapParameter_.RoughpixelSize);
        }
        else
        {
            result_temp =
                mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridFiner, mapbuilder_.MapParameter_.fastResolution, voxel_filtered_scan, pose_guess, mapbuilder_.MapParameter_.RoughpixelSize);
        }

        // 用高分辨率进行细匹配
        mapbuilder_.matchresult =
            mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridFiner, mapbuilder_.MapParameter_.fastResolution / 2, voxel_filtered_scan, result_temp.pose_matched, mapbuilder_.MapParameter_.RoughpixelSize / 2);


        // 根据相对运动判断 是否添加关键帧
        new_pose = mapbuilder_.matchresult.pose_matched;
        if (abs(new_pose(0) - CurrentPose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - CurrentPose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            0 > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            if (abs(new_pose(0) - CurrentPose(0)) > 0.5 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.5)
                cout << "BIGBIGBIG ESTIMATION!!!!!" << endl;
            else
            {
                mapbuilder_.MapParameter_.miniUpdated = true;
                mapbuilder_.AddKeyScan(mapbuilder_.matchresult);
            }
        }
        else
        {
            mapbuilder_.MapParameter_.miniUpdated = false;
        }
        // 添加当前位姿到路径
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);

    }
    mapbuilder_.Drawer(mapbuilder_.AllPoints, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);
}
void TEST_SLAM_V01()
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
    CSM csmer_;
    vector<PGMBuilder> pgm_set;
    vector<PGMBuilder> pgmpool;
    // 主体
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // 读取单帧数据
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter 对数据做一次滤波 
        MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

        // 对第一帧进行初始化
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            continue;
        }

        // 从全局地图点中抽取 部分点来做局部地图用来匹配
        if (mapbuilder_.MapParameter_.miniUpdated == true)
        {

            localMap = mapbuilder_.ExtractLocalMap(mapbuilder_.AllPoints, voxel_filtered_scan, CurrentPose, mapbuilder_.MapParameter_.BorderSize);
            // 把局部地图从点云地图变成栅格地图
            mapbuilder_.OccGridRough =
                mapbuilder_.ComputeOccGrid(localMap, mapbuilder_.MapParameter_.RoughpixelSize);
            // 产生更高分辨率的地图
            mapbuilder_.OccGridFiner =
                mapbuilder_.ComputeOccGrid(localMap, mapbuilder_.MapParameter_.RoughpixelSize / 2);
        }

        // 以匀速运动模型，对当前帧提供一个位姿的估计值
        pose_guess = mapbuilder_.RoughPoseEstimation(CurrentPose, ScanIdx);

        // 进行粗匹配
        MapBuilder::MatchResult result_temp;
        if (mapbuilder_.MapParameter_.miniUpdated == true)
        {
            result_temp =
                mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridRough, mapbuilder_.MapParameter_.fastResolution, voxel_filtered_scan, pose_guess, mapbuilder_.MapParameter_.RoughpixelSize);
        }
        else
        {
            result_temp =
                mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridFiner, mapbuilder_.MapParameter_.fastResolution, voxel_filtered_scan, pose_guess, mapbuilder_.MapParameter_.RoughpixelSize);
        }

        // 用高分辨率进行细匹配
        mapbuilder_.matchresult =
            mapbuilder_.HillClimbMatcher(mapbuilder_.OccGridFiner, mapbuilder_.MapParameter_.fastResolution / 2, voxel_filtered_scan, result_temp.pose_matched, mapbuilder_.MapParameter_.RoughpixelSize / 2);


        // 根据相对运动判断 是否添加关键帧
        new_pose = mapbuilder_.matchresult.pose_matched;
        if (abs(new_pose(0) - CurrentPose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - CurrentPose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            0 > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            if (abs(new_pose(0) - CurrentPose(0)) > 0.5 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.5)
                cout << "BIGBIGBIG ESTIMATION!!!!!" << endl;
            else
            {
                mapbuilder_.MapParameter_.miniUpdated = true;
                mapbuilder_.AddKeyScan(mapbuilder_.matchresult);
            }
        }
        else
        {
            mapbuilder_.MapParameter_.miniUpdated = false;
        }
        // 添加当前位姿到路径
        CurrentPose = new_pose;
        mapbuilder_.path.push_back(CurrentPose);

    }
    mapbuilder_.Drawer(mapbuilder_.AllPoints, 0.2, mapbuilder_.ScalarSet_);
    cv::waitKey(0);

}
void TEST_SLAM_V02()
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
    CSM csmer_;
    deque<PGMBuilder*> pgm_set;
    vector<PGMBuilder*> pgmpool;
    int submap_scale = 10;
    // 主体
    for (int ScanIdx = 0; ScanIdx < mapbuilder_.AllLaserData.size(); ScanIdx++)
    {
        cout << "the ScanIdx is " << ScanIdx << endl;
        mapbuilder_.ScanIdx = ScanIdx;
        // 读取单帧数据
        MapBuilder::LaserData CurrentScan = mapbuilder_.ReadAScan(ScanIdx, mapbuilder_.AllLaserData, mapbuilder_.MapParameter_);
        // Voxel Filter 对数据做一次滤波 
        MapBuilder::LaserData voxel_filtered_scan = mapbuilder_.VoxelFilter(CurrentScan, mapbuilder_.MapParameter_.FinerpixelSize);

        // 对第一帧进行初始化
        if (ScanIdx == 0)
        {
            mapbuilder_.Initialization(voxel_filtered_scan, CurrentPose);
            mapbuilder_.MapParameter_.miniUpdated = true;

            // 为地图插入第一帧
            LaserFan temp(voxel_filtered_scan, Array3d(0, 0, 0));
            pgm_set.push_back(new PGMBuilder());
            pgm_set.back()->InsertFirstScan(temp, double(0.1));
            csmer_.GenerateSearchWindow(temp,
                0.4 // CSM 平移搜素范围
                , 0.5, // CSM 旋转搜索范围
                mapbuilder_.MapParameter_.FinerpixelSize);

            continue;
        }

        // 以匀速运动模型，对当前帧提供一个位姿的估计值
        pose_guess = mapbuilder_.RoughPoseEstimation(CurrentPose, ScanIdx);

        // CSM匹配,粗匹配
        Array3d new_pose = csmer_.Match(voxel_filtered_scan, pose_guess, *pgm_set.front());

        // TODO CERES 精匹配

        // 根据相对运动判断 是否添加关键帧 这些if 函数也就是MotionFilter
        if (abs(new_pose(0) - CurrentPose(0)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            abs(new_pose(1) - CurrentPose(1)) > mapbuilder_.MapParameter_.miniUpdatedDT ||
            0 > mapbuilder_.MapParameter_.miniUpdateDR)
        {
            // 这里跟踪丢失的参数和CSM的匹配参数相等
            if (abs(new_pose(0) - CurrentPose(0)) > 0.4 ||
                abs(new_pose(1) - CurrentPose(1)) > 0.4 ||
                abs(new_pose(2) - CurrentPose(2)) > 0.5)
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
                }
                else
                {
                    pgmpool.push_back(pgm_set.front());
                    pgm_set.pop_front();
                    pgm_set.push_back(new PGMBuilder());
                    pgm_set.front()->Insert(temp_insert);
                    pgm_set.back()->InsertFirstScan(temp_insert, double(0.1));
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


}