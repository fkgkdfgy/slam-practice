#include "slam_2d_v01.h"
// #include "common.h"
using namespace std;
using namespace Eigen;
using namespace cv;
// the function running 
void MapBuilder::Start()
{

    SetLaserParameter();
    SetMapParameter();
    GetAllData("b.txt", AllLaserData);
}

// main funcitons
// 设置雷达参数
void MapBuilder::SetLaserParameter()
{
    LaserParameter_.angle_min = -2.351831;
    LaserParameter_.angle_max = 2.351831;
    LaserParameter_.angle_increment = 0.004363;
    LaserParameter_.npoints = 1079;
    LaserParameter_.range_min = 0.023;
    LaserParameter_.range_max = 60;
    LaserParameter_.time_increment = 1.7361 * 0.00001;
}
// 设置地图参数
void MapBuilder::SetMapParameter()
{
    MapParameter_.BorderSize = 1;
    MapParameter_.RoughpixelSize = 0.2;
    MapParameter_.FinerpixelSize = 0.2 / 2;
    MapParameter_.miniUpdated = false;
    MapParameter_.miniUpdatedDT = 0.01;
    MapParameter_.miniUpdateDR = 0.02;
    MapParameter_.range_eff_max = 24;
    MapParameter_.range_eff_min = 0.5;
    MapParameter_.fastResolution = Array3d(0.05, 0.05, 0.008726646259972);
}
// 读取所有帧的数据 用正则表达式实现 
bool MapBuilder::GetAllData(const string str_temp, map<unsigned int, MapBuilder::LaserData>& AllLaserData)
{
    const string str_ = str_temp;
    fstream file_;
    file_.open(str_);
    if (!file_.is_open())
    {
        cout << "there is something wrong with the file opening " << endl;
        return false;
    }
    string str_temp_;
    string temp;
    regex pattern("(\\d+)\\.(\\d+)");
    string::const_iterator iterStart;
    string::const_iterator iterEnd;
    unsigned int laser_count = 0;

    while (getline(file_, str_temp_))
    {
        iterStart = str_temp_.begin();
        iterEnd = str_temp_.end();
        smatch result;
        MapBuilder::LaserData temp_laser;
        temp_laser.reserve(LaserParameter_.npoints);
        double angle_increment = LaserParameter_.angle_increment;
        double angle_temp = LaserParameter_.angle_min;
        while (regex_search(iterStart, iterEnd, result, pattern))
        {
            temp = result[0];
            temp_laser.push_back(Array2d(stod(temp), angle_temp));
            // cout << temp << " ";
            iterStart = result[0].second;
            angle_temp = angle_temp + angle_increment;
        }
        AllLaserData.insert(pair<unsigned int, MapBuilder::LaserData>(laser_count, temp_laser));
        // TEST Part
        cout << AllLaserData.size() << endl;
        laser_count++;
        if (laser_count == 500)
            break;
    }
    return true;

}
// 读取单个帧的数据
MapBuilder::LaserData MapBuilder::ReadAScan(unsigned int ScanIdx, map<unsigned int, MapBuilder::LaserData>& AllLaserData, MapBuilder::MapParameter MapParameter_)
{

    MapBuilder::LaserData CurrentScan;

    // Get one Scan from the whole information set.
    // Moreover, transfer the polar coordinate to Cartesian coordinate system
    for (auto& element : AllLaserData[ScanIdx])
    {
        if (element(0) < MapParameter_.range_eff_max)
        {
            double ct = std::cos(element(1));
            double st = std::sin(element(1));
            CurrentScan.push_back(Eigen::Array2d(element(0) * ct, element(0) * st));
            // TEST PART
         //      cout<<CurrentScan.back()<<endl;
        }
    }

    return CurrentScan;
}

bool MapBuilder::Initialization(const MapBuilder::LaserData& CurrentScan, const MapBuilder::Pose CurrentPose)
{
    if (keyscans_.size() != 0)
    {
        std::cout << "there is something wrong with this calling " << endl;
        return false;
    }
    MapBuilder::KeyScan keyscan_temp;
    keyscan_temp.keyindex = 1;
    keyscan_temp.iBegin = 0;
    keyscan_temp.iEnd = CurrentScan.size() - 1;
    AllPoints.insert(AllPoints.end(), CurrentScan.begin(), CurrentScan.end());
    PointsPool.push_back(CurrentScan);
    keyscan_temp.pose_ = CurrentPose;
    keyscans_.push_back(keyscan_temp);
    return true;
}

// 体素滤波过滤点云
MapBuilder::LaserData MapBuilder::VoxelFilter(const MapBuilder::LaserData& Laser_data_, double pixelSize)
{
    using keyType = std::bitset<64 * 2>;
    unordered_set<keyType> voxel_set_;
    MapBuilder::LaserData filtered_pointCloud;

    int count = 1;
    for (auto& element : Laser_data_)
    {
        Array2d temp(element / pixelSize);
        
        int x_test = std::round(temp(0));
        int y_test = std::round(temp(1));
        
        keyType x = static_cast<uint64>(std::round(temp(0)));
        keyType y = static_cast<uint64>(std::round(temp(1)));
        //temp_file << x_test << "   " << y_test << endl;
        keyType temp_bitset((x << 64) | y);
        auto it_inserted = voxel_set_.insert(temp_bitset);
        // 如果曾经插入过相同的值 if 条件不会通过
        if (it_inserted.second)
        {
            filtered_pointCloud.push_back(element);
        }
        count++;

    }
    return filtered_pointCloud;
}

// 抽取局地图
MapBuilder::LocalMap_ MapBuilder::ExtractLocalMap(const MapBuilder::LaserData& Allpoints, const MapBuilder::LaserData& filtered_pointCloud
                                                , MapBuilder::Pose pose_, unsigned int BorderSize)
{

    // Here, the uint for BorderSize is Meter
    MapBuilder::LaserData localMap;
    MapBuilder::LaserData transformed_pointCloud = Transform(filtered_pointCloud, pose_);
    std::vector<Array2d>::iterator iterCurrent = transformed_pointCloud.begin();
    std::vector<Array2d>::iterator iterEnd = transformed_pointCloud.end();

    // find the mininum XY and maxium XY of the LaserData Cartesian Coordinates
    Eigen::Array2d minXY(*iterCurrent);
    Eigen::Array2d maxXY(*iterCurrent);
    for (; iterCurrent != iterEnd; iterCurrent++)
    {
        minXY = minXY.min(*iterCurrent - BorderSize);
        maxXY = maxXY.max(*iterCurrent + BorderSize);
    }

    // Extract LocalMap to registration.
    // BorderSize is limited in the range between minXY and maxXY
    for (auto& temp : Allpoints)
    {
        if (temp(0) > minXY(0) && temp(1) > minXY(1) && temp(0) < maxXY(0) && temp(1) < maxXY(1))
        {
            localMap.push_back(temp);
        }
    }
    return MapBuilder::LocalMap_(localMap, minXY, maxXY);
}

MapBuilder::LocalMap_ MapBuilder::ExtractLocalMap_Kdtree(const MapBuilder::LaserData& Allpoints, const MapBuilder::LaserData& filtered_pointCloud
                                                       , MapBuilder::Pose pose_, unsigned int BorderSize)
{
    // 寻找质心
    MapBuilder::LaserData localMap;
    MapBuilder::LaserData transformed_pointCloud = Transform(filtered_pointCloud, pose_);
    Array2d mean_xy = Array2d::Zero();
    Eigen::Array2d minXY(transformed_pointCloud[0] - BorderSize);
    Eigen::Array2d maxXY(transformed_pointCloud[0] + BorderSize);
    for (auto& element : transformed_pointCloud)
    {
        mean_xy += element;
        minXY = minXY.min(element - BorderSize);
        maxXY = maxXY.max(element + BorderSize);
    }
    mean_xy = mean_xy / filtered_pointCloud.size();
    kdtree_.k = 2000>kdtree_.ptr_pool.size()?kdtree_.ptr_pool.size():2000;
    
    vector<Vector2d> temp_result =  tree::findKNN(kdtree_, Vector2d(mean_xy));
    
    //vector<Vector2d> final_result = BasicFunction::VoxelFilter<Vector2d>(temp_result, 0.1);
    for (auto& temp : temp_result)
    {
        if (temp(0) > minXY(0) && temp(1) > minXY(1) && temp(0) < maxXY(0) && temp(1) < maxXY(1))
        {
            localMap.push_back(temp);
        }
    }
    return MapBuilder::LocalMap_(localMap, minXY, maxXY);
}
// local map just mean a smaller map .but is still represented by the world frame
MapBuilder::OccGrid MapBuilder::ComputeOccGrid(const MapBuilder::LocalMap_& local_map_, double pixelSize)
{
    MapBuilder::LaserData temp_point_cloud = local_map_.localmap_points;

    Array2d minXY(temp_point_cloud[0]);
    Array2d maxXY(temp_point_cloud[0]);
    for (auto& element : temp_point_cloud)
    {
        minXY = minXY.min(element);
        maxXY = maxXY.max(element);
    }
    //occMap frame
    Array2d occMap_minXY = minXY - 3 * pixelSize;
    Array2d occMap_maxXY = maxXY + 3 * pixelSize;
    // the world frame parameteres represented in the local frame 
    MapBuilder::Pose pose_temp(-occMap_minXY(0), -occMap_minXY(1), 0);

    // transfer the laserdate representation from world frame to localmap frame,
    // the all points xy >=0
    temp_point_cloud = Transform(temp_point_cloud, pose_temp);
    
    // the distance transform to generate the registration map 
    vector<Array2i> discretePointCloud = DiscretePoints(temp_point_cloud, pixelSize);
    // find the maximum XY point in the localmap frame
    Array2d scale_undicre((occMap_maxXY - occMap_minXY) / pixelSize);
    Array2i scale_discre(round(scale_undicre(0)), round(scale_undicre(1)));

    // use this scale to create a map of suitable scale
    // ie maxminum point is 7.56 7.42  => (discretation) => 8,7 because we count from 0,0 thus the row and cols of the map should be 9 and 8
    cv::Mat occgrid_map( scale_discre(1) + 1 , scale_discre(0)+1, CV_8UC1, cv::Scalar(1));
    // 

    for (auto& element : discretePointCloud)
    {
        Array2i temp_check = element;
        occgrid_map.at<unsigned char>(temp_check(1), temp_check(0)) = 0;

    }

    Mat result;
    Mat result_middle_2;
    distanceTransform(occgrid_map, result, DIST_L2, 5);
    Mat result_middle_1(result.rows, result.cols, CV_32F, Scalar(10));
    threshold(result, result, 10, 10, CV_THRESH_TRUNC);

    //cout << result << endl;
    return MapBuilder::OccGrid(result, occMap_minXY, scale_discre, pixelSize);
}

MapBuilder::Pose MapBuilder::RoughPoseEstimation(Pose& CurrentPose, const unsigned int& ScanIdx) const
{
    if (ScanIdx > 1)
        return DiffPose(path, CurrentPose) + CurrentPose;
    else
        return CurrentPose;
}

MapBuilder::MatchResult MapBuilder::HillClimbMatcher(const MapBuilder::OccGrid& OccGridMap,
    const Eigen::Array3d& MatcherResolution, MapBuilder::LaserData& ScanData, MapBuilder::Pose pose_guess,
    double pixelSize)
{


 // transfer the scan points from the scan frame into world frame
    double ipixel = 1 / OccGridMap.pixelSize;
    Array2d minXY = OccGridMap.minXY;
    cv::Mat metricMap = OccGridMap.occgrid;
    int rows = metricMap.rows;
    int cols = metricMap.cols;

    // the main part of HillClimb
    double trans_ = MatcherResolution(0);
    double rot_angle = MatcherResolution(2);
    Array3d trans_range(-1 * trans_, 0, trans_);
    Array3d rot_range(-1 * rot_angle, 0, rot_angle);
    unsigned int IterMax = 50;
    unsigned int depthMax = 3;
    unsigned int Iter = 0;
    unsigned int depth = 0;
    float bestscore = 1000000;
    std:vector<float> pointsValueIdex;
    std::vector<float> ValueIdexForBest;
    bool noChange = 1;

    MapBuilder::Pose dBestPose(0, 0, 0);
    MapBuilder::Pose middlePose = pose_guess;
    vector<int> BestScoreIndex;
    for (; Iter < IterMax; Iter++)
    {
        noChange = true;
        for (unsigned int angle_ = 0; angle_ < 3; angle_++)
        {
            for (unsigned int trans_x = 0; trans_x < 3; trans_x++)
            {
                for (unsigned int trans_y = 0; trans_y < 3; trans_y++)
                {
                    float score = 0;
                    // 变换的部分
                    Array3d a = MapBuilder::Pose(trans_range(trans_x), trans_range(trans_y), rot_range(angle_))+middlePose;
                    MapBuilder::Pose temp_pose = a + Array3d(-minXY(0),-minXY(1),0) ;
                    MapBuilder::LaserData scan_temp = Transform(ScanData, temp_pose);
                    vector<float> pointsValueIdex;
                    vector<int> scoreIdex;
                    pointsValueIdex.clear();
                    scoreIdex.clear();
                    // Score Part
                    int i = 0;
                    for (auto& element : scan_temp)
                    {
                        int temp_y = round(element(1) * ipixel);
                        int temp_x = round(element(0) * ipixel);
                        if (temp_x > 0 && temp_y > 0 && temp_y < rows && temp_x < cols)
                        {
                            float score_temp = metricMap.at<float>(temp_y, temp_x);
                            score += score_temp;
                            pointsValueIdex.push_back(score_temp);
                            scoreIdex.push_back(i);
                        }
                        i++;
                    }
                    if (score < bestscore)
                    {
                        bestscore = score;
                        dBestPose = a;
                        ValueIdexForBest = pointsValueIdex;
                        BestScoreIndex = scoreIdex;
                        noChange = 0;
                    }
                }
            }
        }
        if (noChange == 1)
        {
            depth++;
            trans_range = trans_range / 2;
            rot_range = rot_range / 2;

            if (depth > depthMax)
            {
                break;
            }

        }
        middlePose = dBestPose;
    }
    if (pixelSize == 0.1)
    {
        cout << "the oorginal pose is " << pose_guess(0) << "," << pose_guess(1) << endl;
        cout << Iter << endl;
    }

    return MapBuilder::MatchResult(dBestPose, ValueIdexForBest, ScanData, BestScoreIndex);
}
// 加入关键帧，主要目的是把 新的点加入到 点集里 因为还没有写回环 作用有限
void MapBuilder::AddKeyScan(const MapBuilder::MatchResult& matchresult_)
{
    MapBuilder::KeyScan keyscan_temp;
    keyscan_temp.pose_ = matchresult_.pose_matched;
    int size = matchresult_.pointsValueIdex_.size();
    int size_points = AllPoints.size();
    MapBuilder::LaserData laser_temp_add;
    MapBuilder::LaserData Laser_temp_world = Transform(matchresult_.ScanPoints_, matchresult_.pose_matched);
    // find the new points
    for (int i = 0; i < size; i++)
    {
        // 这里的if 就是用来检测是不是一个新点的
        if (matchresult_.pointsValueIdex_[i] > 1.1)
        {
            Eigen::Array2d temp_point = Laser_temp_world[i];
            laser_temp_add.push_back(temp_point);
        }
    }
    // KDtree Part
    //for (auto& element : laser_temp_add)
    //{
    //    kdtree_.InsertNode(Vector2d(element));
    //}
    PointsPool.push_back(laser_temp_add);
    AllPoints.insert(AllPoints.end(), laser_temp_add.begin(), laser_temp_add.end());
    keyscan_temp.keyindex = keyscans_.size() + 1;
    keyscan_temp.iBegin = size_points;
    keyscan_temp.iEnd = size_points + laser_temp_add.size() - 1;
    keyscans_.push_back(keyscan_temp);
}

// AddKeyScan for V02  为了模拟Cartographer的状态来进行submap和node的插入
void MapBuilder::AddKeyScan_v02(MapBuilder::MatchResult_v02 & matching_result) 
{



}




// Allpoints should have included the points in Scan_Data
void MapBuilder::Drawer(const MapBuilder::LaserData& Allpoints, double pixelSize, const MapBuilder::ScalarSet scalar_set_)
{

    Eigen::Array2d minXY(Allpoints[0]);
    Eigen::Array2d maxXY(Allpoints[0]);
    for (auto& element : Allpoints)
    {
        minXY = minXY.min(element);
        maxXY = maxXY.max(element);
    }
    Array2i scaleXY(round((maxXY(0) - minXY(0)) / pixelSize), round((maxXY(1) - minXY(1)) / pixelSize));
    cv::Mat map(scaleXY(1) + 50, scaleXY(0) + 50, CV_8UC3, cv::Scalar(255, 255, 255));
    imshow("TEST TEST", map);
    DrawPoint(map, AllPoints, pixelSize, scalar_set_.WholeMapColor, minXY -5);
    imshow("TEST TEST", map);
}

// Basic funtion supported
MapBuilder::LaserData MapBuilder::Transform(const LaserData& UnprocessScanData, Pose pose_)
{
    MapBuilder::LaserData transformed_scandata;
    double theta = pose_(2);
    double x = pose_(0);
    double y = pose_(1);
    double ct = std::cos(theta);
    double st = std::sin(theta);
    for (auto& temp : UnprocessScanData)
    {
        transformed_scandata.push_back(Eigen::Array2d(temp(0) * ct - temp(1) * st + x, temp(0) * st + temp(1) * ct + y));
    }
    return transformed_scandata;
}
MapBuilder::Pose MapBuilder::DiffPose(const vector<Pose>& path, Pose Current) const
{
    Array3d temp = Current - path.back();
    temp(2) = 0;
    return temp;
}
vector<Array2i> MapBuilder::DiscretePoints(const LaserData& PointCloud, double pixelSize) const
{
    vector<Eigen::Array2i> discretePointCloud;
    for (auto& point : PointCloud)
    {
        Array2d temp = point / pixelSize;
        int x = std::round(temp(0));
        int y = std::round(temp(1));
        discretePointCloud.push_back(Array2i(x, y));
    }
    return discretePointCloud;
}
// Drawer Basic supported function
// Draw points
void MapBuilder::DrawPoint(cv::Mat& map_, const LaserData& laserdata_, double pixelSize, const cv::Scalar& scalar_, const Array2d & _offsetXY)
{
    LaserData temp_laser = Transform(laserdata_,Array3d(-1*_offsetXY(0),-1*_offsetXY(1),0));
    vector<Array2i> discretePointCloud = DiscretePoints(temp_laser, pixelSize);
    for (auto& element : discretePointCloud)
    {

        Array2i temp = element;
        map_.at<cv::Vec3b>(temp(1), temp(0))[0] = scalar_(0);
        map_.at<cv::Vec3b>(temp(1), temp(0))[1] = scalar_(1);
        map_.at<cv::Vec3b>(temp(1), temp(0))[2] = scalar_(2);
    }
    //cout << "the map test is " << map_ << endl;
}
void MapBuilder::DrawLine(cv::Mat& map_, Eigen::Array2i start_, Eigen::Array2i end_, const cv::Scalar& scalar_)
{
    // DarkOrange for Line
    cv::Point start(start_(0), start_(1));
    cv::Point end(end_(0), end_(1));
    cv::line(map_, start, end, scalar_, 1);
}
void MapBuilder::DrawLine(cv::Mat& map_, cv::Point start_, cv::Point end_, const cv::Scalar& scalar_)
{
    cv::line(map_, start_, end_, scalar_, 1);
}
void MapBuilder::DrawPath(cv::Mat& map_, const  vector<MapBuilder::Pose>& path, double pixelSize, const cv::Scalar& scalar_)
{
    // transfor Pose type into Point
    vector<Array2d> path_arrary2d;
    for (auto& element : path)
    {
        path_arrary2d.push_back(Array2d(element(0), element(1)));
    }
    vector<Array2i> path_array2i = DiscretePoints(path_arrary2d, pixelSize);
    vector<Array2i>::iterator iBegin = path_array2i.begin();
    vector<Array2i>::iterator iEnd = path_array2i.end();
    iEnd--;
    for (; iBegin < iEnd; iBegin++)
    {
        DrawLine(map_, *iBegin, *(iBegin++), scalar_);
    }
}
void MapBuilder::DrawCar(cv::Mat& map_, const MapBuilder::Pose& pose_, const cv::Scalar& scalar_, double pixelSize)
{
    // length 50 cm width 30cm with blue color
    double length = 0.5;
    double width = 0.3;
    double directionLineLength = 4;
    vector<Array2d> CarPoint;
    CarPoint.push_back(Array2d(0, 0));
    CarPoint.push_back(Array2d(width, length));
    CarPoint.push_back(Array2d(0, length));
    CarPoint.push_back(Array2d(width, 0));
    CarPoint.push_back(Array2d(width / 2, length / 2));
    CarPoint.push_back(Array2d(width / 2, length / 2 + directionLineLength));
    vector<Array2d> transformed_carpoints = Transform(CarPoint, Pose(-width / 2 + pose_(0), -length / 2 + pose_(1), pose_(2)));
    vector<Array2i> discre_trans_carpoints = DiscretePoints(transformed_carpoints, pixelSize);
    DrawLine(map_, discre_trans_carpoints[0], discre_trans_carpoints[2], scalar_);
    DrawLine(map_, discre_trans_carpoints[0], discre_trans_carpoints[3], scalar_);
    DrawLine(map_, discre_trans_carpoints[1], discre_trans_carpoints[2], scalar_);
    DrawLine(map_, discre_trans_carpoints[1], discre_trans_carpoints[3], scalar_);
    DrawLine(map_, discre_trans_carpoints[4], discre_trans_carpoints[5], scalar_);
}
Mat MapBuilder::MatsAndOperator(Mat& mat1, Mat& mat2)
{
    assert((mat1.rows == mat2.rows) && (mat1.cols == mat2.cols));
    assert(mat1.type() == CV_32FC1);
    assert(mat2.type() == CV_16SC1);

    Mat result(mat1.rows, mat1.cols, CV_32FC1, Scalar(0));
    Mat_<float>::iterator iBegin_mat1 = mat1.begin<float>();
    Mat_<float>::iterator iEnd_mat1 = mat1.end<float>();
    Mat_<short>::iterator iBegin_mat2 = mat2.begin<short>();
    Mat_<float>::iterator iBegin_result = result.begin<float>();
    for (; iBegin_mat1 != iEnd_mat1; iBegin_mat1++)
    {
        if (*iBegin_mat2 == 1)
        {
            *iBegin_result = *iBegin_mat1;
        }
        iBegin_result++;
        iBegin_mat2++;
    }
    return result;
}
Mat MapBuilder::self_distanceTransform(Mat& pic_, int size_)
{
    assert((size_ % 2) == 1);
    assert(pic_.type() == CV_16SC1);
    //define a mask for DT Complete
    Mat mask(size_, size_, CV_32F);
    int center_x = ceil(size_ / 2);
    int center_y = center_x;
    int border_x = size_ - center_x - 1;
    int border_y = border_x;

    // give values to the mask
    for (int rows = 0; rows < size_; rows++)
    {
        float* data = mask.ptr<float>(rows);
        for (int cols = 0; cols < size_; cols++)
        {
            data[cols] = calculateL2distance(center_x, center_y, cols, rows);
        }
    }

    // Main Part
    Mat DTmap(pic_.rows + 2 * border_y, pic_.cols + 2 * border_x, CV_16SC1, Scalar(0));
    Mat DT(pic_.rows, pic_.cols, CV_32FC1, Scalar(0));
    Mat result(size_, size_, CV_32FC1, Scalar(0));
    Mat main_DTmap = DTmap(Rect(2, 2, pic_.cols, pic_.rows));
    pic_.copyTo(main_DTmap);
    Mat_<float>::iterator iBegin;
    Mat_<float>::iterator iEnd;
    Mat_<short>::iterator iBegin_temp;
    Mat_<short>::iterator iEnd_temp;
    Mat_<float>::iterator iBegin_DT = DT.begin<float>();
    float min = 10000;
    for (int i = 0; i < pic_.rows; i++)
    {
        for (int j = 0; j < pic_.cols; j++)
        {
            min = 10000;
            Mat temp = DTmap(Rect(j, i, size_, size_));
            result = MatsAndOperator(mask, temp);
            iBegin = result.begin<float>();
            iEnd = result.end<float>();
            iBegin_temp = temp.begin<short>();
            iEnd_temp = temp.end<short>();
            for (; iBegin != iEnd; iBegin++)
            {
                short temp_1 = *iBegin_temp;
                if (*iBegin_temp == 1 && *iBegin < min)
                {
                    min = *iBegin;
                }
                iBegin_temp++;
            }
            if (min == 10000)
                *iBegin_DT = 10;
            else
                *iBegin_DT = min;
            iBegin_DT++;
        }
    }
    return DT;
}
double MapBuilder::calculateL2distance(int x1, int y1, int x2, int y2)
{
    return double(sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)));
}





