#ifndef _SLAM_2D_V01_H_
#define _SLAM_2D_V01_H_
#include "common.h"


using namespace Eigen;
using namespace std;
using namespace cv;
class MapBuilder
{
public:
    typedef vector<Array2d> LaserData;
    // pose(0) ->x pose(1) ->y pose(2)->theta
    typedef Eigen::Array3d Pose;
    
    struct LaserParameters
    {
        double angle_min;
        double angle_max;
        double angle_increment;
        unsigned int npoints;
        double range_min;
        double range_max;
        double scan_time;
        double time_increment;
        double angles;
    };
    
    struct KeyScan
    {
        MapBuilder::Pose pose_;
        unsigned int keyindex;
        unsigned int iBegin;
        unsigned int iEnd;
    };
        
    struct OccGrid
    {
        OccGrid ( const cv::Mat & occgrid_, const Array2d & minXY_ , const Array2i & scaleXY_,const double & pixelSize_)
        :occgrid(occgrid_),minXY(minXY_),scaleXY(scaleXY_),pixelSize(pixelSize_){}
        OccGrid() {}
        cv::Mat occgrid;
        Array2d minXY;
        Array2i scaleXY;
        double pixelSize;
    };
    
    struct MatchResult
    {
        MatchResult(const MapBuilder::Pose & pose,const std::vector<float> & pointsValueIdex,const MapBuilder::LaserData & pointsdata,vector<int> scoreidex_)
        : pose_matched(pose),pointsValueIdex_(pointsValueIdex){ ScanPoints_ = pointsdata; ScoreIdex = scoreidex_;}
        MatchResult(){}
        MapBuilder::Pose pose_matched;
        std::vector<float> pointsValueIdex_;
        MapBuilder::LaserData ScanPoints_;
        vector<int> ScoreIdex;
    };

    struct MatchResult_v02
    {
        MatchResult_v02(const MapBuilder::Pose& pose,const vector<Array2d> & laserdata)
            : pose_matched(pose){
            LaserScanData_ptr = make_shared<vector<Array2d>>(vector<Array2d>(laserdata));
        }
        MatchResult_v02() {}
        MapBuilder::Pose pose_matched;
        shared_ptr<vector<Array2d>> LaserScanData_ptr;
        
    };
    // I just Add a _ at the end of the struct name. Then everything become ok.
    struct LocalMap_
    {
        LocalMap_ ( const LaserData & localmap,const  Array2d &  minXY_ , const Array2d & maxXY_)
                :localmap_points(localmap),minXY(minXY_),maxXY(maxXY_){}
        LocalMap_(){}
        LaserData localmap_points;
        Array2d minXY;
        Array2d maxXY;
    };
    // the set collect all the colors that will be used in the Drawer function
    struct ScalarSet
    {
        ScalarSet( const cv::Scalar & wholeampcolor = cv::Scalar(0,0,0),const cv::Scalar & scancolor = cv::Scalar(255,165,0), const cv::Scalar & pathcolor = cv::Scalar(0,250,0), const cv::Scalar & carcolor = cv::Scalar(0,0,200), const cv::Scalar & backgroundcolor = cv::Scalar(255,255,255))
        :WholeMapColor(wholeampcolor),ScanColor(scancolor),PathColor(pathcolor),CarColor(carcolor),BackGroundColor(backgroundcolor){}
        cv::Scalar WholeMapColor;
        cv::Scalar ScanColor;
        cv::Scalar PathColor;
        cv::Scalar CarColor;
        cv::Scalar BackGroundColor;
    };
    // the map Parameters
    struct MapParameter
    {
        int BorderSize;
        double RoughpixelSize;
        double FinerpixelSize;
        bool miniUpdated;
        double miniUpdatedDT;
        double miniUpdateDR;
        double range_eff_max;
        double range_eff_min;
        Eigen::Array3d fastResolution;
    };
    MapBuilder () {}
    
    // Set the simulating Laser's Parameter  
    // 设置参数
    void SetLaserParameter (); 
    
    // 设置参数 
    void SetMapParameter(); 
    
    // 得到所有Scan的数据
    bool GetAllData(const string,map<unsigned int,MapBuilder::LaserData> & AllLaserData);
    // ReadAScam is used to extract the LaserData from the file.
    // The raw data consist of range and angle;
    // And this funciton consist of some data preprocess 
    //i.e. extract the data in specific range and Transform them into a Cartesian Coordinate 
    MapBuilder::LaserData ReadAScan(unsigned int ScanIdx,map<unsigned int,MapBuilder::LaserData> & AllLaserData, MapBuilder::MapParameter);
        
    // VoxelFilter() is used to reduce the quantity of sacn points;
    // this is just a approximate voxel filter. the real voxel filter should extract the mass mean of points in the voxel cell
    MapBuilder::LaserData &  VoxelFilter (const MapBuilder::LaserData &, double);

    // this function is trying to build a set of scan's near points
    MapBuilder::LocalMap_ ExtractLocalMap(const MapBuilder::LaserData & Allpoints, const MapBuilder::LaserData & filtered_pointCloud, MapBuilder::Pose pose_ ,unsigned int BorderSize);
    
    MapBuilder::LocalMap_ ExtractLocalMap_Kdtree(const MapBuilder::LaserData& Allpoints, const MapBuilder::LaserData& filtered_pointCloud
        , MapBuilder::Pose pose_, unsigned int BorderSize);
    // ComputeOccGrid is used to generate a Occupied Map of the local map
    MapBuilder::OccGrid ComputeOccGrid( const MapBuilder::LocalMap_ & local_map_, double pixelSize);
    
    // prepare for the first keyscan
    bool Initialization( const MapBuilder::LaserData & CurrentScan, const MapBuilder::Pose CurrentPose);
    
    // guess a new pose
    MapBuilder::Pose RoughPoseEstimation(Pose & CurrentPose, const unsigned int & ScanIdx) const;
    
    MapBuilder::MatchResult HillClimbMatcher(const MapBuilder::OccGrid & OccGridMap, const Eigen::Array3d & MatcherResolution, 
                                             MapBuilder::LaserData & ScanData, MapBuilder::Pose pose_guess,
                                             double pixelSize);
    
    void AddKeyScan( const MapBuilder::MatchResult & matchresult_);
    void AddKeyScan_v02(MapBuilder::MatchResult_v02& matching_result);
    void Start();
    // this function is based on Opencv to draw something
    void Drawer(const MapBuilder::LaserData& Allpoints, double pixelSize, const MapBuilder::ScalarSet scalar_set_);

    // 回环的内容
    void updatedata(const vector<Array2d>& voxel_filtered_scan_, Array3d pose_, Array2i map_index, int ScanIdx);

    //private:
    vector<MapBuilder::LaserData> LaserDataPool_;
    vector<vector<Array2d>> V02_keyscan_pool;
    vector<Array3d> PosePool;
    vector<KeyScan>  keyscans_;
    map<unsigned int,LaserData> AllLaserData;
    LaserParameters LaserParameter_;
    MapParameter MapParameter_;
    vector<LaserData> PointsPool;
    LaserData AllPoints;
    vector<Pose> path;
    cv::Mat whole_map;
    MapBuilder::OccGrid OccGridRough;
    MapBuilder::OccGrid OccGridFiner;
    MapBuilder::MatchResult matchresult;
    ScalarSet ScalarSet_;
    tree::Kdtree kdtree_;
    int ScanIdx;
    // this function just return distance difference of pose
    Pose DiffPose(const vector<Pose> & path, Pose Current) const;
    
    // this function is used to transfor the points' from lcoal frame to world frame
    LaserData Transform(const LaserData & UnprocessScanData, Pose pose_);
    // used in ComputeOccGrids to discrete the points coordiante
    vector<Array2i> DiscretePoints(const LaserData & PointCloud, double pixelSize) const ;
    
    // Drawer Basic Funtion 画图的函数 不是都用到了
    void DrawPoint(cv::Mat& map_, const LaserData& laserdata_, double pixelSize, const cv::Scalar& scalar_, const Array2d& _offsetXY);
    void DrawLine(cv::Mat & map_, Eigen::Array2i start_,Eigen::Array2i end_, const cv::Scalar & scalar_);
    void DrawLine(cv::Mat & map_, cv::Point start_,cv::Point end_, const cv::Scalar & scalar_);
    void DrawPath(cv::Mat & map_,const  vector<MapBuilder::Pose> & path , double pixelSize, const cv::Scalar & scalar_);
    void DrawCar(cv::Mat& map_, const MapBuilder::Pose & pose_, const cv::Scalar & scalar_,double pixelSize);
    Mat self_distanceTransform(Mat& pic_, int size_ = 5);
    Mat MatsAndOperator(Mat& mat1, Mat& mat2);
    double calculateL2distance(int x1, int y1,int x2,int y2);
};





#endif
