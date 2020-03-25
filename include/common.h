#ifndef _COMMON_H_
#define _COMMON_H_
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include "Eigen/Eigen"
#include <fstream>
#include <regex>
#include <map>
#include <iterator>
#include <math.h>
#include <set>
#include <unordered_set>
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/core/core.hpp"        
#include "opencv2/highgui/highgui.hpp"
#include <assert.h>
#include <bitset>
#include <deque>
#include <map>
#include <array>
#include "Tree.h"

using namespace Eigen;
using namespace std;
using namespace cv;
typedef vector<Array2d> LaserData;
typedef vector<Array2i> discreLaserData;
// pose(0) ->x pose(1) ->y pose(2)->theta
typedef Eigen::Array3d Pose;

// �Ժ�����ĸ�ֵ����ĳ�ָ������������
struct LaserFan
{
    LaserFan() {}
    LaserFan( std::vector<Array2d> ends, const Array3d& pose_) : pose(pose_) { end = ends; }
     std::vector<Array2d>  end;
    Array3d pose;
};
struct MapLimits
{
    MapLimits() {}
    MapLimits(double pixelsize_, Array2i scaleXY_)
    {
        pixelSize = pixelsize_;
        scaleXY = scaleXY_;
    }

    // ���ﴫ�������Ǵӣ�0��0����ʼ��������������
    bool Contain(Eigen::Array2i XY)  const
    {
        return XY(0) >= 0 && XY(1) >= 0 && XY(0) < scaleXY(0) && XY(1) < scaleXY(1);
    }
    bool Contain(Eigen::Vector2i XY) const
    {
        return XY(0) >= 0 && XY(1) >= 0 && XY(0) < scaleXY(0) && XY(1) < scaleXY(1);
    }

    double pixelSize;
    // �����scaleXY ����һ������index ����һ��һ�л���һ�а�����������
    Array2i scaleXY;

    // from ZERO to the max_x and max_y cell
    // eg. x,y = 4.3,4.3 in the cell which index is 5,5

};
struct Submap_Scan
{
    struct Node_info
    {
        enum status { inside = 1, outside =0 };
        Node_info(status status__, const Array3d& world_pose_,float score_) :status_(status__), 
            world_pose(world_pose_),score(score_) {}
        status status_;
        Array3d world_pose;
        float score;
    };
    Submap_Scan() {}
    explicit Submap_Scan(int SubmapIdx_) : SubmapIdx(SubmapIdx_)
    {
        // ��ǰԤ����λ��
        ScanIdxPool.reserve(10);
        mean_pose = Array3d::Zero();
        total_pose = Array3d::Zero();
    }
    bool Insert(int SubmapIdx_,int ScanIdx, Submap_Scan::Node_info::status status_,
        const Array3d& world_pose,const float score_)
    {
        assert(SubmapIdx_ == SubmapIdx);
        auto temp = this->Node_info_pool.insert(pair<int, Node_info>(ScanIdx, Node_info(status_, world_pose,score_)));
        if (temp.second == true)
        {
            ScanIdxPool.push_back(ScanIdx);
            total_pose += world_pose;
            mean_pose = total_pose / ScanIdxPool.size();
        }
        return temp.second;
    }
    vector<int> ScanIdxPool;
    bool UpdateScale_Min(const Array2i& scale_min_) { scale_min = scale_min_; return true; }
    int SubmapIdx;
    Array2i scale_min;
    double pixelSize;
    map<int, Node_info> Node_info_pool;
    Array3d mean_pose;
    Array3d total_pose;
};
namespace BasicFunction
{
    vector<Array2i> DiscretePoints(const LaserData& PointCloud, double pixelSize);
    LaserData Transform(const LaserData& UnprocessScanData, Pose pose_);
    uint32_t FlatToIndex(const Array2i temp, const MapLimits& temp_limits);
    template <typename T>
    void FindMaxAndMin(T& min_, T& max_, const vector<T>& table_)
    {
        min_ = table_[0];
        max_ = table_[0];
        assert(table_.size() != 0);
        for (auto& element : table_)
        {
            min_ = min_.min(element);
            max_ = max_.max(element);
        }
    }
    template<typename T>
    vector<T> VoxelFilter(const vector<T>& Laser_data_, double pixelSize)
    {
        using keyType = std::bitset<64 * 2>;
        unordered_set<keyType> voxel_set_;
        vector<T> filtered_pointCloud;

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
            // ��������������ͬ��ֵ if ��������ͨ��
            if (it_inserted.second)
            {
                filtered_pointCloud.push_back(element);
            }
            count++;

        }
        return filtered_pointCloud;
    }
}


#endif
