#include "common.h"

namespace BasicFunction
{
vector<Array2i> DiscretePoints(const LaserData & PointCloud, double pixelSize)
{
    vector<Eigen::Array2i> discretePointCloud;
    discretePointCloud.reserve(PointCloud.size());

    for (auto& point : PointCloud)
    {
        Array2d temp = Array2d(point) / pixelSize;
        int x = std::round(temp(0));
        int y = std::round(temp(1));
        // 大场景下可能触发以下的assert
        assert(x<32767 && x>-32767);
        assert(y<32767 && y>-32767);
        discretePointCloud.push_back(Array2i(x, y));
    }
    return discretePointCloud;
}
LaserData Transform(const LaserData& UnprocessScanData, Pose pose_)
{
    LaserData transformed_scandata;
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
uint32_t FlatToIndex(const Array2i temp, const MapLimits & temp_limits)
{
    return temp(1)* temp_limits.scaleXY(0) + temp(0);
}
}
