#ifndef _PGM_H_
#define _PGM_H_
#include "common.h"

// All function have been test;
// 还没有进行整合调试

using namespace Eigen;
using namespace cv;
using namespace std;
// odd occupied probability grid map
class PGMBuilder
{
public:

    struct ProbabilityAParameter
    {
        float lower_bound;
        float high_bound;
        float HitOdd;
        float MissOdd;
        float unknown;
    };

    PGMBuilder() {}
    void  InsertFirstScan(const LaserFan& laserfan_, double pixelsize_);

    // Main part of PGM class Odd 
    bool Insert(LaserFan& laserfan_);
    bool CastRays(const Array2i origin, vector<Array2i> ends);
    bool CastRay(const Eigen::Array2i origin, const Eigen::Array2i end, const vector<uint16_t>& table_);
    bool ApplyOdd(const Eigen::Array2i XY, const vector<uint16_t>& table_);
    void GrowLimitsIfNeeded(const Vector2i minXY, const Vector2i maxXY);
    void Finished();

    // return the value of PGM
    // These Four functions are test passed;
    const vector<uint16_t>& GetProbabilityValue() const { return ProbabilityValue; }
    MapLimits GetMapLimits()const { return maplimits_; }
    AlignedBox2i GetAlignedBox()const { return scale_; }
    vector<float> GetProbabilityMap() const;
    double GetPixelSize()const { return pixelSize_; }
    // the XY is in the local frame
    float GetProbability_LocalCandidates(const Eigen::Array2i& XY_) const
    {
        assert(maplimits_.Contain(XY_));
        return 1.f - ValueToCorrespondenceCost[ProbabilityValue[XY_(1) * maplimits_.scaleXY(0) + XY_(0)]];
    }
    float GetProbability_WorldCandidates(const Eigen::Array2i& XY_) const
    {
        Eigen::Array2i temp = XY_ - Array2i(scale_.min());
        assert(maplimits_.Contain(temp));
        return 1.f - ValueToCorrespondenceCost[ProbabilityValue[temp(1) * maplimits_.scaleXY(0) + temp(0)]];
    }
    float GetProbability_WorldCandidates(const Eigen::Array2d& XY_) const
    {
        Eigen::Array2i temp = Array2i(std::round(XY_(0)), std::round(XY_(1))) - Array2i(scale_.min());
        assert(maplimits_.Contain(temp));
        return 1.f - ValueToCorrespondenceCost[ProbabilityValue[temp(1) * maplimits_.scaleXY(0) + temp(0)]];
    }
    int GetSubmapIdx()const { return SubmapIdx; }
    const AlignedBox2i& GetScale()const { return scale_; }

    static vector<float> ValueToProbabilityChart;
    static vector<float> ValueToCorrespondenceCost;

    // TEST PASS
    // Probability Basic Function for PGM
    // the functions below have been check that there is no vector::resize is called;
    // 这里的这个initialization应该是一个静态函数的
    static bool Initialization();
    // These functions Test Pass 
    static float Odd (float ProbabilityValue) {return ProbabilityValue/(1.f-ProbabilityValue);}
    static float ProbabilityFromOdd (float Odds) {return Odds/(1.f +Odds);}
    static float ProbabilityToCorrespondenceCost(const float probability) {return 1.f - probability;}
    static float CorrespondenceCostToProbability(const float correspondence_cost) {return 1.f - correspondence_cost;}
    // This Function Test Pass
    static float ClampProbability(const float probability);
    // This function Test Pass
    static uint16_t ProbabilityToValue(const float probability);
    // This function is same to the function above
    static uint16_t CorrespondenceCostToValue( const float correspondence_cost);
    // These three functions are pass;
    static vector<uint16_t> ComputetoApplyProbabilityChart(float odds);
    static vector<uint16_t> ComputetoApplyCorrespondenceCostChart(float odds);
    static void ComputeValueToProbabilityAndCorrespondCostChart();



    vector<uint16_t> ProbabilityValue;
    vector<uint16_t> CellIndex;
    Array2d origin;
    MapLimits maplimits_;
    uint16_t laserfan_num;
    AlignedBox2i scale_;
    int SubmapIdx;
    double pixelSize_;
    static int SubmapIdx_count;
    static ProbabilityAParameter probabilityparameter_;
    static bool InitializationFlag;
    static vector<uint16_t> HitAddChartProbability;
    static vector<uint16_t> MissAddChartProbability;
    static vector<uint16_t> HitAddChartCorrespondenceCost;
    static vector<uint16_t> MissAddChartCorrespondenceCost;
};

// probability using the Correspond Value to describ



#endif
