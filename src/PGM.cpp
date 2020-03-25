#include "PGM.h"
uint16_t kUpdateMaker = uint16_t(1)<<15;
PGMBuilder::ProbabilityAParameter PGMBuilder::probabilityparameter_;
// 其实这两个Probability有关的函数，之后用不到
vector<uint16_t> PGMBuilder::HitAddChartProbability;
vector<uint16_t> PGMBuilder::MissAddChartProbability;

vector<uint16_t> PGMBuilder::HitAddChartCorrespondenceCost;
vector<uint16_t> PGMBuilder::MissAddChartCorrespondenceCost;
vector<float> PGMBuilder::ValueToProbabilityChart;
vector<float> PGMBuilder::ValueToCorrespondenceCost;
bool PGMBuilder::InitializationFlag;
int PGMBuilder::SubmapIdx_count = 0;



// laserfan_ should been represented by the world frame
// TEST PASS
void PGMBuilder::InsertFirstScan(const LaserFan & laserfan_, double pixelsize_)
{
    SubmapIdx = SubmapIdx_count;
    assert(InitializationFlag == true);
    origin = Array2d(laserfan_.pose(0),laserfan_.pose(1));
    pixelSize_ = pixelsize_;
    // give value to PGMBuilder origin to generate a Localmap
    
    // Discrete the point cloud and origin;
    vector<Array2i> discreend = BasicFunction::DiscretePoints(laserfan_.end,pixelsize_);
    Vector2i dis_origin = Vector2i(round(origin(0)/pixelsize_),round(origin(0)/pixelsize_));
    Array2i offsetXYtoOrigin (0,0);
    // generate the Maplimits
    Array2i minXY(discreend[0]);
    Array2i maxXY(discreend[0]);
    for (auto& element : discreend)
    {
        minXY = minXY.min(element);
        maxXY = maxXY.max(element);
    }
    laserfan_num = 1;
    scale_.extend(dis_origin);
    scale_.extend(Vector2i(minXY));
    scale_.extend(Vector2i(maxXY));
    
    Vector2i BottomLeft = scale_.min();
    Vector2i TopRight = scale_.max(); 
    int CellNumber = (TopRight(0) - BottomLeft(0) + 1) * (TopRight(1) - BottomLeft(1) + 1);
    ProbabilityValue.resize(CellNumber,uint16_t(0));
    maplimits_ = MapLimits(pixelsize_, Array2i(TopRight - BottomLeft +Vector2i(1,1)));
    
    // the main part of the PGM Initialization
    CastRays(dis_origin,discreend);
    SubmapIdx_count++;
    Finished();
}

// laserfan_ should been represented by the world frame
// Insert the data into the probability map 
bool PGMBuilder::Insert(LaserFan & laserfan_ )
{
    //
    assert(InitializationFlag == true);
    double pixelsize_ = maplimits_.pixelSize;
    vector<Array2i> discreend = BasicFunction::DiscretePoints(laserfan_.end,pixelsize_);
    Vector2i dis_origin = Vector2i(round(laserfan_.pose(0)/pixelsize_),round(laserfan_.pose(1)/pixelsize_));
    // generate the Maplimits
    // figure out if new maxXY and minXY exist  in the submap
    Array2i minXY(discreend[0]);
    Array2i maxXY(discreend[0]);
    for (auto& element : discreend)
    {
        minXY = minXY.min(element);
        maxXY = maxXY.max(element);
    }
    minXY = minXY.min(Array2i(dis_origin));
    maxXY = maxXY.max(Array2i(dis_origin));
    laserfan_num ++;
    GrowLimitsIfNeeded(Vector2i(minXY),Vector2i(maxXY));
    CastRays(dis_origin,discreend);
    Finished();
    return true;
}

// minXY maxXY origin are in the world frame now.
// 基本上逻辑是对的
void PGMBuilder::GrowLimitsIfNeeded(const Eigen::Vector2i minXY, const Eigen::Vector2i maxXY)
{
    // transfer these two points into the localmap frame
    Vector2i minXY_localmap, maxXY_localmap;
    minXY_localmap = minXY - scale_.min();
    maxXY_localmap = maxXY - scale_.min();
    Vector2i MinExtendArea(0,0);
    Vector2i MaxExtendArea(0,0);
    bool SmallerMin = 0 ;
    bool LargerMax = 0;
    // Check whether all these three points are in the localmap
    if(!maplimits_.Contain(minXY_localmap))
    {
        Array2i temp_min(minXY);
        temp_min = temp_min.min(Array2i(scale_.min()));
        MinExtendArea = 2 * (scale_.min()- Vector2i(temp_min));
        SmallerMin = 1;
    }
    if(!maplimits_.Contain(maxXY_localmap))
    {
        Array2i temp_max(maxXY);
        temp_max = temp_max.max(Array2i(scale_.max()));
        MaxExtendArea = 2 * (Vector2i(temp_max) - scale_.max());
        LargerMax = 1; 
    }
    if(SmallerMin ==0 && LargerMax == 0)
        return;
    MapLimits new_limits;
    vector<uint16_t> newtable;

    // Generate the new Probability Value 
    if(SmallerMin == 1 || LargerMax == 1)
    {    
        scale_.extend(scale_.min()- MinExtendArea);
        scale_.extend(scale_.max()+MaxExtendArea);
        new_limits.scaleXY = Array2i(scale_.max() - scale_.min()) + Array2i(1,1);
        new_limits.pixelSize = maplimits_.pixelSize;
        newtable.resize((new_limits.scaleXY(0))*(new_limits.scaleXY(1)), uint16_t(0));
        uint32_t stride = new_limits.scaleXY(0);
        uint32_t offset_bottom = stride * MinExtendArea(1);
        uint32_t offset_leftMargin = MinExtendArea(0);
        uint32_t Old_stride = maplimits_.scaleXY(0);

        // Update the whole localMap
        for(uint32_t row = 0;row <maplimits_.scaleXY(1);row++)
        {
            for (uint32_t col = 0; col<maplimits_.scaleXY(0);col++)
            {
                uint32_t temp_x =offset_leftMargin + col;
                newtable[offset_bottom + row*stride + temp_x ] = ProbabilityValue[row*Old_stride + col];
            }
        }
    }
    ProbabilityValue = newtable;
    maplimits_ = new_limits;
    
}

bool PGMBuilder::CastRays(const Eigen::Array2i origin, vector<Eigen::Array2i> ends )
{
    //
    Vector2i BottomLeft = scale_.min();
    Array2i minValue(BottomLeft);
    // Update the Value for the MAP ends first
    for (auto & element : ends)
    {
        assert(maplimits_.Contain(Array2i(element-minValue)));
        ApplyOdd(element-minValue,HitAddChartCorrespondenceCost);
    }
    
    // Update the Value of cell which is on the CastRays
    for (auto & element :ends)
    {
        CastRay(origin-minValue,element-minValue,MissAddChartCorrespondenceCost);
    }
    return true;
}

// orgin and end should in the frame of localmap
// orgin's x should smaller than end's x
// both xy of orgin and end are positive as these two points have been transfered into the loca;map frame
bool PGMBuilder::CastRay(const Eigen::Array2i origin, const Eigen::Array2i end, const vector<uint16_t> & table_)
{
    assert(maplimits_.Contain(origin)&&maplimits_.Contain(end));
    ApplyOdd(origin,table_);
    
    //cout << PGMBuilder::ValueToCorrespondenceCost[ProbabilityValue[BasicFunction::FlatToIndex(origin, maplimits_)]-32768] << endl;
    
    // 保证这个x 计数方式是从左到右的
    if (origin(0)>end(0))
    {
        CastRay(end,origin,table_);
        return true;
    }
    /*cout << " the origin is " << origin(0) << "," << origin(1) << endl;*/
    // 为了保证代码的鲁棒性 防止ratio->无穷
    if (end(0) == origin(0))
    {
        int temp = end(1) - origin(1);
        if (temp == 0)
            return true;
        int signal = temp / abs(temp);
        for (int i = 0; i <= abs(temp); i++)
        {
            ApplyOdd(Array2i(origin(0), origin(1) + i * signal), MissAddChartCorrespondenceCost);
        }
        return true;
    }


    float ratio = float(end(1)-origin(1))/float(end(0) - origin(0));
    float length_x = 0.5;
    float length_y = ratio/2;
    //cout << "The origin is " << origin(0) << "," << origin(1) << endl;
    //cout << "The end is " << end(0) << "," << end(1) << endl;
    //cout << "The axis go by point " << endl;
    while(length_x<(end(0)-origin(0)))
    {
        float coordinate_x = std::round(length_x - 1e-7);// Here is the function of KPadding in the cartographer
        int sub_y = std::round(length_y - ratio/2);
        int up_y = std::round(length_y);

        // 还是有一些不必要的计算，但是效果是可以实现的
        for(;sub_y <= up_y; sub_y++)
        {
            ApplyOdd(Array2i(std::round(coordinate_x + origin(0)),sub_y + origin(1)),table_);
            
        }
        length_x += 0.5;
        length_y += ratio/2;
    }
    return true;
    
}

// XY 是点云再localmap 中的坐标
bool PGMBuilder::ApplyOdd(const Eigen::Array2i XY, const vector<uint16_t> &  table_)
{
    assert(maplimits_.Contain(XY));
    assert(InitializationFlag == true);
    int stride = maplimits_.scaleXY(0);
    assert(maplimits_.Contain(XY));
    int index = XY(1) * stride + XY(0);
    if (ProbabilityValue[index]<kUpdateMaker)
    {
    ProbabilityValue[index] = table_[ProbabilityValue[index]];
    return  true;
    }
    else
    {
        return false;
    }
}
void PGMBuilder::Finished()
{
    for(auto & element: ProbabilityValue)
    {
        if (element >kUpdateMaker)
            element -= kUpdateMaker;
    }
}

vector<float> PGMBuilder::GetProbabilityMap() const
{
    vector<float> temp;
    temp.reserve(ProbabilityValue.size());
    for (auto & element : ProbabilityValue)
    {
        temp.push_back(1.f - ValueToCorrespondenceCost[element]);
    }
    return temp;
}

// use to trunc the probability
inline  float PGMBuilder::ClampProbability(const float probability) {
    if (probability > probabilityparameter_.high_bound)
        return probabilityparameter_.high_bound;
    if (probability<probabilityparameter_.lower_bound)
        return probabilityparameter_.lower_bound;
    return probability;
}
// transfer the probability to integer value
inline uint16_t PGMBuilder::ProbabilityToValue(const float probability)
{
    uint16_t value =
      round(
          (ClampProbability(probability) - probabilityparameter_.lower_bound) *
          (32766.f / (probabilityparameter_.high_bound - probabilityparameter_.lower_bound))) +
      1;
      return value;
}
// transfer the correspondence_cost to value
inline uint16_t PGMBuilder::CorrespondenceCostToValue(const float correspondence_cost)
{
     uint16_t value =
      round(
          (ClampProbability(correspondence_cost) - probabilityparameter_.lower_bound) *
          (32766.f / (probabilityparameter_.high_bound - probabilityparameter_.lower_bound))) +
      1;
      return value;
}
// Based on Corresponding represent i.e. the possibility of freedom;
vector<uint16_t> PGMBuilder::ComputetoApplyProbabilityChart(float odds)
{
    vector<uint16_t> temp_result;
    // include the possibility unkonwn so there are 1<<15 possibilities
    temp_result.reserve(uint16_t(1)<<15);
    // in order to get a uniform form we multiply a 1 behind odds
    temp_result.push_back(ProbabilityToValue(ProbabilityFromOdd(odds*1)) + kUpdateMaker) ;
    for (uint16_t ProbabilityToValue_ = 1; ProbabilityToValue_!= (uint16_t(1)<<15);ProbabilityToValue_++)
    {
        temp_result.push_back(ProbabilityToValue(ProbabilityFromOdd(odds * Odd(ValueToProbabilityChart[ProbabilityToValue_])))+ kUpdateMaker);
    }
    assert(temp_result.size() == (uint16_t(1)<<15));
    return temp_result;
}

vector<uint16_t> PGMBuilder::ComputetoApplyCorrespondenceCostChart(float odds)
{
    vector<uint16_t> temp_result;
    temp_result.reserve(uint16_t(1)<<15);
    temp_result.push_back(CorrespondenceCostToValue(
				ProbabilityToCorrespondenceCost(ProbabilityFromOdd(
					odds ))) + kUpdateMaker);
    int cell = 1;
    for (; cell!=(uint16_t(1)<<15);cell ++)
    {
        
        temp_result.push_back(CorrespondenceCostToValue(
				ProbabilityToCorrespondenceCost(ProbabilityFromOdd(
					odds * Odd(CorrespondenceCostToProbability(
					ValueToCorrespondenceCost[cell]))))) + kUpdateMaker);
    }
    return temp_result;
}
// get the relationship between value and probability and correspondence_cost 
void PGMBuilder::ComputeValueToProbabilityAndCorrespondCostChart()
{
    vector<uint16_t> temp_result;
    temp_result.reserve(uint16_t(1)<<15);
    ValueToProbabilityChart.reserve(32768);
    ValueToProbabilityChart.push_back(0.5f);
    for (int cell = 1;cell!= (int(1)<<15);cell++)
    {
        float temp = (probabilityparameter_.high_bound - probabilityparameter_.lower_bound)*((cell-1)/32766.f)+ probabilityparameter_.lower_bound;
        ValueToProbabilityChart.push_back(temp);
    }
    ValueToCorrespondenceCost = ValueToProbabilityChart;
}

 bool PGMBuilder::Initialization()
{
    // set the parameters
   
    probabilityparameter_.high_bound = 0.9;
    probabilityparameter_.lower_bound = 0.1;
    probabilityparameter_.HitOdd = 0.55/0.45;
    probabilityparameter_.MissOdd = 0.45/0.55;
    probabilityparameter_.unknown = 2.f;
    
    //Build the Table or Chart for the following using 
    
    ComputeValueToProbabilityAndCorrespondCostChart();
    HitAddChartProbability = ComputetoApplyProbabilityChart(probabilityparameter_.HitOdd);
    MissAddChartProbability = ComputetoApplyProbabilityChart(probabilityparameter_.MissOdd);
    HitAddChartCorrespondenceCost = ComputetoApplyCorrespondenceCostChart(probabilityparameter_.HitOdd);
    MissAddChartCorrespondenceCost = ComputetoApplyCorrespondenceCostChart(probabilityparameter_.MissOdd);

    InitializationFlag = 1;
    return true;
}
 