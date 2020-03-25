#include "M2MPGM.h"

int M2MPGM::max_height = 8;

M2MPGM::M2MPGM(shared_ptr<PGMBuilder> pgmbuilder_,int height)
{
    pixelSize_ = pgmbuilder_->GetPixelSize();
    // 寻找合适的层数
    maplimits_ = pgmbuilder_->GetMapLimits();
    scale_ = pgmbuilder_->GetAlignedBox();
    int scale_min = maplimits_.scaleXY(0) > maplimits_.scaleXY(1) ? maplimits_.scaleXY(1) : maplimits_.scaleXY(0);
    int count = 0;
    while (int(scale_min / 2) != 0)
    {
        count++;
        scale_min = scale_min / 2;
    } 

    if (height > count)
        height = count;

    if (height>max_height)

    PGMset_.reserve(height);
    PGMset_.push_back(pgmbuilder_->GetProbabilityMap());

    for (int Curr_height = 1; Curr_height<=height;Curr_height ++)
    {
        PGMset_.push_back(ManyPixelPGMCreate(Curr_height, pgmbuilder_->GetProbabilityValue()
                          ,pgmbuilder_->ValueToCorrespondenceCost));
    }
    height_ = height;
    CreatedFinished = true;
}
M2MPGM::M2MPGM(const PGMBuilder& pgmbuilder_, int height)
{
    pixelSize_ = pgmbuilder_.GetPixelSize();
    // 寻找合适的层数
    maplimits_ = pgmbuilder_.GetMapLimits();
    scale_ = pgmbuilder_.GetAlignedBox();
    int scale_min = maplimits_.scaleXY(0) > maplimits_.scaleXY(1) ? maplimits_.scaleXY(1) : maplimits_.scaleXY(0);
    int count = 0;
    while (int(scale_min / 2) != 0)
    {
        count++;
        scale_min = scale_min / 2;
    }

    if (height > count)
        height = count;

    if (height > max_height)
        height = max_height;

    PGMset_.reserve(height);
    PGMset_.push_back(pgmbuilder_.GetProbabilityMap());

    for (int Curr_height = 1; Curr_height <= height; Curr_height++)
    {
        PGMset_.push_back(ManyPixelPGMCreate(Curr_height, pgmbuilder_.GetProbabilityValue()
            , pgmbuilder_.ValueToCorrespondenceCost));
    }
    height_ = height;
    CreatedFinished = true;

}

// 这个地方应该要有优化一下存储结构
// Origin_table_ is the value table_, Corres_table_ is a table from value to probability
vector<float> M2MPGM::ManyPixelPGMCreate(int height,const vector<uint16_t> & Origin_table_ ,
    const vector<float> & Corres_table_)
{
    height_ = height;
    assert(height>=0);
    int stride = int (1)<<uint16_t(height);
    int max_x_border = maplimits_.scaleXY(0);
    int max_y_border = maplimits_.scaleXY(1);
    // Limit the scale of stride;
    //assert(stride<=floor(maplimits_.scaleXY(0)/10));
    //assert(stride<=floor(maplimits_.scaleXY(1)/10));
    SlideWindow.clear();
    vector<float> new_table_;
    vector<float> new_table_result_;
    // 扩大地盘 在地图的左边和上边各自扩张stride-1 个格子
    new_table_.resize((maplimits_.scaleXY(0) + stride - 1 )*(maplimits_.scaleXY(1)+ stride -1 ), 0.f);
    new_table_result_.resize((maplimits_.scaleXY(0) + stride - 1) * (maplimits_.scaleXY(1) + stride - 1), 0.f);
    int new_scale_x = maplimits_.scaleXY(0) + stride - 1;
    int new_scale_y = maplimits_.scaleXY(1) + stride - 1;
    
    
    // handle the row first 
    // Step1:find the max in the x0->x0 + stride - 1 in the X-asix
    // Step2:find the max in the y0 ->y0 +stride -1 in the Y-asix
    // This Cycle is for Step 1
    int stride_x = stride > max_x_border ? max_x_border : stride;
    for (int row = 0; row < maplimits_.scaleXY(1);row++)
    {
        // 把滑窗填满
        for (int col_count = 1 ; (col_count <= stride) ;col_count++)
        {
            uint32_t temp = Origin_table_[row * max_x_border + max_x_border - col_count];
            AddValue(1.f - (Corres_table_)[temp]);
            new_table_[row * new_scale_x + new_scale_x - col_count] = GetMaxProb();
        }
        
        // TODO 地图非常小的时候会出现，max_x_border<stride 的情况会报错 这里要处理一下
        for (int col_count = stride + 1 ; col_count <= max_x_border + stride -1; col_count++ )
        {
            uint32_t temp_1 = Origin_table_[row *  max_x_border + max_x_border - col_count +stride];
            RemoveValue(1.f - (Corres_table_)[temp_1]);
            if (col_count>max_x_border)
            {
                AddValue(0.f);
            }
            else
            {
                uint32_t temp = Origin_table_[row * max_x_border + max_x_border - col_count];
                AddValue(1.f - (Corres_table_)[temp]);
            }
            new_table_[row * new_scale_x + new_scale_x - col_count] = GetMaxProb();
        }
        SlideWindow.clear();
    }

    SlideWindow.clear();
    
    // This cycle is for step 2
    for (int col = 0; col< new_scale_x;col++)
    {
       // 填滑窗
        for (int row_count = 1; row_count <= stride ; row_count ++)
        {
            float temp = new_table_[(row_count-1)*new_scale_x + col];
            AddValue(temp);
            new_table_result_[(row_count - 1)* new_scale_x + col] = GetMaxProb();
            //cout << GetMaxProb() << endl;
        }
        
        //TODO 还是一样需要进行检测
        for (int row_count = stride + 1; row_count <= new_scale_y ; row_count++)
        {
            RemoveValue(new_table_[(row_count - 1 - stride) * new_scale_x + col]);
            float temp = new_table_[(row_count - 1) * new_scale_x + col];
            AddValue(temp);
            new_table_result_[(row_count - 1) * new_scale_x + col] = GetMaxProb();
        }
        SlideWindow.clear();
    }
    return new_table_result_;
}

void M2MPGM::AddValue(float temp)
{
    while(!SlideWindow.empty()&& SlideWindow.back()<temp)
    {
        SlideWindow.pop_back();
    }
    SlideWindow.push_back(temp);
}

void M2MPGM::RemoveValue(float temp)
{
    if (temp == SlideWindow.front())
        SlideWindow.pop_front();
}

float M2MPGM::GetMaxProb()
{
    return SlideWindow.front();
}

