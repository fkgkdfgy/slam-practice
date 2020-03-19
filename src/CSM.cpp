#include "CSM.h"

// Scan data should be represent in the local frame
Array3d CSM::Match(const LaserData & laserdata_,Array3d pose_guess_,const PGMBuilder & PGM_ptr )
{
    // imply the complete of GenertateSearchWindow
    assert(initializationFlag == true);
    
    // Preparing for the Scoring Part
    vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans(laserdata_,pose_guess_);
    
    //cout << " generate the pose angle " << endl;
    //for (auto& element : rotateDiscreteResults)
    //{
    //    cout << element.pose(2) << "   ";
    //}
    //cout << endl;
    vector<Array2i> candidates = GenerateCandidate();
    CSM::MatchResult BestResult(0.f,Array3d(0,0,0));
    double pixelSize_ = SearchWindowParameter_.translation_step_size;
    
    // Inplement the Scoring Part
    for (CSM::DiscreteRotatedResult & rotateDiscreteResult : rotateDiscreteResults )
    {
         CSM::MatchResult temp = ScoreAllCandidates(PGM_ptr.ProbabilityValue,rotateDiscreteResult.DiscreteLaserData,
                                                    PGM_ptr.GetMapLimits(), PGM_ptr.GetAlignedBox(),candidates);
         if (temp.bestscore > BestResult.bestscore)
         {
             BestResult.bestscore = temp.bestscore;
             BestResult.bestpose = rotateDiscreteResult.pose + Array3d(temp.bestpose(0),temp.bestpose(1),0);
         }
    }
    // 真实微小运动和估值相反
    cout << " the best average score is " << endl;
    cout << BestResult.bestscore/laserdata_.size() << endl;
    cout << " the matched pose is " << endl;
    cout << BestResult.bestpose(0) << "  " << BestResult.bestpose(1) << "  " << BestResult.bestpose(2) << endl;
    return BestResult.bestpose;
}
void CSM::GenerateSearchWindow(const LaserFan & laserfan_, double linear_window, double rotation_window, double pixelSize_ )
{
    
    Array3d pose = laserfan_.pose;
    vector<Array2d> ends = laserfan_.end;
    SearchWindowParameter_.rotation_lower_bound = -1*abs(rotation_window);
    SearchWindowParameter_.translation_lower_bound = -1*abs(linear_window);
    SearchWindowParameter_.rotation_step_num = std::ceil(2*rotation_window/pixelSize_);
    SearchWindowParameter_.translation_step_size = pixelSize_;

    // Solve the step size of rotation
    // Find the longest range from origin to ends
    double max_range = 0;
    vector<Array2d> endsInOriginFrame = BasicFunction::Transform(laserfan_.end,-1*pose);
    for(auto & element : endsInOriginFrame)
    {
        max_range = std::max(max_range,Vector2d(element).norm());
    }
    SearchWindowParameter_.rotation_step_size = std::acos(1.- (pixelSize_)*(pixelSize_)/(2*max_range*max_range));
    SearchWindowParameter_.rotation_step_num = 2*ceil(abs(rotation_window)/SearchWindowParameter_.rotation_step_size);
    
    initializationFlag = true;
}

vector<CSM::DiscreteRotatedResult> CSM::GenerateRotationDiscreteScans(const LaserData& laserdata_,Array3d pose_guess_)
{
    assert(initializationFlag == true);
    vector<CSM::DiscreteRotatedResult> temp;
    double pixelSize_ = SearchWindowParameter_.translation_step_size;
    // 此处从resize 改成了reserve
    temp.reserve(SearchWindowParameter_.rotation_step_num);
    for (int  offset_rotation = 0; offset_rotation < SearchWindowParameter_.rotation_step_num; offset_rotation++)
    {
        // 得到旋转的offset
        Array3d pose_offset = Array3d(0,0,offset_rotation*SearchWindowParameter_.rotation_step_size + SearchWindowParameter_.rotation_lower_bound);
        // 将旋转offset 混合pose_guess 变换到世界坐标系中
        vector<Array2d> worldFrameLaserData = BasicFunction::Transform(laserdata_,pose_guess_+pose_offset);
        // 离散化并记录位姿
        temp.push_back(CSM::DiscreteRotatedResult(BasicFunction::DiscretePoints(worldFrameLaserData,pixelSize_),pose_guess_+pose_offset));
    }
    return temp;
}

CSM::MatchResult CSM::ScoreAllCandidates(const vector<uint16_t>& ProbabilityTable, const vector<Eigen::Array2i>& discretelaserdata,
                                         const MapLimits& maplimits_, const Eigen::AlignedBox2i& scale_,
                                         const vector<Array2i> & candidates)
{
   Array2i offset = Array2i(scale_.min());
   float bestscore = 0.f;
   Array2i BestCandidate(0,0);
   double pixelSize = SearchWindowParameter_.translation_step_size;
   for (auto & candidate: candidates)
   {
        float score = 0.f;
        for (auto & element: discretelaserdata)
        {
            Array2i temp = element -offset + candidate;
            if (maplimits_.Contain(temp))
            {
                score += (1.f - PGMBuilder::ValueToCorrespondenceCost[ProbabilityTable[int(temp(1) * maplimits_.scaleXY(0)) + temp(0)]]);
            }
            else
            continue;
        }
        if (bestscore < score)
        {
            bestscore = score;
            BestCandidate = candidate;
        }
        score = 0;
   }
    return CSM::MatchResult(bestscore,Array3d(BestCandidate(0)*pixelSize,BestCandidate(1)*pixelSize,0));
}
vector<Array2i> CSM::GenerateCandidate()
{
    vector<Array2i> candidates;
    int lower_pixel = floor(SearchWindowParameter_.translation_lower_bound/SearchWindowParameter_.translation_step_size);
    int uppper_pixel = -lower_pixel;
    
    for (int cur_x = lower_pixel ; cur_x <= uppper_pixel;cur_x++ )
    {
        for (int cur_y = lower_pixel ; cur_y <= uppper_pixel;cur_y++  )
        {
            candidates.push_back(Array2i(cur_x,cur_y));
        }
    }
    return candidates;
}

// laserdata_ is representede by the local map
// pose_guess 应该是相对于submap坐标的坐标系
CSM::MatchResult_BBM CSM::BranchAndBoundMatch(shared_ptr<M2MPGM> m2mpgm_, const LaserData& laserdata_, const Array3d& pose_gues)
{
    assert(m2mpgm_->GetFinished() == true);
    double pixelSize = m2mpgm_->GetPixelSize();
    LaserFan laserfan_(laserdata_, pose_gues);
    GenerateSearchWindow(laserfan_, 10e6 * pixelSize, 3.1415926, pixelSize);
    vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans(laserdata_, pose_gues);
    vector<Array4i> limits = ShrinktoFit(m2mpgm_->GetMapLimits(), rotateDiscreteResults, m2mpgm_->GetScale(), m2mpgm_->GetHeight());
    vector<CSM::Candidates_BBM> candidates = GenerateCandidate_BBM(limits, m2mpgm_->GetHeight(), m2mpgm_->GetMapLimits());
    double pixelSize_ = SearchWindowParameter_.translation_step_size;
    CSM::MatchResult BestResult(0.f, Array3d(0, 0, 0));
    
    // 最上层的打分和排序
    assert(candidates.size() > 0);
    vector<CSM::MatchResult_BBM> upper_level_result = ScoreAllCandidates_BBM(m2mpgm_->GetMap(m2mpgm_->GetHeight()),
        rotateDiscreteResults, m2mpgm_->GetMapLimits(), m2mpgm_->GetScale(),
        candidates, m2mpgm_->GetHeight());

    // 分支定界的递归部分
    CSM::MatchResult_BBM result_BBM;
    result_BBM = BranchAndBoundLoop(m2mpgm_, m2mpgm_->GetHeight() - 1, upper_level_result,
        rotateDiscreteResults);
   
    return result_BBM;



    //// TEST PART
    //assert(m2mpgm_->GetFinished() == true);
    //double pixelSize = m2mpgm_->GetPixelSize();
    //LaserFan laserfan_(laserdata_, pose_gues);
    //GenerateSearchWindow(laserfan_, 10e6 * pixelSize, 3.1415926, pixelSize);
    ////分支定界最上层的处理
    //// 先为最上层做分支
    //vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans(laserdata_, pose_gues);
    //// Candidates 对于BBM的划分还要再考虑一下
    //int height = m2mpgm_->GetHeight();
    //// TEST PART
    //// 抽取 349 350
    //vector<CSM::DiscreteRotatedResult> scan_collect;
    //
    //vector<int> index_set;
    //int i = 0;
    //for (auto& element : rotateDiscreteResults)
    //{
    //    if (element.pose(2) < 1.607)
    //    {
    //        if (element.pose(2) > 1.595)
    //        {
    //            scan_collect.push_back(element);
    //            index_set.push_back(i);
    //            Array2i min_, max_;
    //            // 以下主要是为了找到大范围不能用的原因
    //            BasicFunction::FindMaxAndMin(min_, max_, element.DiscreteLaserData);
    //            cout << "min_ is " << min_(0) << "   " << min_(1) << endl;
    //            cout << "max_ is " << max_(0) << "   " << max_(1) << endl;
    //        }
    //    }
    //    i++;
    //}
    //vector<Array4i> limits = ShrinktoFit(m2mpgm_->GetMapLimits(), scan_collect, m2mpgm_->GetScale(),m2mpgm_->GetHeight());
    //vector<CSM::Candidates_BBM> candidates_temp = GenerateCandidate_BBM(limits, height, m2mpgm_->GetMapLimits());
    //vector<CSM::MatchResult_BBM> upper_level_result = ScoreAllCandidates_BBM(m2mpgm_->GetMap(m2mpgm_->GetHeight()),
    //    scan_collect, m2mpgm_->GetMapLimits(), m2mpgm_->GetScale(),
    //    candidates_temp, m2mpgm_->GetHeight());
    //CSM::MatchResult_BBM result_BBM;
    //result_BBM = BranchAndBoundLoop(m2mpgm_, m2mpgm_->GetHeight() - 1, upper_level_result,
    //    scan_collect);
    //return result_BBM;
}

CSM::MatchResult_BBM CSM::BranchAndBoundLoop(shared_ptr<M2MPGM> m2mpgm, int height,
    const vector<CSM::MatchResult_BBM>& upper_level_result, 
    const vector<CSM::DiscreteRotatedResult>& RotDiscreData)
{
    assert(height >= 0);
    vector<MatchResult_BBM> all_choose;
    all_choose.reserve(upper_level_result.size() * 4);
    int stride = 1 << height;
    // 之后修改一下
    vector<Array2i> new_candidates;
    new_candidates.reserve(4);
    for (auto& element : upper_level_result)
    {
        //这里应该有两个score 的判断条件
        if (element.bestscore < cur_maxscore_BBM || element.bestscore < score_threshold)
            continue;
        
        // Generate Candidates of lower height
        // 这里stride 的加减号是由 多层的新时候进行
        for (auto candidate : { Array2i(0,0),Array2i(0,-stride),Array2i(+stride,0),Array2i(+stride,-stride) })
        {
            new_candidates.push_back(candidate + element.candidate);
        }
        vector<CSM::MatchResult_BBM> temp_result = ScoreAllCandidates(m2mpgm->GetMap(height),
            RotDiscreData[element.RotDiscreIndex], m2mpgm->GetMapLimits(), m2mpgm->GetScale(),
            new_candidates,element.RotDiscreIndex,height);
        
        // 基线条件
        if (height == 0)
        {
            if (temp_result.front().bestscore > cur_maxscore_BBM)
            {
                cur_maxscore_BBM = temp_result.front().bestscore;
                return temp_result.front();
            }
            else
            {
                return CSM::MatchResult_BBM();
            }
        }
        all_choose.push_back( BranchAndBoundLoop(m2mpgm, height - 1, temp_result, RotDiscreData));
        new_candidates.clear();
    }
    sort(all_choose.begin(), all_choose.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b)
        {return a.bestscore > b.bestscore; });
    if (all_choose.size() == 0)
        return CSM::MatchResult_BBM();
    return all_choose.front();
}


vector<CSM::MatchResult_BBM> CSM::ScoreAllCandidates_BBM(const vector<float>& ProbabilityTable,
    const vector<CSM::DiscreteRotatedResult>& discretelaserdata, const MapLimits& maplimits_,
    const Eigen::AlignedBox2i& scale_, const  vector<CSM::Candidates_BBM>& candidates,int height)
{
    int stride = 1 << height;
    Array2i offset_stride(stride - 1, 0);
    vector<CSM::MatchResult_BBM> result;
    result.reserve(discretelaserdata.size() * candidates.size());
    int RotDiscreIndex = 0;
    Array2i area;
    int total_area = maplimits_.scaleXY(1) * maplimits_.scaleXY(0);
    // 第一层把 DiscreteRotationResukt 抽出来
    for (const CSM::DiscreteRotatedResult& discreteElement : discretelaserdata)
    {
        Array2i offset = Array2i(scale_.min());
        float bestscore = 0.f;
        Array2i BestCandidate(0, 0);
        double pixelSize = SearchWindowParameter_.translation_step_size;
        for (auto& candidate : candidates[RotDiscreIndex].candidates)
        {
            /*area = max_collector[RotDiscreIndex] + candidate;*/
            //if (area(0) > 0 && area(1) > 0)
            //{
            //    /*if (area(0) * area(1) < (total_area / 10))
            //        continue;*/
            //}
            //else 
            //    continue;
           /* area = min_collector[RotDiscreIndex] + candidate;*/
            //if (area(0) < maplimits_.scaleXY(0) && area(0) < maplimits_.scaleXY(1))
            //{
            ///*    int area_ = (maplimits_.scaleXY(0) - area(0)) * (maplimits_.scaleXY(1) - area(1));
            //    if (area_ < total_area / 10)
            //        continue;*/

            //}
            //else
            //    continue;
            float score = 0.f;
            for (auto& element : discreteElement.DiscreteLaserData)
            {
                Array2i temp = element - offset + candidate + offset_stride;
                if ((temp(0) > (maplimits_.scaleXY(0) + stride - 1 - 1)) || (temp(1) > (maplimits_.scaleXY(1) + stride - 1 - 1)))
                    continue;
                else if (temp(0) < 0 || temp(1) < (0))
                    continue;
                else
                    score += ProbabilityTable[int(temp(1) * (maplimits_.scaleXY(0) + stride - 1)) + temp(0)];
            }
            result.push_back(CSM::MatchResult_BBM(score, discreteElement.pose, candidate, height,RotDiscreIndex));
        }
        RotDiscreIndex++;
    }
    // 降序排列
    sort(result.begin(), result.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b) 
        {return a.bestscore > b.bestscore; });
    return result;
}

vector<CSM::MatchResult_BBM> CSM::ScoreAllCandidates(const vector<float>& ProbabilityTable, 
    const DiscreteRotatedResult& discretelaserdata, const MapLimits& maplimits_,
    const Eigen::AlignedBox2i& scale_, const  vector<Array2i>& candidates,int RotDiscreIndex,int height)
{
    int stride = 1 << height;
    Array2i offset_stride(stride - 1, 0);
    Array2i offset = Array2i(scale_.min());
    float bestscore = 0.f;
    Array2i BestCandidate(0, 0);
    double pixelSize = SearchWindowParameter_.translation_step_size;
    vector<CSM::MatchResult_BBM> result;
    result.reserve(candidates.size());
    for (auto& candidate : candidates)
    {
        float score = 0.f;
        for (auto& element : discretelaserdata.DiscreteLaserData)
        {
            Array2i temp = element - offset + candidate + offset_stride;
            if ((temp(0) > (maplimits_.scaleXY(0) + stride - 1 - 1)) || (temp(1) > (maplimits_.scaleXY(1) + stride - 1 - 1)))
                continue;
            else if (temp(0) < 0 || temp(1) < (0))
                continue;
            else
            score += ProbabilityTable[int(temp(1) * (maplimits_.scaleXY(0)+stride-1)) + temp(0)];
        }
        result.push_back(CSM::MatchResult_BBM(score, discretelaserdata.pose, candidate, 
            height, RotDiscreIndex));   
    }
    sort(result.begin(), result.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b)
    {return a.bestscore > b.bestscore; });
    return result;
}


// 这里的DiscreteRotScans应该是相对于submap的坐标
// 这里钻了2d submap 坐标系相对世界坐标系没有旋转的空子
// 这里没有用到拓展的地图 之后回来确认一下是不是这个的问题 再进行修改
vector<Array4i> CSM::ShrinktoFit(const MapLimits& limits, const vector<DiscreteRotatedResult>& DiscreRotScans
, const Eigen::AlignedBox2i& scale,int height)
{
    int stride_offset = int(1 << height) - 1;
    Array2i offset_map(int(1 << height) - 1, 0);
    Array2i offset = scale.min();
    vector<Array4i> result;
    result.reserve(DiscreRotScans.size());
    for (auto& element : DiscreRotScans)
    {
        Eigen::Array2i min = Array2i::Zero();
        Eigen::Array2i max = Array2i::Zero();
        
        BasicFunction::FindMaxAndMin(min, max, element.DiscreteLaserData);
        // 在这里坐标系进行了更换
        min = min - offset + offset_map;
        max = max - offset + offset_map;

        Array2i scale_ = max - min;
        Array2i min_border = Array2i::Zero();
        min = min.max(min_border);
        Array2i max_border = limits.scaleXY + Array2i(stride_offset,stride_offset) - Array2i(1,1);
        max = max.min(max_border);
        
        // 这里不保险啊 不像是PGM 还可以增大地图 
        // 可能会存在超出地图的情况
        // 那就是必须在子图进行拼接的时候把地图扩大到合适的大小

        //assert(limits.Contain(min));
        //assert(limits.Contain(max));

        // 需要让范围进一步缩小
        // 防止出现巨大偏差的情况
        
        // TODO 进一步缩小范围
        // 排除没有接触的情况
        
        max_collector.push_back(max);
        min_collector.push_back(min);

        result.push_back(Array4i(-min(0), limits.scaleXY(0) + stride_offset - 1 - max(0),
            -min(1), limits.scaleXY(1) + stride_offset - 1 - max(1)));
    }
    return result;
}


vector<CSM::Candidates_BBM> CSM::GenerateCandidate_BBM(const vector<Array4i>& limits_,
    int height,const MapLimits & maplimits_)
{
    int stride = 1 << height;
    vector<CSM::Candidates_BBM> result;
    result.reserve(limits_.size());
    vector<Array2i> candidates;
    // 这个地方之后还是要再调节

    // the lower_pixel must be a negative value;
    int lower_pixel = floor(SearchWindowParameter_.translation_lower_bound / SearchWindowParameter_.translation_step_size);
    // thus the upper_pixel must be a positive value;
    int uppper_pixel = -lower_pixel;
    
    int cur_shift_x = 0;
    int cur_shift_y = 0;

    int index = 0;
    // 这里考虑到在涉及打分函数的时候，会先对点在不在这个Limits 内进行一次判断再决定是否打分
    // 所以可以选择一个比原地图稍大的打分范围，不用担心会不会出线访问地图外的点的问题
    for (auto& limit : limits_)
    {
        /*while (cur_shift_x > limit(0)-stride)
        {
            while (cur_shift_y > limit(2) - stride)
            {
                candidates.push_back(Array2i(cur_shift_x, cur_shift_y));
                cur_shift_y -= stride;
            }
            cur_shift_y = stride;
            while (cur_shift_y < limit(3) + stride)
            {
                candidates.push_back(Array2i(cur_shift_x, cur_shift_y));
                cur_shift_y += stride;
            }
            cur_shift_x -= stride;
        }
        cur_shift_y = 0;
        cur_shift_x = stride;

        while (cur_shift_x < limit(1) + stride)
        {
            while (cur_shift_y > limit(2) - stride)
            {
                candidates.push_back(Array2i(cur_shift_x, cur_shift_y));
                cur_shift_y -= stride;
            }
            cur_shift_y = stride;
            while (cur_shift_y < limit(3) + stride)
            {
                candidates.push_back(Array2i(cur_shift_x, cur_shift_y));
                cur_shift_y += stride;
            }
            cur_shift_x += stride;
        }*/

        Array2i min_ = Array2i(limit(0), limit(2));
        Array2i max_ = Array2i(limit(1), limit(3));
        candidates.reserve(int((max_(1) - min_(0) + stride) / stride) * int((max_(1) - min_(1) + stride) / stride));
        for (int cur_x = min_(0); cur_x < (max_(0) + stride); cur_x += stride)
        {
            for (int cur_y = min_(1); cur_y < (max_(0) + stride); cur_y += stride)
            {
                candidates.push_back(Array2i(cur_x, cur_y));
            }
        }

        result.push_back(CSM::Candidates_BBM(candidates, index));
        index++;
        candidates.clear();
    }
    return result;

}