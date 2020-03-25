#include "CSM.h"

float CSM::BBM_min_score = 0.6;
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
    // ��ʵ΢С�˶��͹�ֵ�෴
    cout << " the best average score is " << endl;
    cout << BestResult.bestscore << endl;
    cout << " the matched pose is " << endl;
    cout << BestResult.bestpose(0) << "  " << BestResult.bestpose(1) << "  " << BestResult.bestpose(2) << endl;
    return BestResult.bestpose;
}
void CSM::GenerateSearchWindow(const LaserFan& laserfan_, double linear_window, double rotation_window, double pixelSize_)
{

    Array3d pose = laserfan_.pose;
    vector<Array2d> ends = laserfan_.end;
    SearchWindowParameter_.rotation_lower_bound = -1 * abs(rotation_window);
    SearchWindowParameter_.translation_lower_bound = -1 * abs(linear_window);
    SearchWindowParameter_.rotation_step_num = std::ceil(2 * rotation_window / pixelSize_);
    SearchWindowParameter_.translation_step_size = pixelSize_;

    // Solve the step size of rotation
    // Find the longest range from origin to ends
    double max_range = 0;
    vector<Array2d> endsInOriginFrame = BasicFunction::Transform(laserfan_.end, -1 * pose);
    for (auto& element : endsInOriginFrame)
    {
        max_range = std::max(max_range, Vector2d(element).norm());
    }
    SearchWindowParameter_.rotation_step_size = std::acos(1. - (pixelSize_) * (pixelSize_) / (2 * max_range * max_range));
    SearchWindowParameter_.rotation_step_num = 2 * ceil(abs(rotation_window) / SearchWindowParameter_.rotation_step_size);

    initializationFlag = true;
}
void CSM::GenerateSearchWindow_BBM(const LaserFan & laserfan_, double linear_window, double rotation_window, double pixelSize_ )
{
    
    Array3d pose = laserfan_.pose;
    vector<Array2d> ends = laserfan_.end;
    SearchWindowParameter_BBM.rotation_lower_bound = -1*abs(rotation_window);
    SearchWindowParameter_BBM.translation_lower_bound = -1*abs(linear_window);
    SearchWindowParameter_BBM.rotation_step_num = std::ceil(2*rotation_window/pixelSize_);
    SearchWindowParameter_BBM.translation_step_size = pixelSize_;

    // Solve the step size of rotation
    // Find the longest range from origin to ends
    double max_range = 0;
    vector<Array2d> endsInOriginFrame = BasicFunction::Transform(laserfan_.end,-1*pose);
    for(auto & element : endsInOriginFrame)
    {
        max_range = std::max(max_range,Vector2d(element).norm());
    }
    SearchWindowParameter_BBM.rotation_step_size = std::acos(1.- (pixelSize_)*(pixelSize_)/(2*max_range*max_range));
    SearchWindowParameter_BBM.rotation_step_num = 2*ceil(abs(rotation_window)/SearchWindowParameter_BBM.rotation_step_size);
    
    initializationFlag_BBM = true;
}

vector<CSM::DiscreteRotatedResult> CSM::GenerateRotationDiscreteScans(const LaserData& laserdata_,Array3d pose_guess_)
{
    assert(initializationFlag == true);
    vector<CSM::DiscreteRotatedResult> temp;
    double pixelSize_ = SearchWindowParameter_.translation_step_size;
    // �˴���resize �ĳ���reserve
    temp.reserve(SearchWindowParameter_.rotation_step_num);
    for (int  offset_rotation = 0; offset_rotation < SearchWindowParameter_.rotation_step_num; offset_rotation++)
    {
        // �õ���ת��offset
        Array3d pose_offset = Array3d(0,0,offset_rotation*SearchWindowParameter_.rotation_step_size + SearchWindowParameter_.rotation_lower_bound);
        // ����תoffset ���pose_guess �任����������ϵ��
        vector<Array2d> worldFrameLaserData = BasicFunction::Transform(laserdata_,pose_guess_+pose_offset);
        // ��ɢ������¼λ��
        temp.push_back(CSM::DiscreteRotatedResult(BasicFunction::DiscretePoints(worldFrameLaserData,pixelSize_),pose_guess_+pose_offset));
    }
    return temp;
}
vector<CSM::DiscreteRotatedResult> CSM::GenerateRotationDiscreteScans_BBM(const LaserData& laserdata_, Array3d pose_guess_)
{
    assert(initializationFlag == true);
    vector<CSM::DiscreteRotatedResult> temp;
    double pixelSize_ = SearchWindowParameter_BBM.translation_step_size;
    // �˴���resize �ĳ���reserve
    temp.reserve(SearchWindowParameter_BBM.rotation_step_num);
    for (int offset_rotation = 0; offset_rotation < SearchWindowParameter_BBM.rotation_step_num; offset_rotation++)
    {
        // �õ���ת��offset
        Array3d pose_offset = Array3d(0, 0, offset_rotation * SearchWindowParameter_BBM.rotation_step_size + SearchWindowParameter_BBM.rotation_lower_bound);
        // ����תoffset ���pose_guess �任����������ϵ��
        vector<Array2d> worldFrameLaserData = BasicFunction::Transform(laserdata_, pose_guess_ + pose_offset);
        // ��ɢ������¼λ��
        temp.push_back(CSM::DiscreteRotatedResult(BasicFunction::DiscretePoints(worldFrameLaserData, pixelSize_), pose_guess_ + pose_offset));
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
   int border = discretelaserdata.size()*0.8;
   for (auto & candidate: candidates)
   {
        float score = 0.f;
        int count = 0;
        for (auto & element: discretelaserdata)
        {
            Array2i temp = element -offset + candidate;
            if (maplimits_.Contain(temp))
            {
                score += (1.f - PGMBuilder::ValueToCorrespondenceCost[ProbabilityTable[int(temp(1) * maplimits_.scaleXY(0)) + temp(0)]]);
                count++;
            }
            else
            continue;
        }
        
        if (count < border)
            continue;

        score = score / count;
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
// pose_guess Ӧ���������submap���������ϵ
CSM::MatchResult_BBM CSM::BranchAndBoundMatch(shared_ptr<M2MPGM> m2mpgm_, const LaserData& laserdata_, const Array3d& pose_gues)
{
    assert(m2mpgm_->GetFinished() == true);
    double pixelSize = m2mpgm_->GetPixelSize();
    LaserFan laserfan_(laserdata_, pose_gues);
    GenerateSearchWindow_BBM(laserfan_, 10e6 * pixelSize, 3.1415926, pixelSize);
    vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans_BBM(laserdata_, pose_gues);
    vector<Array4i> limits = ShrinktoFit(m2mpgm_->GetMapLimits(), rotateDiscreteResults, m2mpgm_->GetScale(), m2mpgm_->GetHeight());
    vector<CSM::Candidates_BBM> candidates = GenerateCandidate_BBM(limits, m2mpgm_->GetHeight(), m2mpgm_->GetMapLimits());
    double pixelSize_ = SearchWindowParameter_BBM.translation_step_size;
    CSM::MatchResult BestResult(0.f, Array3d(0, 0, 0));
    
    // ���ϲ�Ĵ�ֺ�����
    assert(candidates.size() > 0);
    vector<CSM::MatchResult_BBM> upper_level_result = ScoreAllCandidates_BBM(m2mpgm_->GetMap(m2mpgm_->GetHeight()),
        rotateDiscreteResults, m2mpgm_->GetMapLimits(), m2mpgm_->GetScale(),
        candidates, m2mpgm_->GetHeight());

    // ��֧����ĵݹ鲿��
    CSM::MatchResult_BBM result_BBM;
    result_BBM = BranchAndBoundLoop(m2mpgm_, m2mpgm_->GetHeight() - 1, upper_level_result,
        rotateDiscreteResults);
   
    return result_BBM;

    //// TEST PART
    //assert(m2mpgm_->GetFinished() == true);
    //double pixelSize = m2mpgm_->GetPixelSize();
    //LaserFan laserfan_(laserdata_, pose_gues);
    //GenerateSearchWindow(laserfan_, 10e6 * pixelSize, 3.1415926, pixelSize);
    ////��֧�������ϲ�Ĵ���
    //// ��Ϊ���ϲ�����֧
    //vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans(laserdata_, pose_gues);
    //// Candidates ����BBM�Ļ��ֻ�Ҫ�ٿ���һ��
    //int height = m2mpgm_->GetHeight();
    //// TEST PART
    //// ��ȡ 349 350
    //vector<CSM::DiscreteRotatedResult> scan_collect;
    //
    //vector<int> index_set;
    //int i = 0;
    //for (auto& element : rotateDiscreteResults)
    //{
    //    if (element.pose(2) < 2.006)
    //    {
    //        if (element.pose(2) > 1.996)
    //        {
    //            scan_collect.push_back(element);
    //            index_set.push_back(i);
    //            Array2i min_, max_;
    //            // ������Ҫ��Ϊ���ҵ���Χ�����õ�ԭ��
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

CSM::MatchResult_BBM CSM::BranchAndBoundMatch(M2MPGM & m2mpgm_, const LaserData& laserdata_, const Array3d& pose_gues)
{
    assert(m2mpgm_.GetFinished() == true);
    double pixelSize = m2mpgm_.GetPixelSize();
    LaserFan laserfan_(laserdata_, pose_gues);
    GenerateSearchWindow_BBM(laserfan_, 10e6 * pixelSize, 3.1415926, pixelSize);
    vector<CSM::DiscreteRotatedResult> rotateDiscreteResults = GenerateRotationDiscreteScans_BBM(laserdata_, pose_gues);
    vector<Array4i> limits = ShrinktoFit(m2mpgm_.GetMapLimits(), rotateDiscreteResults, 
        m2mpgm_.GetScale(), m2mpgm_.GetHeight());
    vector<CSM::Candidates_BBM> candidates = GenerateCandidate_BBM(limits, m2mpgm_.GetHeight(), m2mpgm_.GetMapLimits());
    double pixelSize_ = SearchWindowParameter_.translation_step_size;
    CSM::MatchResult BestResult(0.f, Array3d(0, 0, 0));

    // ���ϲ�Ĵ�ֺ�����
    assert(candidates.size() > 0);
    vector<CSM::MatchResult_BBM> upper_level_result = ScoreAllCandidates_BBM(m2mpgm_.GetMap(m2mpgm_.GetHeight()),
        rotateDiscreteResults, m2mpgm_.GetMapLimits(), m2mpgm_.GetScale(),
        candidates, m2mpgm_.GetHeight());

    // ��֧����ĵݹ鲿��
    CSM::MatchResult_BBM result_BBM;
    result_BBM = BranchAndBoundLoop(m2mpgm_, m2mpgm_.GetHeight() - 1, upper_level_result,
        rotateDiscreteResults);

    return result_BBM;
}

CSM::MatchResult_BBM CSM::BranchAndBoundLoop(M2MPGM & m2mpgm, int height,
    const vector<CSM::MatchResult_BBM>& upper_level_result,
    const vector<CSM::DiscreteRotatedResult>& RotDiscreData)
{
    assert(height >= 0);
    vector<MatchResult_BBM> all_choose;
    all_choose.reserve(upper_level_result.size() * 4);
    int stride = 1 << height;
    // ֮���޸�һ��
    vector<Array2i> new_candidates;
    new_candidates.reserve(4);
    for (auto& element : upper_level_result)
    {
        //����Ӧ��������score ���ж�����
        if (element.bestscore < cur_maxscore_BBM || element.bestscore < score_threshold)
            continue;

        // Generate Candidates of lower height
        // ����stride �ļӼ������� ������ʱ�����
        for (auto candidate : { Array2i(0,0),Array2i(0,-stride),Array2i(+stride,0),Array2i(+stride,-stride) })
        {
            new_candidates.push_back(candidate + element.candidate);
        }
        vector<CSM::MatchResult_BBM> temp_result = ScoreAllCandidates(m2mpgm.GetMap(height),
            RotDiscreData[element.RotDiscreIndex], m2mpgm.GetMapLimits(), m2mpgm.GetScale(),
            new_candidates, element.RotDiscreIndex, height);
        if (temp_result.empty())
            continue;
        // ��������
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
        all_choose.push_back(BranchAndBoundLoop(m2mpgm, height - 1, temp_result, RotDiscreData));
        new_candidates.clear();
    }
    std::sort(all_choose.begin(), all_choose.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b)
        {return a.bestscore > b.bestscore; });
    if (all_choose.size() == 0)
        return CSM::MatchResult_BBM();
    return all_choose.front();
}


CSM::MatchResult_BBM CSM::BranchAndBoundLoop(shared_ptr<M2MPGM> m2mpgm, int height,
    const vector<CSM::MatchResult_BBM>& upper_level_result, 
    const vector<CSM::DiscreteRotatedResult>& RotDiscreData)
{
    assert(height >= 0);
    vector<MatchResult_BBM> all_choose;
    all_choose.reserve(upper_level_result.size() * 4);
    int stride = 1 << height;
    // ֮���޸�һ��
    vector<Array2i> new_candidates;
    new_candidates.reserve(4);
    for (auto& element : upper_level_result)
    {
        //����Ӧ��������score ���ж�����
        if (element.bestscore < cur_maxscore_BBM || element.bestscore < score_threshold)
            continue;
        
        // Generate Candidates of lower height
        // ����stride �ļӼ������� ������ʱ�����
        for (auto candidate : { Array2i(0,0),Array2i(0,-stride),Array2i(+stride,0),Array2i(+stride,-stride) })
        {
            new_candidates.push_back(candidate + element.candidate);
        }
        vector<CSM::MatchResult_BBM> temp_result = ScoreAllCandidates(m2mpgm->GetMap(height),
            RotDiscreData[element.RotDiscreIndex], m2mpgm->GetMapLimits(), m2mpgm->GetScale(),
            new_candidates,element.RotDiscreIndex,height);
        if (temp_result.empty())
            continue;
        // ��������
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
    std::sort(all_choose.begin(), all_choose.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b)
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
    // ��һ��� DiscreteRotationResukt �����
    for (const CSM::DiscreteRotatedResult& discreteElement : discretelaserdata)
    {
        Array2i offset = Array2i(scale_.min());
        float bestscore = 0.f;
        Array2i BestCandidate(0, 0);
        double pixelSize = SearchWindowParameter_.translation_step_size;
        for (auto& candidate : candidates[RotDiscreIndex].candidates)
        {
            int hit_count = 0;
            float score = 0.f;
            for (auto& element : discreteElement.DiscreteLaserData)
            {
                Array2i temp = element - offset + candidate + offset_stride;
                if ((temp(0) > (maplimits_.scaleXY(0) + stride - 1 - 1)) || (temp(1) > (maplimits_.scaleXY(1) + stride - 1 - 1)))
                    continue;
                else if (temp(0) < 0 || temp(1) < (0))
                    continue;
                else
                {
                    float score_temp = ProbabilityTable[int(temp(1) * (maplimits_.scaleXY(0) + stride - 1)) + temp(0)];
                    if (score_temp == 0.5)
                        continue;
                    else
                    {
                        score += score_temp;
                        hit_count++;
                    }
                }
            }
            score = score / hit_count;
            if (hit_count > (discreteElement.DiscreteLaserData.size() * 0.47)&& score>BBM_min_score)
            {
                result.push_back(CSM::MatchResult_BBM(score, discreteElement.pose, candidate, height, RotDiscreIndex));
            }
        }
        RotDiscreIndex++;
    }
    // ��������
    if (result.empty())
        return result;
    std::sort(result.begin(), result.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b) 
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
        int hit_count = 0;
        for (auto& element : discretelaserdata.DiscreteLaserData)
        {
            Array2i temp = element - offset + candidate + offset_stride;
            if ((temp(0) > (maplimits_.scaleXY(0) + stride - 1 - 1)) || (temp(1) > (maplimits_.scaleXY(1) + stride - 1 - 1)))
                continue;
            else if (temp(0) < 0 || temp(1) < (0))
                continue;
            else
            {
                float score_temp = ProbabilityTable[int(temp(1) * (maplimits_.scaleXY(0) + stride - 1)) + temp(0)];
                if (score_temp == 0.5)
                    continue;
                else
                {
                    score += score_temp;
                    hit_count++;
                }
            }
        }
        score = score / hit_count;
        if (hit_count > (discretelaserdata.DiscreteLaserData.size() * 0.47) && score > BBM_min_score)
        {    
            result.push_back(CSM::MatchResult_BBM(score, discretelaserdata.pose, candidate,
                height, RotDiscreIndex));
        }
    }
    if (result.empty())
        return result;
    sort(result.begin(), result.end(), [](CSM::MatchResult_BBM& a, CSM::MatchResult_BBM& b)
    {return a.bestscore > b.bestscore; });
    return result;
}


// �����DiscreteRotScansӦ���������submap������
// ��������2d submap ����ϵ�����������ϵû����ת�Ŀ���
// ����û���õ���չ�ĵ�ͼ ֮�����ȷ��һ���ǲ������������ �ٽ����޸�
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
        // ����������ϵ�����˸���
        min = min - offset + offset_map;
        max = max - offset + offset_map;

        Array2i scale_ = max - min;
        Array2i min_border = Array2i::Zero();
        min = min.max(min_border);
        Array2i max_border = limits.scaleXY + Array2i(stride_offset,stride_offset) - Array2i(1,1);
        max = max.min(max_border);
        
        // ���ﲻ���հ� ������PGM �����������ͼ 
        // ���ܻ���ڳ�����ͼ�����
        // �Ǿ��Ǳ�������ͼ����ƴ�ӵ�ʱ��ѵ�ͼ���󵽺��ʵĴ�С

        //assert(limits.Contain(min));
        //assert(limits.Contain(max));

        // ��Ҫ�÷�Χ��һ����С
        // ��ֹ���־޴�ƫ������
        
        // TODO ��һ����С��Χ
        // �ų�û�нӴ������
        
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
    // the lower_pixel must be a negative value;
    int lower_pixel = floor(SearchWindowParameter_BBM.translation_lower_bound / SearchWindowParameter_BBM.translation_step_size);
    // thus the upper_pixel must be a positive value;
    int uppper_pixel = -lower_pixel;
    
    int cur_shift_x = 0;
    int cur_shift_y = 0;

    int index = 0;
    // ���￼�ǵ����漰��ֺ�����ʱ�򣬻��ȶԵ��ڲ������Limits �ڽ���һ���ж��پ����Ƿ���
    // ���Կ���ѡ��һ����ԭ��ͼ�Դ�Ĵ�ַ�Χ�����õ��Ļ᲻����߷��ʵ�ͼ��ĵ������
    for (auto& limit : limits_)
    {


        Array2i min_ = Array2i(limit(0), limit(2));
        Array2i max_ = Array2i(limit(1), limit(3));
        candidates.reserve(int((max_(1) - min_(0) + stride) / stride) * int((max_(1) - min_(1) + stride) / stride));
        for (int cur_x = min_(0); cur_x < (max_(0) + stride); cur_x += stride)
        {
            for (int cur_y = min_(1); cur_y < (max_(1) + stride); cur_y += stride)
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