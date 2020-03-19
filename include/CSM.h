#ifndef _CSM_H_
#define _CSM_H_

#include "common.h"
#include "slam_2d_v01.h"
#include "PGM.h"
#include "M2MPGM.h"
class CSM
{
public:
    struct SearchWindowParameter
    {
        double translation_step_size;
        int translation_step_num;
        double rotation_step_size;
        int rotation_step_num;
        double translation_lower_bound;
        double rotation_lower_bound;
    };
    struct MatchResult
    {
      MatchResult(const float & bestscore_, const Array3d & bestpose_):bestscore(bestscore_),bestpose(bestpose_){}
      MatchResult(){}
      float bestscore ;
      Array3d bestpose;
    };
    struct DiscreteRotatedResult
    {
      DiscreteRotatedResult(){}
      DiscreteRotatedResult(const vector<Array2i> & discretelaserdata, const Array3d pose_):DiscreteLaserData(discretelaserdata),pose(pose_){}
      vector<Array2i> DiscreteLaserData;
      Array3d pose;
    };
    struct MatchResult_BBM:MatchResult
    {
        MatchResult_BBM():MatchResult(0,Array3d(0,0,0)) {}
        MatchResult_BBM(const float& bestscore_, const Array3d& bestpose_,
            const Array2i& candidate_, const int& height_, int RotDiscreIndex_)
            :MatchResult(bestscore_, bestpose_),
            candidate(candidate_), height(height_), RotDiscreIndex(RotDiscreIndex_){}
        Eigen::Array2i candidate;
        int height;
        int RotDiscreIndex;
    };
    struct Candidates_BBM
    {
        Candidates_BBM(const vector<Array2i> & candidates_, 
            int RotDiscreIndex_ ): candidates(candidates_),RotDiscreIndex(RotDiscreIndex_){}
        vector<Eigen::Array2i> candidates;
        int RotDiscreIndex;
    };
    // CSM part
    CSM(){ initializationFlag = false;}
    Eigen::Array3d Match(const LaserData& laserdata_, Eigen::Array3d pose_guess_, const PGMBuilder & PGM_ptr);
    // �õ���DiscreteRotationResult �еĵ㣬���Ѿ������˱任�ĵ�
    void GenerateSearchWindow(const LaserFan & laserfan_, double linear_window, double rotation_window, double pixelSize_);
    void UpdateSearchWindow(double longest_range_);
    vector<CSM::DiscreteRotatedResult> GenerateRotationDiscreteScans(const LaserData& laserdata_, 
        Array3d pose_guess_);
    CSM::MatchResult ScoreAllCandidates(const vector<uint16_t>& ProbabilityTable, const vector<Eigen::Array2i>& discretelaserdata, const MapLimits& maplimits_,
                                        const Eigen::AlignedBox2i& scale_,const  vector<Array2i> & candidates);
    vector<Array2i> GenerateCandidate();
    
    // Branch and Bound Part
    CSM::MatchResult_BBM BranchAndBoundMatch(shared_ptr<M2MPGM>, const LaserData& laserdata_, const Array3d& pose_guess);
    
    CSM::MatchResult_BBM BranchAndBoundLoop(shared_ptr<M2MPGM> m2mpgm, int height,
       const vector<CSM::MatchResult_BBM> & , const vector<CSM::DiscreteRotatedResult> &);

    // BBM ��������ֺ���
    vector<CSM::MatchResult_BBM> ScoreAllCandidates_BBM(const vector<float>& ProbabilityTable,
        const vector<CSM::DiscreteRotatedResult>& discretelaserdata, const MapLimits& maplimits_,
        const Eigen::AlignedBox2i& scale_, const  vector<CSM::Candidates_BBM>& candidates, int height);
    
    vector<CSM::MatchResult_BBM> ScoreAllCandidates(const vector<float>& ProbabilityTable, const DiscreteRotatedResult& discretelaserdata, const MapLimits& maplimits_,
        const Eigen::AlignedBox2i& scale_, const  vector<Array2i>& candidates, int RotDiscreIndex, int height);

    vector<Candidates_BBM> GenerateCandidate_BBM(const vector<Array4i> & limits, int height ,
        const MapLimits & );
    

    // ������������С�͵�ͼ��Сƥ�䵽һ�����ʵĳߴ�
    // ��ֹ������Χ���󣬵��������ϵ��˷�
    // ȥ�������candidate;
    // �����漰��candidate��ɾ����vector�Ѿ���������
    // ���Ծ�ֱ�Ӳ���candiddate ���Ǹĳ���bound ����candidate ��������
    // Array4i �洢˳�� ---> X-left X-right Y-down Y-up ��˼�� discan�ܽ���ƽ�Ƶķ�Χ left ֵ�Ǹ�ֵ right ֵ����ֵ
    vector<Array4i> ShrinktoFit(const MapLimits& limits, const vector<DiscreteRotatedResult>& DiscreRotScan
        ,const Eigen::AlignedBox2i & scale,int height);
private:
    SearchWindowParameter SearchWindowParameter_;
    bool initializationFlag;
    float BBMscore;
    float score_threshold;
    float cur_maxscore_BBM;
    vector<Array2i> max_collector;
    vector<Array2i> min_collector;
};


#endif //