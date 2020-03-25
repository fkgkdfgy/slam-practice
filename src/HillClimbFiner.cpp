#include "HillClimbFiner.h"


Array3d HillClimbFiner::FinerEstimaiton(const vector<Array2d>& laser_data_, const Array3d& initial_pose_, float & score)
{
    // transfer the scan points from the scan frame into world frame
    double ipixel = 1 / pgm_ref_.pixelSize_;

    // the main part of HillClimb
    double trans_ = csmer_ref_.GetWindowParameter().translation_step_size*0.8;
    double rot_angle = csmer_ref_.GetWindowParameter().rotation_step_size*0.8;
    Array3d trans_range(-1 * trans_, 0, trans_);
    Array3d rot_range(-1 * rot_angle, 0, rot_angle);
    unsigned int IterMax = 50;
    unsigned int depthMax = 4;
    unsigned int Iter = 0;
    unsigned int depth = 0;
    float bestscore = -1;
    bool noChange = 1;

    MapBuilder::Pose dBestPose(0, 0, 0);
    MapBuilder::Pose middlePose = initial_pose_;
    vector<int> BestScoreIndex;
    for (; Iter < IterMax; Iter++)
    {
        noChange = true;
        for (unsigned int angle_ = 0; angle_ < 3; angle_++)
        {
            Array3d temp_pose = Array3d(0, 0, rot_range(angle_)) + middlePose;
            MapBuilder::LaserData scan_tem_rot = BasicFunction::Transform(laser_data_, temp_pose);
            for (unsigned int trans_x = 0; trans_x < 3; trans_x++)
            {
                for (unsigned int trans_y = 0; trans_y < 3; trans_y++)
                {
                    // 变换的部分
                    Array3d temp_pose = Array3d(trans_range(trans_x), trans_range(trans_y), rot_range(angle_)) + middlePose;
                    MapBuilder::LaserData scan_temp= BasicFunction::Transform(laser_data_, temp_pose);
                    // 打分系统 统计未知的点
                    float score;
                    ScoreAll(scan_temp, score);
                    if (score > bestscore)
                    {
                        bestscore = score;
                        dBestPose = temp_pose;
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
       cout << Iter << endl;
       score = bestscore;
       return dBestPose;
}

void HillClimbFiner::ScoreAll(const vector<Array2d>& laser_data_ ,float & score)
{
    // 离散化
    score = 0.f;
    int half_count = 0;
    vector<Array2i> discre_data_ = BasicFunction::DiscretePoints(laser_data_, pgm_ref_.pixelSize_);
    for (auto & element:discre_data_)
    { 
        element = element - Array2i(pgm_ref_.scale_.min());
        if (pgm_ref_.GetMapLimits().Contain(element))
        {
            float temp_score = pgm_ref_.GetProbability_LocalCandidates(element);
            if (temp_score != float(0.5))
            {
                score += temp_score;
                half_count++;
            }
        }
        continue;
    }
    
    if (half_count <= laser_data_.size() *0.45)
        score = 0.f;
    else
        score = score / half_count;

}