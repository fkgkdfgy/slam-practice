#pragma once
#include "common.h"
#include "PGM.h"
#include "CSM.h"
class HillClimbFiner
{
public:
	HillClimbFiner(const PGMBuilder& pgm_ref, const CSM & csmer_ref) :pgm_ref_(pgm_ref),
	csmer_ref_(csmer_ref){}
	// laser_data shoule be represented by the local frame
	Array3d FinerEstimaiton(const vector<Array2d>& laser_data_, const Array3d& initial_pose_,float & score);
	void ScoreAll(const vector<Array2d>& laser_data,float & score);
private:
	const PGMBuilder& pgm_ref_;
	const CSM& csmer_ref_;
};

