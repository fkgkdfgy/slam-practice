#ifndef _TEST_
#define _TEST_
#include "common.h"
#include "PGM.h"
#include "slam_2d_v01.h"
#include "CSM.h"
#include "M2MPGM.h"
#include "HillClimbFiner.h"
void TEST_PGM_INITIALIZATION();
void TEST_PGM_INSERTSCAN();
void TEST_CSM();
void TEST_M2MPGM();
void TEST_BBM_SIMPLE(Array3d target_);
void TEST_CERES();
void TEST_KDtree();
void TEST_SLAM_V01();
void TEST_SLAM_V02();
void TEST_FINERMATCHER();
#endif