#pragma once
#include "ceres/ceres.h"
#include "PGM.h"
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
class GridAdapter
{
public:
    enum { DATA_DIMENSION = 1 };
    explicit GridAdapter(const PGMBuilder & grid) : grid_(grid) {}
    void GetValue(const int x,const int y, double* const value) const
    {
        if(!grid_.maplimits_.Contain(Array2i(x,y)))
            *value = 0.9;
        else
            *value = double(1) - grid_.GetProbability_LocalCandidates(Array2i(x,y));
    }
private:
    const PGMBuilder & grid_;
};
// Error Function
class OccupiedSpaceError
{
public:
    OccupiedSpaceError(shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr, const double & scale_factor,
                       shared_ptr<PGMBuilder> pgm_ptr_):laserdata_(Laserdata_ptr),pgm_ptr(pgm_ptr_),scale_factor_(scale_factor) {}
    static ceres::CostFunction* CreateAutoDiffCostFunction(shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr,
                                                           const double scale_factor,shared_ptr<PGMBuilder> pgm_ptr_)
    {
        return new ceres::AutoDiffCostFunction<OccupiedSpaceError, ceres::DYNAMIC, 3>
                (
                        new OccupiedSpaceError(Laserdata_ptr, scale_factor, pgm_ptr_)
                );
    }

    template <typename T>
    bool operator () (const T* const pose, T * residual) const
    {

        Eigen::Rotation2D<T> rotation_(pose[2]);
        Eigen::Matrix<T,2,2> rotation_matrix = rotation_.toRotationMatrix();
        Eigen::Matrix<T,2,1> translation_(pose[0], pose[1]);
        Eigen::Matrix<T, 3, 3> translation_matrix_;
        translation_matrix_ << rotation_matrix,translation_,T(0),T(0),T(1);
        int i = 0;
        const   GridAdapter adapter(*pgm_ptr);
        ceres::BiCubicInterpolator<GridAdapter> interpolator(adapter);
        for (auto& point : *laserdata_)
        {
            const Matrix<T, 3, 1> point_(T(point.x()), T(point.y()), T(1));
            const Matrix<T, 3, 1> point_world_ = translation_matrix_ * point_;
            interpolator.Evaluate(point_world_(0) - T(pgm_ptr->scale_.min()(0)) +T(0.5) ,
                    point_world_(1)- T(pgm_ptr->scale_.min()(1)) + T(0.5) ,&residual[i]);
            i++;
        }
        return true;
    }
private:

    shared_ptr<PGMBuilder> pgm_ptr;
    double scale_factor_;
    shared_ptr<vector<Eigen::Array2d>> laserdata_;
};

class OptimizationProblem
{
public:
	OptimizationProblem(double scale_factor,
	shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr,
	shared_ptr<PGMBuilder> pgm_ptr, const Eigen::Array3d initial_pose): 
	scale_factor_(scale_factor),Laserdata_ptr_(Laserdata_ptr),pgm_ptr_(pgm_ptr)
	,initial_pose_(initial_pose) {}
	void Solve();
private:
	double scale_factor_;
	shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr_;
	shared_ptr<PGMBuilder> pgm_ptr_;
	Eigen::Array3d initial_pose_;
};
