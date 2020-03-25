#pragma once
#include "common.h"
#include "ceres/ceres.h"
#include "PGM.h"
#include "ceres/cubic_interpolation.h"

// make the PGMBuilder can be used by the Ceres Frame
// to fit the BicubicInterpolation
// we made an adapter for the PGM

class GridAdapter
{
public:
    enum { DATA_DIMENSION = 1 };
    explicit GridAdapter(const PGMBuilder & grid) : grid_(grid) {}

    // Here ,we assume the coordinate (input) is represented by the local frame
    void GetValue(const int row,const int col, double* const value) const
    {
        if(!grid_.maplimits_.Contain(Array2i(col,row)))
            *value = 0.9;
        else
            *value = double(1) - grid_.GetProbability_LocalCandidates(Array2i(col,row));
    }
private:
    const PGMBuilder & grid_;
};

// Error Function
class OccupiedSpaceError
{
public:
    // Here, Currently, we do not use the scale_factor
    OccupiedSpaceError(const vector<Eigen::Array2d> & Laserdata_ptr,
            const double & scale_factor,
            const PGMBuilder & pgm_ptr_):laserdata_(Laserdata_ptr),pgm_ptr(pgm_ptr_),
            scale_factor_(scale_factor) {}
    static ceres::CostFunction* CreateAutoDiffCostFunction(const vector<Eigen::Array2d> & Laserdata_ptr,
                                                           const double scale_factor,const PGMBuilder & pgm_ptr_)
    {
        return new ceres::AutoDiffCostFunction<OccupiedSpaceError, ceres::DYNAMIC, 3>
                (
                        new OccupiedSpaceError(Laserdata_ptr, scale_factor, pgm_ptr_),
                        Laserdata_ptr.size()
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
        const   GridAdapter adapter(pgm_ptr);
        ceres::BiCubicInterpolator<GridAdapter> interpolator(adapter);
        for (auto& point : laserdata_)
        {
            const Matrix<T, 3, 1> point_(T(point.x()), T(point.y()), T(1));
            const Matrix<T, 3, 1> point_world_ = translation_matrix_ * point_;
            // the coordinate here we used is in the localmap frame
            // 这里确实是有问题
            interpolator.Evaluate(point_world_(0)/T(pgm_ptr.pixelSize_)  - T(pgm_ptr.scale_.min()(0)) +T(0.5) ,
                                  point_world_(1)/T(pgm_ptr.pixelSize_)  - T(pgm_ptr.scale_.min()(1))+ T(0.5) ,&residual[i]);
            i++;
        }
        return true;
    }
private:

    const PGMBuilder & pgm_ptr;
    double scale_factor_;
    const vector<Eigen::Array2d> & laserdata_;
};

class OptimizationProblem
{
public:
	explicit OptimizationProblem(double scale_factor,
	        const vector<Eigen::Array2d> & Laserdata_ptr,
	const PGMBuilder & pgm_ptr, const Eigen::Array3d initial_pose):
	scale_factor_(scale_factor),Laserdata_ptr_(Laserdata_ptr),pgm_ptr_(pgm_ptr)
	,initial_pose_(initial_pose) {}
	void Solve();
	const Eigen::Array3d & GetResult() const  {return result_pose;}
private:
	double scale_factor_;
	const vector<Eigen::Array2d> & Laserdata_ptr_;
	const PGMBuilder & pgm_ptr_;
	Eigen::Array3d initial_pose_;
	Eigen::Array3d result_pose;
};
