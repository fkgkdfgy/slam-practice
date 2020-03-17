#pragma once
#include "ceres/ceres.h"
#include "PGM.h"

class OptimizationProblem
{
public:
	OptimizationProblem(double scale_factor,
	shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr,
	shared_ptr<PGMBuilder> pgm_ptr, const Eigen::Array3d initial_pose): 
	scale_factor_(scale_factor),Laserdata_ptr_(Laserdata_ptr),pgm_ptr_(pgm_ptr)
	,initial_pose_(initial_pose)
	{}
	
	void Solve();
private:
	double scale_factor_;
	shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr_;
	shared_ptr<PGMBuilder> pgm_ptr_;
	Eigen::Array3d initial_pose_;
};

class OccupiedSpaceError
{
public:
	OccupiedSpaceError(shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr, const double & scale_factor,
		shared_ptr<PGMBuilder> pgm_ptr_):laserdata_(Laserdata_ptr),pgm_ptr(pgm_ptr_),scale_factor_(scale_factor)
	{
		interpolator = make_shared<BicubicInerpolator>(new BicubicInerpolator());
	}
	static ceres::CostFunction* CreateAutoDiffCostFunction(shared_ptr<vector<Eigen::Array2d>> Laserdata_ptr,
		const double scale_factor,shared_ptr<PGMBuilder> pgm_ptr_)
	{
		return new ceres::AutoDiffCostFunction<OccupiedSpaceError, ceres::DYNAMIC, 3>
			(
				new OccupiedSpaceError(Laserdata_ptr, scale_factor, pgm_ptr_)
			);
	}
	template <typename T>
	bool operator ()(const T* const pose, T * residual) 
	{
		Eigen::Rotation2D<T> rotation_(pose[2]);
		Eigen::Matrix<T,2,2> rotation_matrix = rotation_.toRotationMatrix();
		Eigen::Matrix<T,2,1> translation_(pose[0], pose[1]);
		Eigen::Matrix<T, 3, 3> translation_matrix_;
		translation_matrix_ << rotation_matrix,translation_matrix,T(0),T(0),T(1);
		int i = 0;
		for (auto& point : *laserdata_)
		{
			const Matrix<T, 3, 1> point_(T(point.x()), T(point.y()), T(1));
			const Matrix<T, 3, 1> point_world_ = translation_matrix_ * point_;
			interpolator->Evaluate(point_world_, residual[i],NULL);
			i++;
		}
	}
private:
	shared_ptr<BicubicInerpolator> interpolator;
	shared_ptr<PGMBuilder> pgm_ptr;
	double scale_factor_;
	shared_ptr<vector<Eigen::Array2d>> laserdata_;
};

class BicubicInerpolator
{
public:
	BicubicInerpolator() {}
	template <typename T,typename T1>
	void Evaluate(const T & ,T1 * f_, T1 * df_) 
	{
		const int x = std::round(T(0));
		const int y = std::round(T(1));
		double p0, p1, p2, p3, dx, f0,f1,f2,f3,f;
		double df0, df1, df2, df3, df;
		
		p0 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x - 1, y-1)));
		p1 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x , y-1)));
		p2 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 1, y-1)));
		p3 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 2, y-1)));
		CubicInterpolation<double>(p0, p1, p2, p3, T(0) - x + 0.5, f0, df0);
		
		p0 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x - 1, y )));
		p1 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x, y )));
		p2 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 1, y )));
		p3 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 2, y )));
		CubicInterpolation<double>(p0, p1, p2, p3, T(0) - x + 0.5, f1, df1);

		p0 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x - 1, y+1)));
		p1 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x, y+1)));
		p2 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 1, y+1)));
		p3 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 2, y+1)));
		CubicInterpolation<double>(p0, p1, p2, p3, T(0) - x + 0.5, f2, df2);

		p0 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x - 1, y + 2)));
		p1 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x, y + 2)));
		p2 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 1, y + 2)));
		p3 = double(1.f-pgm_ptr->GetProbability_WorldCandidates(Eigen::Array2i(x + 2, y + 2)));
		CubicInterpolation<double>(p0, p1, p2, p3, T(0) - x + 0.5, f3, df3);

		CubicInterpolation<double>(f0, f1, f2, f3, T(1) - y + 0.5, &f_, &df_);
		if (dfdc != NULL) 
		{
			// Interpolate vertically the derivative along the columns.
			CubicHermiteSpline<Grid::DATA_DIMENSION>(df0, df1, df2, df3,
				T(1) - y + 0.5, &df_, NULL);
		}
	}
	template<typename T>
	void CubicInterpolation(const T& t1, const T& t2, const T& t3, 
	const T& t4, const T & td, T * estimation,T* estimation_derivation);
private:
	shared_ptr<PGMBuilder> pgm_ptr;
};