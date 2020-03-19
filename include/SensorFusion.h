#pragma once
#include "common.h"
class SensorFusion;
class imu_tracker;

class imu_tracker
{
public:
	imu_tracker() {}
	void Advance(uint64_t timestamp,
		const deque<SensorFusion::IMUdata>& imupool);
	void ComplementryFilter(imu_tracker* temp_,
		const deque<SensorFusion::IMUdata>& imupool, const int& front_search
		, const int& back_search, uint64_t timestamp);
	uint64_t last_time_;
	Vector3d pose;
	Vector3d gravity_direction;
	Vector3d rotation_velocity;
	static float gravity_norm;
};
class SensorFusion
{
public:
	// ��Ϊ���ﲻ���ڶԵ�ͼ���в���
	// ���Բ���ʹ��Array3d ����ʹ��Vector3d
	// IMU gyro������roll pitch yew ����Ԫ��ת������
	// IMU acc ���ݿ�������� pitch ��yew
	// imu_tracker pose ������roll pitch yew ����ʾ
	// imu_tracker ʹ����Ӧ�����˲�����
	struct IMUdata
	{
		Vector3d gyro_data;
		Vector3d acc_data;
		uint64_t timestamp;
	};
	struct Odometrydata
	{
		Vector3d pose;
		uint64_t timestamp;
	};
	struct Posedata
	{
		Vector3d pose;
		uint64_t timestamp;
	};
	SensorFusion(){}
	void AddPose(const Vector3d& temp_);
	void AddOdometryData(const Odometrydata& temp_);
	void AddImuData(const IMUdata& temp_);
	void TrimIMUdata(uint64_t timestamp);
	void TrimeOdometrydata(uint64_t timestamp);
	Vector3d PoseEstimation(uint64_t now_);
	const deque<IMUdata>& GetImupool() const { return imupool_; }
private:
	deque<IMUdata> imupool_;
	deque<Odometrydata> odometrypool_;
	deque<Posedata> posepool_;
	Vector3d linear_velocity;
	imu_tracker * pose_tracker;
	imu_tracker * imu_tracker_;
	imu_tracker * odometry_tracker;
};




