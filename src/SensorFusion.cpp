#include "SensorFusion.h"
float imu_tracker::gravity_norm = 9.81;

void imu_tracker::Advance(uint64_t timestamp,
	const deque<SensorFusion::IMUdata>& imupool)
{
	// 如果IMU 组件出问题或者多线程的时候出现问题
	// 或者AddPose pose的时间比较靠前的话
	if (imupool.empty())
		return;
	// find all imu data that satisfy the timestamp condition
	int back_search = imupool.size()-1;
	bool flag_change = false;
	int front_search = 0;
	for (;front_search!=back_search;)
	{
		if (imupool[front_search].timestamp < last_time_)
		{
			front_search++;
			flag_change = true;
		}
		if (imupool[back_search].timestamp > timestamp)
		{
			back_search--;
			flag_change = true;
		}
		if (flag_change == false)
			break;
		flag_change = false;
	}
	if (front_search == back_search)
		return;
	// 找到合适的IMU数据集之后对imu_tracker进行互补滤波
	ComplementryFilter( this, imupool, front_search, back_search, timestamp);
}

//互补滤波主体
void imu_tracker::ComplementryFilter(imu_tracker* temp_,
	const deque<SensorFusion::IMUdata> & imupool, const int & front_search
	,const int & back_search,uint64_t timestamp)
{
	assert(front_search < back_search);
	// 处理数据
	// gyro_data_ 进行角度修正
	// roll 无法估计且紧贴地面 roll一直为0;
	// 这里只对pitch yew 来进行
    // 先对pose 的首信息进行处理
	int front = front_search;
	float gyro_coeff = 0.8;
	float acc_coeff = 0.2;
	float gyro_coeff_high_acc = 0.995;
	float acc_coeff_high_acc = 0.005;
	float acc_roll = atan2(imupool[front_search].acc_data(1),-1*imupool[front_search].acc_data(2));
	float acc_pitch = atan2(-imupool[front_search].acc_data(0),-imupool[front_search].acc_data(2));
	uint64_t delta_t = last_time_ - imupool[front_search].timestamp;
	pose[0] = gyro_coeff * ( pose[0] + rotation_velocity(0) * delta_t * 1e-6 )
			  + acc_coeff * acc_roll;
	pose[1] = gyro_coeff * ( pose[1] + rotation_velocity(1) * delta_t * 1e-6)
		      + acc_coeff * acc_pitch;
	pose[2] = (pose[2] + rotation_velocity(1) * delta_t * 1e-6);
	front++;

	// 处理之后的数据
	while (front < back_search)
	{
		acc_roll = atan2(imupool[front].acc_data(1), -1 * imupool[front].acc_data(2));
		acc_pitch = atan2(-1 * imupool[front].acc_data(0), -1 * imupool[front].acc_data(2));
		delta_t = imupool[front].timestamp - imupool[front - 1].timestamp;
		if (imupool[front].acc_data.norm() < 1.15 * gravity_norm)
		{
			pose[0] = gyro_coeff * (pose[1] + imupool[front - 1].gyro_data(0) * delta_t * 1e-6)
				+ acc_coeff * acc_roll;
			pose[1] = gyro_coeff * (pose[1] + imupool[front - 1].gyro_data(1) * delta_t * 1e-6)
				+ acc_coeff * acc_pitch;
			pose[2] =  (pose[2] + imupool[front - 1].gyro_data(2) * delta_t * 1e-6);
		}
		else
		{
			pose[0] = gyro_coeff_high_acc * (pose[1] + imupool[front - 1].gyro_data(0) * delta_t * 1e-6)
				+ acc_coeff_high_acc * acc_roll;
			pose[1] = gyro_coeff_high_acc * (pose[1] + imupool[front - 1].gyro_data(1) * delta_t * 1e-6)
				+ acc_coeff_high_acc * acc_pitch;
			pose[2] = (pose[2] + imupool[front - 1].gyro_data(2) * delta_t * 1e-6);
		}
		front++;
	}
	
	// 处理末尾数据
	rotation_velocity = imupool[front].gyro_data;
	last_time_ = imupool[front].timestamp;
}

void SensorFusion::AddImuData(const IMUdata& data_)
{
	if (imupool_.empty())
	{
		imupool_.push_back(data_);
	}
	assert(imupool_.back().timestamp < data_.timestamp);
	imupool_.push_back(data_);
}

// 这里和Carto一样还是只是用两个Odo 信息 
void SensorFusion::AddOdometryData(const Odometrydata& data_)
{
	if (odometrypool_.size() < 2)
		odometrypool_.push_back(data_);
	return;
	
	// 如果还没有为 odometry建立imu_tracker 
	odometry_tracker = new imu_tracker(*imu_tracker_);
	
	// 超过两个就开始正常的预测过程
	odometrypool_.pop_front();
	odometrypool_.push_back(data_);
	Vector3d delta_trans = odometrypool_.front().pose - odometrypool_.back().pose;
	
	// 处理imu 数据
	odometry_tracker->Advance(odometrypool_.back().timestamp, imupool_);
	
	// 速度大小
	Vector3d  linear_velocity_from_odometry = delta_trans / (1e-6 * 
		(odometrypool_.back().timestamp - odometrypool_.front().timestamp));
	Vector3d  direction = odometry_tracker->pose + odometry_tracker->rotation_velocity * (odometrypool_.back().timestamp
		- odometry_tracker->last_time_) * 1e-6;
}
