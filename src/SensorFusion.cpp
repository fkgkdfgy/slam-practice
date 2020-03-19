#include "SensorFusion.h"
float imu_tracker::gravity_norm = 9.81;

void imu_tracker::Advance(uint64_t timestamp,
	const deque<SensorFusion::IMUdata>& imupool)
{
	// ���IMU �����������߶��̵߳�ʱ���������
	// ����AddPose pose��ʱ��ȽϿ�ǰ�Ļ�
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
	// �ҵ����ʵ�IMU���ݼ�֮���imu_tracker���л����˲�
	ComplementryFilter( this, imupool, front_search, back_search, timestamp);
}

//�����˲�����
void imu_tracker::ComplementryFilter(imu_tracker* temp_,
	const deque<SensorFusion::IMUdata> & imupool, const int & front_search
	,const int & back_search,uint64_t timestamp)
{
	assert(front_search < back_search);
	// ��������
	// gyro_data_ ���нǶ�����
	// roll �޷������ҽ������� rollһֱΪ0;
	// ����ֻ��pitch yew ������
    // �ȶ�pose ������Ϣ���д���
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

	// ����֮�������
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
	
	// ����ĩβ����
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

// �����Cartoһ������ֻ��������Odo ��Ϣ 
void SensorFusion::AddOdometryData(const Odometrydata& data_)
{
	if (odometrypool_.size() < 2)
		odometrypool_.push_back(data_);
	return;
	
	// �����û��Ϊ odometry����imu_tracker 
	odometry_tracker = new imu_tracker(*imu_tracker_);
	
	// ���������Ϳ�ʼ������Ԥ�����
	odometrypool_.pop_front();
	odometrypool_.push_back(data_);
	Vector3d delta_trans = odometrypool_.front().pose - odometrypool_.back().pose;
	
	// ����imu ����
	odometry_tracker->Advance(odometrypool_.back().timestamp, imupool_);
	
	// �ٶȴ�С
	Vector3d  linear_velocity_from_odometry = delta_trans / (1e-6 * 
		(odometrypool_.back().timestamp - odometrypool_.front().timestamp));
	Vector3d  direction = odometry_tracker->pose + odometry_tracker->rotation_velocity * (odometrypool_.back().timestamp
		- odometry_tracker->last_time_) * 1e-6;
}
