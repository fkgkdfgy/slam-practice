#pragma once
#include "PGM.h"
// �����Ӧ�����ڷ�֧�����������д
// ���������Ӧ�ý��������Լ���ͽ����Ż�����Ӧ���������ٲ���Ķ���
class PoseGraph
{
public:
	void AddScan();
	void ComputeConstraint();
	void Solve();


};

