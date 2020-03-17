#pragma once
#include "PGM.h"
// 这个类应该是在分支定界测试完再写
// 这个类做的应该仅仅是添加约束和进行优化，不应该在里面再掺别的东西
class PoseGraph
{
public:
	void AddScan();
	void ComputeConstraint();
	void Solve();


};

