#pragma once
#include <iostream>
#include "Eigen/Eigen"
#include <algorithm>
using namespace std;
using namespace Eigen;
namespace tree
{
	// represent in the firstchild and rightsib method and in order to find the parent quickly,
	// Add a pointer for parent 
	struct Data
	{
		Data(Vector2d data_,bool flag_ ,
			Data* leftchild_ = nullptr, Data* rightchild_ = nullptr, Data* parent_ = nullptr) :data(data_),
			parent(parent_), leftchild(leftchild_),
			rightchild(rightchild_) ,XYflag(flag_),checkflag(0){}
		Vector2d data;
		Data* parent;
		Data* leftchild;
		Data* rightchild;
		bool XYflag;
		bool checkflag;
	};
	
	// Kd-tree main part
	class Kdtree
	{
	public:
		Kdtree(const int & num):k(num) { downhill = true; }
		Kdtree(){}
		~Kdtree() {}
		Kdtree(const vector<Vector2d>& temp,const int & num );
		void Initialization(const vector<Vector2d>& temp, const int num =3);
		typedef vector<Vector2d>::iterator Iter;
		vector<Vector2d> waiting_set;
		Vector2d CPoint;
		Data* root_ptr;
		bool XYflag;
		int k;
		vector<Vector2d> ascendingX;
		vector<Vector2d> ascendingY;
		vector<Data*> ptr_pool;
		vector<Data*> checklist_;
		bool downhill;
		bool ReorderFinished;
	    void Reorder(const vector<Vector2d>& temp);
		Data* CreateTree( uint32_t start, uint32_t end ,bool flag_);
		void InsertNode(Vector2d temp_);
		void Finish();
		void DestroyTree();
		double CalculateDistance(Vector2d p1, Vector2d p2);
	    static void Preordershow(const Data* ptr);
		void kNN( vector<Vector2d>& temp, int k,  Data* ptr);
		friend vector<Vector2d> findKNN( Kdtree& temp, const Vector2d& cpoint);
	};
	vector<Vector2d> findKNN(Kdtree& temp, const Vector2d& cpoint);
}

