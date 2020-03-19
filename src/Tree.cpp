#include  "Tree.h"

namespace tree
{
	void Kdtree::InsertNode(Vector2d temp_)
	{
		// 找到合适的位置
		// checkflag = 0 对比X checkflag = 1 对比Y
		Data* temp_ptr = root_ptr;
		Data* parent_ptr = root_ptr->parent;
		while (temp_ptr != nullptr)
		{
			if (temp_ptr->XYflag == 0)
			{
				if (temp_(0) < temp_ptr->data(0))
				{
					parent_ptr = temp_ptr;
					temp_ptr = temp_ptr->leftchild;
				}
				else
				{
					parent_ptr = temp_ptr;
					temp_ptr = temp_ptr->rightchild;
				}
			}
			else
			{
				if (temp_(1) < temp_ptr->data(1))
				{
					parent_ptr = temp_ptr;
					temp_ptr = temp_ptr->leftchild;
				}
				else
				{
					parent_ptr = temp_ptr;
					temp_ptr = temp_ptr->rightchild;
				}
			}
		}
		temp_ptr = new Data(temp_, !(parent_ptr->XYflag), nullptr, nullptr, parent_ptr);
		if (parent_ptr->XYflag == 0)
		{
			if (temp_(0) < parent_ptr->data(0))
			{
				parent_ptr->leftchild = temp_ptr;

			}
			else
			{
				parent_ptr->rightchild = temp_ptr;
			}
		}
		else
		{
			if (temp_(1) < parent_ptr->data(1))
			{
				parent_ptr->leftchild = temp_ptr;

			}
			else
			{
				parent_ptr->rightchild = temp_ptr;
			}
		}
		ptr_pool.push_back(temp_ptr);
	}


	void Kdtree::DestroyTree()
	{
		//TEST PART
		cout << "show the quantity of vector pool " << ptr_pool.size() << endl;
		int size = ptr_pool.size();
		while (size != 0)
		{
			delete ptr_pool.back();
			ptr_pool.pop_back();
			size--;
		}
		cout << "now the quantity of vector pool is " << ptr_pool.size() << endl;
	}

	void Kdtree::Initialization(const vector<Vector2d>& temp, const int num)
	{
		downhill = true;
		assert(temp.size() >= 2);
		Reorder(temp);
		ptr_pool.reserve(temp.size());
		this->root_ptr = CreateTree(0, temp.size() - 1, 0);
	}

	Kdtree::Kdtree(const vector<Vector2d>& temp, const int& num = 3) :k(num), downhill(true)
	{
		assert(temp.size() >= 2);
		Reorder(temp);
		ptr_pool.reserve(temp.size());
	}
	// Create the kdtree
	Data* Kdtree::CreateTree(uint32_t start, uint32_t end
		, bool flag_)
	{
		bool XYflag = flag_;
		uint32_t temp = (end - start + 1) / 2;
		// Baseline condition
		// end -start = 0 就是说只有一个数进来了
		// end -start = 1 就是说有两个数进来了
		if (end - start == 0)
		{
			if (XYflag == 1)
			{
				auto temp = new Data(ascendingY[start], flag_);
				ptr_pool.push_back(temp);
				return temp;
			}
			else
			{
				auto temp = new Data(ascendingX[start], flag_);
				ptr_pool.push_back(temp);
				return temp;
			}
		}
		else if (end - start == 1)
		{
			if (XYflag == 1)
			{
				auto temp = new Data(ascendingY[end], flag_, CreateTree(start, start, !flag_));
				temp->leftchild->parent = temp;
				ptr_pool.push_back(temp);
				return temp;
			}
			else
			{
				auto temp = new Data(ascendingX[end], flag_, CreateTree(start, start, !flag_));
				temp->leftchild->parent = temp;
				ptr_pool.push_back(temp);
				return temp;
			}
		}

		// Based on the ascending order
		// Reordering the Vector2d
		// if the XYflag is equal to 0 ,sorting the x 
		// if the XYflag is equal to 1，sorting the y
		if (XYflag == 0)
		{
			Data* temp_1 = new Data(ascendingX[start + temp], flag_, CreateTree(start, (start + temp - uint32_t(1)), !flag_), CreateTree(start + temp + uint32_t(1), end, !flag_));
			if (temp_1->leftchild != nullptr)
				temp_1->leftchild->parent = temp_1;
			if (temp_1->rightchild != nullptr)
				temp_1->rightchild->parent = temp_1;
			ptr_pool.push_back(temp_1);
			return temp_1;
		}
		else
		{
			Data* temp_1 = new Data(ascendingY[start + temp], flag_, CreateTree(start, (start + temp - uint32_t(1)), !flag_), CreateTree(start + temp + uint32_t(1), end, !flag_));
			if (temp_1->leftchild != nullptr)
				temp_1->leftchild->parent = temp_1;
			if (temp_1->rightchild != nullptr)
				temp_1->rightchild->parent = temp_1;
			ptr_pool.push_back(temp_1);
			return temp_1;
		}
	}

	// 重新排序加速树的构建
	void Kdtree::Reorder(const vector<Vector2d>& temp)
	{

		ascendingX.assign(temp.begin(), temp.end());
		ascendingY.assign(temp.begin(), temp.end());
		assert(temp.size() >= 2);
		sort(ascendingX.begin(), ascendingX.end(), [](Vector2d a, Vector2d b) {return a.x() < b.x(); });
		sort(ascendingY.begin(), ascendingY.end(), [](Vector2d a, Vector2d b) {return a.y() < b.y(); });
		ReorderFinished = 1;
	}

	// 前序遍历
	void Kdtree::Preordershow(const Data* ptr)
	{
		if (ptr == nullptr)
			return;
		cout << " This is " << endl;
		cout << ptr->data(0) << "   " << ptr->data(1) << endl;
		if (ptr->leftchild != nullptr)
		{
			cout << " the Left child is " << endl;
			cout << ptr->leftchild->data(0) << "   " << ptr->leftchild->data(1) << endl;
		}
		else
		{
			cout << " Dont have left child" << endl;
		}
		if (ptr->rightchild != nullptr)
		{
			cout << "the Right child is" << endl;
			cout << ptr->rightchild->data(0) << "   " << ptr->rightchild->data(1) << endl;
		}
		else
		{
			cout << "Dont have right child" << endl;
		}
		Preordershow(ptr->leftchild);
		Preordershow(ptr->rightchild);
	}
	vector<Vector2d> findKNN(Kdtree& temp, const Vector2d& cpoint)
	{
		vector<Vector2d> v1;
		temp.CPoint = cpoint;
		v1.reserve(temp.k);
		temp.downhill = true;
		temp.kNN(v1, temp.k, temp.root_ptr);
		temp.Finish();
		return v1;
	}
	// 寻找KNN的主体函数
	void Kdtree::kNN(vector<Vector2d>& temp, int size, Data* ptr)
	{

		// 问题是 ptr 如果是空指针的话 ptr->checkflg  根本就不会存在所以会报错
		// 所以选择把 两者顺序交换成如下的顺序)
		if (ptr == nullptr || ptr->checkflag == true)
		{
			downhill = false;
			return;
		}

		// 根节点的处理基本结束
		if (ptr->leftchild == nullptr && ptr->rightchild == nullptr)
		{
			ptr->checkflag = true;
			checklist_.push_back(ptr);
			downhill = false;
			if (temp.size() < size)
			{
				temp.push_back(ptr->data);
				if (temp.size() > 1)
					sort(temp.begin(), temp.end(), [this](Vector2d a, Vector2d b) {return Vector2d(this->CPoint - a).norm() < Vector2d(this->CPoint - b).norm(); });
				return;
			}
			else
			{
				if (Vector2d(this->CPoint - ptr->data).norm() < Vector2d(this->CPoint - temp.back()).norm())
				{
					temp.pop_back();
					temp.push_back(ptr->data);
					sort(temp.begin(), temp.end(), [this](Vector2d a, Vector2d b) {return Vector2d(this->CPoint - a).norm() < Vector2d(this->CPoint - b).norm(); });
					return;
				}
				return;
			}
		}

		// 非根节点
		if (downhill == true)
		{
			// compare x
			if (ptr->XYflag == 0)
			{
				if (CPoint.x() < ptr->data.x())
					kNN(temp, size, ptr->leftchild);
				else
					kNN(temp, size, ptr->rightchild);
			}
			else // compare y
			{
				if (CPoint.y() < ptr->data.y())
					kNN(temp, size, ptr->leftchild);
				else
					kNN(temp, size, ptr->rightchild);
			}
			assert(downhill == false);
			ptr->checkflag = true;
			checklist_.push_back(ptr);
			// add points
			if (temp.size() < size)
			{
				temp.push_back(ptr->data);
				sort(temp.begin(), temp.end(), [this](Vector2d a, Vector2d b) {return Vector2d(this->CPoint - a).norm() < Vector2d(this->CPoint - b).norm(); });
			}
			else
			{
				if ((CPoint - ptr->data).norm() < (CPoint - temp.back()).norm())
				{
					temp.pop_back();
					temp.push_back(ptr->data);
					sort(temp.begin(), temp.end(), [this](Vector2d a, Vector2d b) {return Vector2d(this->CPoint - a).norm() < Vector2d(this->CPoint - b).norm(); });
				}
			}

			if (ptr->XYflag == 0)
			{
				// < 则说明应该进行查找	
				if (std::abs(CPoint.x() - ptr->data.x()) < (CPoint - ptr->data).norm())
				{
					// 因为已经是从一端上来的 再次下到上来的一端会直接返回所以这里也就不需要再进行条件判断
					downhill = true;
					kNN(temp, size, ptr->leftchild);
					downhill = true;
					kNN(temp, size, ptr->rightchild);
					return;
				}
				else
				{
					downhill = false;
					return;
				}
			}
			else
			{
				// < 则说明应该进行查找	
				if (std::abs(CPoint.y() - ptr->data.y()) < (CPoint - ptr->data).norm())
				{
					// 因为已经是从一端上来的 再次下到上来的一端会直接返回所以这里也就不需要再进行条件判断
					downhill = true;
					kNN(temp, size, ptr->leftchild);
					downhill = true;
					kNN(temp, size, ptr->rightchild);
					return;
				}
				else
				{
					downhill = false;
					return;
				}
			}
		}
	}

	void Kdtree::Finish()
	{
		for (auto& element : checklist_)
		{
			checklist_.back()->checkflag = false;
			checklist_.pop_back();
		}
	}
}