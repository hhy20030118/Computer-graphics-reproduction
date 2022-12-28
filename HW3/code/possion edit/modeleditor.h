#pragma once
#include <Dense>
#include <Core>
#include <QImage>
#include <vector>
#include <SparseCholesky>
#include <Sparse>
#include <SparseLU>
using namespace Eigen;
using namespace std;
class modeleditor
{
private:

	QImage image;
	QVector<QVector<short>>regionBool;
	VectorXi sizes;
	QVector<QVector<long>>region;
	SparseMatrix<double> a1;
	long n;
public:
	SimplicialLDLT<SparseMatrix<double>> solver;
	void setregionBool(QVector<QVector<short>> b) {
		regionBool = b;
	}
	void setregion(QVector<QVector<long>> b) {
		region = b;
	}
	void setimage(QImage b) {
		image = b;
	}
	void setPosnum(long num) {
		n = num;
	}
	SparseMatrix<double> getmatrix() {
		return a1;
	}
	void initial() {
		long k = 0;
		a1.resize(n, n);
		sizes.resize(n);
		for (register long i = 0; i < n; ++i) {
			sizes(i) = 6;
		}
		a1.reserve(sizes);
		for (register int i = 0; i < regionBool.size(); ++i) {
			for (register int j = 0; j < regionBool.at(0).size(); ++j) {
				if (regionBool.at(i).at(j) == 2) {
					a1.insert(k , k ) = 1.0;   //普通
					k++;
				}
				else if (regionBool.at(i).at(j) == 1) {
					a1.insert(k, k) = -4.0;
					//需要判断邻居是不是在区域内部
					if (i != 0 && regionBool.at(i - 1).at(j) != 2) {
						a1.insert(k, region.at(i - 1).at(j)) = 1.0;
					}
					if (i != regionBool.size() - 1 && regionBool.at(i + 1).at(j) != 2) {
						a1.insert(k, region.at(i + 1).at(j)) = 1.0;
					}
					if (j != 0 && regionBool.at(i).at(j - 1) != 2) {
						a1.insert(k, k - 1) = 1.0;
					}
					if (j != regionBool.at(0).size() - 1 && regionBool.at(i).at(j + 1) != 2) {
						a1.insert(k, k + 1) = 1.0;
					}
					k++;
				}
			}
		}
		a1.makeCompressed();
	}
	
};

