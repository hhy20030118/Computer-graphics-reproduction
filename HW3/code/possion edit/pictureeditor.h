#pragma once
#include <Dense>
#include <Core>
#include <QImage>
#include <vector>
#include <QPainter>
#include <SparseCholesky>
#include <SparseLU>
using namespace Eigen;
using namespace std;
class pictureeditor
{
private:
	long n;
	SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	QVector<QVector<short>>regionBool;
	QImage image1 , image2 , image3;
	SparseMatrix<double> a1;
	VectorXd x1= Vector<double, Dynamic>(), x2= Vector<double, Dynamic>(),x3=Vector<double, Dynamic>(), b1 = Vector<double, Dynamic>(),b2 = Vector<double, Dynamic>(),b3 = Vector<double, Dynamic>();
public:
	inline int correct(double x) {  //对rgb进行纠正
		x = (int)x;
		if (x > 255) {
			x = 255;
		}
		else if (x < 0) {
			x = 0;
		}
		return x;
	}
	void setMatrix(SparseMatrix<double> a) {
		a1 = a;
	}
	void setregionBool(QVector<QVector<short>> b) {
		regionBool = b;
	}
	void setimage1(QImage im) {  //设置初始图象（梯度场来源）
		image1 = im;
		image3 = im;
	}
	void setimage2(QImage im) {  //设置目标图像（边界值来源）
		image2 = im;

	}
	QImage getimage() {
		return image3;
	}
	void setnum(long h) {
		n = h;
	}
	void analysis() {
		long k = 0;
		for (register int i = 0; i < regionBool.size(); ++i) {
			for (register int j = 0; j < regionBool.at(0).size(); ++j) {
				if (regionBool.at(i).at(j) == 2) {
					b1(k) = double(qRed(image2.pixel(j , i)));  //边界值相同
					b2(k) = double(qGreen(image2.pixel(j, i)));
					b3(k) = double(qBlue(image2.pixel(j, i)));
					k++;
				}
				else if (regionBool.at(i).at(j) == 1) {
					b1(k) = 0;
					b2(k) = 0;
					b3(k) = 0;
					
					//red
					if (i != 0 && (abs(qRed(image1.pixel(j, i - 1)) - qRed(image1.pixel(j, i))) < abs(qRed(image2.pixel(j, i - 1)) - qRed(image2.pixel(j, i)))))b1(k) += qRed(image2.pixel(j, i - 1)) - qRed(image2.pixel(j, i));
					else if (i != 0)b1(k) += qRed(image1.pixel(j, i - 1)) - qRed(image1.pixel(j, i));
					
					if (i != regionBool.size() - 1 && (abs(qRed(image1.pixel(j , i+1)) - qRed(image1.pixel(j, i))) < abs(qRed(image2.pixel(j , i+1)) - qRed(image2.pixel(j, i)))))b1(k) += qRed(image2.pixel(j , i+1)) - qRed(image2.pixel(j, i));
					else if (i != regionBool.size() - 1)b1(k) += qRed(image1.pixel(j, i + 1)) - qRed(image1.pixel(j, i));

					if (j != 0 && (abs(qRed(image1.pixel(j - 1, i)) - qRed(image1.pixel(j, i))) < abs(qRed(image2.pixel(j - 1, i)) - qRed(image2.pixel(j, i)))))b1(k) += qRed(image2.pixel(j - 1, i)) - qRed(image2.pixel(j, i));
					else if (j != 0)b1(k) += qRed(image1.pixel(j-1, i)) - qRed(image1.pixel(j, i));
					if (j != regionBool.at(0).size() - 1 && (abs(qRed(image1.pixel(j + 1, i)) - qRed(image1.pixel(j, i))) < abs(qRed(image2.pixel(j + 1, i)) - qRed(image2.pixel(j, i)))))b1(k) += qRed(image2.pixel(j + 1, i)) - qRed(image2.pixel(j, i));
					else if (j != regionBool.at(0).size() - 1)b1(k) += qRed(image1.pixel(j+1, i)) - qRed(image1.pixel(j, i));
					//green
					if (i != 0 && (abs(qGreen(image1.pixel(j, i - 1)) - qGreen(image1.pixel(j, i))) < abs(qGreen(image2.pixel(j, i - 1)) - qGreen(image2.pixel(j, i)))))b2(k) += qGreen(image2.pixel(j, i - 1)) - qGreen(image2.pixel(j, i));
					else if (i != 0)b2(k) += qGreen(image1.pixel(j, i - 1)) - qGreen(image1.pixel(j, i));
					if (i != regionBool.size() - 1 && (abs(qGreen(image1.pixel(j, i + 1)) - qGreen(image1.pixel(j, i))) < abs(qGreen(image2.pixel(j, i + 1)) - qGreen(image2.pixel(j, i)))))b2(k) += qGreen(image2.pixel(j , i+1)) - qGreen(image2.pixel(j, i));
					else if (i != regionBool.size() - 1)b2(k) += qGreen(image1.pixel(j, i + 1)) - qGreen(image1.pixel(j, i));

					if (j != 0 &&( abs(qGreen(image1.pixel(j - 1, i)) - qGreen(image1.pixel(j, i))) < abs(qGreen(image2.pixel(j - 1, i)) - qGreen(image2.pixel(j, i)))))b2(k) += qGreen(image2.pixel(j - 1, i)) - qGreen(image2.pixel(j, i));
					else if (j != 0)b2(k) += qGreen(image1.pixel(j - 1, i)) - qGreen(image1.pixel(j, i));
					if (j != regionBool.at(0).size() - 1 &&(abs(qGreen(image1.pixel(j + 1, i)) - qGreen(image1.pixel(j, i))) < abs(qGreen(image2.pixel(j + 1, i)) - qGreen(image2.pixel(j, i)))))b2(k) += qGreen(image2.pixel(j + 1, i)) - qGreen(image2.pixel(j, i));
					else if (j != regionBool.at(0).size() - 1)b2(k) += qGreen(image1.pixel(j + 1, i)) - qGreen(image1.pixel(j, i));
					//blue
					if (i != 0 &&( abs(qBlue(image1.pixel(j, i - 1)) - qBlue(image1.pixel(j, i))) < abs(qBlue(image2.pixel(j, i - 1)) - qBlue(image2.pixel(j, i)))))b3(k) += qBlue(image2.pixel(j, i - 1)) - qBlue(image2.pixel(j, i));
					else if (i != 0)b3(k) += qBlue(image1.pixel(j, i - 1)) - qBlue(image1.pixel(j, i));
					if (i != regionBool.size() - 1 && (abs(qBlue(image1.pixel(j, i + 1)) - qBlue(image1.pixel(j, i))) < abs(qBlue(image2.pixel(j, i + 1)) - qBlue(image2.pixel(j, i)))))b3(k) += qBlue(image2.pixel(j , i+1)) - qBlue(image2.pixel(j, i));
					else if (i != regionBool.size() - 1)b3(k) += qBlue(image1.pixel(j, i + 1)) - qBlue(image1.pixel(j, i));

					if (j != 0 && (abs(qBlue(image1.pixel(j - 1, i)) - qBlue(image1.pixel(j, i))) < abs(qBlue(image2.pixel(j - 1, i)) - qBlue(image2.pixel(j, i)))))b3(k) += qBlue(image2.pixel(j - 1, i)) - qBlue(image2.pixel(j, i));
					else if (j != 0)b3(k) += qBlue(image1.pixel(j - 1, i)) - qBlue(image1.pixel(j, i));
					if (j != regionBool.at(0).size() - 1 && (abs(qBlue(image1.pixel(j + 1, i)) - qBlue(image1.pixel(j, i))) < abs(qBlue(image2.pixel(j + 1, i)) - qBlue(image2.pixel(j, i)))))b3(k) += qBlue(image2.pixel(j + 1, i)) - qBlue(image2.pixel(j, i));
					else if (j != regionBool.at(0).size() - 1)b3(k) += qBlue(image1.pixel(j + 1, i)) - qBlue(image1.pixel(j, i));
					
					if (i != 0 && regionBool.at(i - 1).at(j) == 2) {
						b1(k) -= qRed(image2.pixel(j, i-1));
						b2(k) -= qGreen(image2.pixel(j, i-1));
						b3(k) -= qBlue(image2.pixel(j, i-1));
					}
					if (i != regionBool.size() - 1 && regionBool.at(i + 1).at(j) == 2) {
						b1(k) -= qRed(image2.pixel(j, i + 1));
						b2(k) -= qGreen(image2.pixel(j, i + 1));
						b3(k) -= qBlue(image2.pixel(j, i + 1));
					}
					if (j != 0 && regionBool.at(i).at(j - 1) == 2) {
						b1(k) -= qRed(image2.pixel(j-1, i));
						b2(k) -= qGreen(image2.pixel(j-1, i ));
						b3(k) -= qBlue(image2.pixel(j-1, i ));
					}
					if (j != regionBool.at(0).size() - 1 && regionBool.at(i).at(j + 1) == 2) {
						b1(k) -= qRed(image2.pixel(j + 1, i));
						b2(k) -= qGreen(image2.pixel(j + 1, i));
						b3(k) -= qBlue(image2.pixel(j + 1, i));
					}
					k++;
				}
			}
		}
		x1 = solver.solve(b1);
		x2 = solver.solve(b2); 
		x3 = solver.solve(b3);
		
		k = 0;
		image3 = image2;
		for (register int i = 0; i < regionBool.size(); ++i) {
			for (register int j = 0; j < regionBool.at(0).size(); ++j) {
				if (regionBool.at(i).at(j) != 0) {//更新图像
					image3.setPixel(j, i, qRgb(correct(x1(k)), correct(x2(k)), correct(x3(k))));
					k++;
				}
			}
		}
	}
	void initial() {  //对solver初始化
		b1.resize(n);
		b2.resize(n);
		b3.resize(n);
		x1.resize(n);
		x2.resize(n);
		x3.resize(n);
		a1.makeCompressed();
		solver.compute(a1);
		if (solver.info() != Eigen::Success)
		{
			throw std::exception("Compute Matrix is error");
			return;
		}
	}
};

