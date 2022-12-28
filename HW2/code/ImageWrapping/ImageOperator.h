#pragma once
#include "qvector.h"
class pairpos {//存储关于操作点（不动点，线等等）的数据
public:
	QPoint bpos , epos;
};
class ImageOperator
{
	
	QVector<QVector<int>> Picture;  //用一个二维向量储存图像
	
public:
	ImageOperator(int h, int w);
	ImageOperator();	
	int height;
	int width;
	QVector<pairpos> ChangedPos; //操作点的集合
};

