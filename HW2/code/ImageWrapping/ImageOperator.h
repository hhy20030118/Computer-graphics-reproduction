#pragma once
#include "qvector.h"
class pairpos {//�洢���ڲ����㣨�����㣬�ߵȵȣ�������
public:
	QPoint bpos , epos;
};
class ImageOperator
{
	
	QVector<QVector<int>> Picture;  //��һ����ά��������ͼ��
	
public:
	ImageOperator(int h, int w);
	ImageOperator();	
	int height;
	int width;
	QVector<pairpos> ChangedPos; //������ļ���
};

