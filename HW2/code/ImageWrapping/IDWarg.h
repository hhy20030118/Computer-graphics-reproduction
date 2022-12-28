#pragma once
#include "ImageOperator.h"
#include<qpoint.h>
#define TYPE float
class IDWarg :
    public ImageOperator
{
private:
    QPointF x , y;
    TYPE R = 5000.0;
    TYPE inline d(QPointF a, QPointF b) {//����֮��ľ���
        return sqrt((((TYPE)a.x() - (TYPE)b.x()) * ((TYPE)a.x() - (TYPE)b.x()) + ((TYPE)a.y() - (TYPE)b.y()) * ((TYPE)a.y() - (TYPE)b.y())));
    }
    TYPE inline getW(int i , QPointF p) {  //��Ȩ��
        TYPE e = 0, sum = 0 , u;
        for (register int j = 0; j < ChangedPos.size(); ++j) {
            u = 1.0 / pow((d(p, QPointF(ChangedPos.at(j).bpos))),1);//��ͨ������
            //u = pow((R- d(p, QPointF(ChangedPos.at(j).bpos)))/(R* d(p, QPointF(ChangedPos.at(j).bpos))), 1); //locally bounded weight ����
            sum += u;
            if (j == i)e = u;
        }
        return e / sum;
    }
    QPointF inline getF(int i, QPointF p) {  //���Ӧ��
        return ChangedPos.at(i).epos + p - ChangedPos.at(i).bpos;
    }
public:
    QPointF inline getpos(int i , int j){//�������꣬�������ĺ������
        y = QPointF(i, j);
        x = QPointF(0, 0);
        for (register int a = 0; a < ChangedPos.size(); ++a) {
            TYPE r = getW(a, y);
            QPointF t= getF(a, y);
            x += r * t;  //��Ӧ�����
        }
        return x;
    }
};

