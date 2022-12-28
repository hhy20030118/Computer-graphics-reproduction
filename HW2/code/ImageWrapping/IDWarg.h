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
    TYPE inline d(QPointF a, QPointF b) {//二者之间的距离
        return sqrt((((TYPE)a.x() - (TYPE)b.x()) * ((TYPE)a.x() - (TYPE)b.x()) + ((TYPE)a.y() - (TYPE)b.y()) * ((TYPE)a.y() - (TYPE)b.y())));
    }
    TYPE inline getW(int i , QPointF p) {  //算权重
        TYPE e = 0, sum = 0 , u;
        for (register int j = 0; j < ChangedPos.size(); ++j) {
            u = 1.0 / pow((d(p, QPointF(ChangedPos.at(j).bpos))),1);//普通径向函数
            //u = pow((R- d(p, QPointF(ChangedPos.at(j).bpos)))/(R* d(p, QPointF(ChangedPos.at(j).bpos))), 1); //locally bounded weight 函数
            sum += u;
            if (j == i)e = u;
        }
        return e / sum;
    }
    QPointF inline getF(int i, QPointF p) {  //算对应点
        return ChangedPos.at(i).epos + p - ChangedPos.at(i).bpos;
    }
public:
    QPointF inline getpos(int i , int j){//传入坐标，传出更改后的坐标
        y = QPointF(i, j);
        x = QPointF(0, 0);
        for (register int a = 0; a < ChangedPos.size(); ++a) {
            TYPE r = getW(a, y);
            QPointF t= getF(a, y);
            x += r * t;  //对应点相乘
        }
        return x;
    }
};

