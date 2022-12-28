#pragma once
#include <vector>
#include <QPainter>
#include <QMouseEvent>
class PointAndLine//关于操作点和操作线，和minidraw有点像
{
public:
    int type;
    QVector<QPoint>path;
    QPen pen;
    PointAndLine() {

    }
    ~PointAndLine() {

    }
    void draw(QPoint pos) { //自由绘制
        path.append(pos);
    }
    void draw(QPoint pos1 , QPoint pos2) {  //绘制边框
        path.replace(0 , pos1);
        QPoint p(pos1.x(), pos2.y());
        path.replace(1, p);//对应长方体的四个点
        QPoint q(pos2.x(), pos1.y());
        path.replace(3, q);
        path.replace(2, pos2);
        path.replace(4, path.at(0));
    }
};


