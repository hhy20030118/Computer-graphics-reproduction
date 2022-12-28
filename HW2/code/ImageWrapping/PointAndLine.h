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
};


