#pragma once
#include <vector>
#include <QPainter>
#include <QMouseEvent>
class PointAndLine//���ڲ�����Ͳ����ߣ���minidraw�е���
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


