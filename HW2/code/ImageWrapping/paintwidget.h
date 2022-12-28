#pragma once
#include <qwidget.h>
#include <QtWidgets/QMainWindow>
#include "PointAndLine.h"
#include <QLabel>
#include <QPoint>
#include <QColor>
#include <QPaintEvent>
#include <QImage>
#include <QPixmap>
#include <QMouseEvent>
class paintwidget : public QLabel
{
public:
    //bool eventFilter(QObject* obj, QEvent* event);
protected:
    //virtual void paintEvent(QPaintEvent* event) override;
    //void paintEvent(QPaintEvent* ev);
    //void mouseMoveEvent(QMouseEvent* ev);
    //void mousePressEvent(QMouseEvent* ev);
public:
    QVector<PointAndLine*>* graphics;
    QPen pen_pub;
};

