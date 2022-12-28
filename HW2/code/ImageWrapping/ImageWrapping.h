#pragma once

#include <QtWidgets/QMainWindow>
#include "ui_ImageWrapping.h"
#include "PointAndLine.h"
#include <QPainter>
#include <QMouseEvent>
#include "ImageOperator.h"
class ImageWrapping : public QMainWindow
{
    Q_OBJECT
protected://事件函数
    void paintEvent(QPaintEvent* ev);
    void mouseMoveEvent(QMouseEvent* ev);
    void mousePressEvent(QMouseEvent* ev);
public:
    ImageWrapping(QWidget *parent = Q_NULLPTR);
    int type;
public slots://触发函数
    void OpenFile();
    void ChangeTypeToPoint();
    void ChangeTypeToLine();
    void undo();
    void IDWoperate();
    void RBFoperate();
    void renew();
    void save();
private:
    Ui::ImageWrappingClass ui;
    QVector<PointAndLine*> graphics;//操作点线
    int pictureopened = 0;
    QPen pen_pub;//公用笔
    QImage myImage;//当前图像
    QImage myImage_origin;//原始图像
    QPixmap map;
    ImageOperator o;//操作器
    
};
