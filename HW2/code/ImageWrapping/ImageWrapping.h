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
protected://�¼�����
    void paintEvent(QPaintEvent* ev);
    void mouseMoveEvent(QMouseEvent* ev);
    void mousePressEvent(QMouseEvent* ev);
public:
    ImageWrapping(QWidget *parent = Q_NULLPTR);
    int type;
public slots://��������
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
    QVector<PointAndLine*> graphics;//��������
    int pictureopened = 0;
    QPen pen_pub;//���ñ�
    QImage myImage;//��ǰͼ��
    QImage myImage_origin;//ԭʼͼ��
    QPixmap map;
    ImageOperator o;//������
    
};
