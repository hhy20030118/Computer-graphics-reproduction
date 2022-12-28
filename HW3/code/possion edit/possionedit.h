#pragma once
#include <Dense>
#include <Core>
#include <QtWidgets/QMainWindow>
#include "ui_possionedit.h"
#include "PointAndLine.h"
#include <SparseCholesky>
#include <Sparse>
#include <SparseLU>
#include "modeleditor.h"
#include "pictureeditor.h"
#include <vector>
using namespace Eigen;
using namespace std;
class possionedit : public QMainWindow
{
    Q_OBJECT
protected://事件函数
    void paintEvent(QPaintEvent* ev);
    void mouseMoveEvent(QMouseEvent* ev);
    void mousePressEvent(QMouseEvent* ev);
    void mouseReleaseEvent(QMouseEvent* ev);
public:
    possionedit(QWidget *parent = Q_NULLPTR);
    int type;  //标记是否在绘图，绘制的是什么
    pictureeditor p;  //对矩阵进行操作(原始区域）
    modeleditor e;   //对向量进行操作（目标区域）
public slots://触发函数
    void OnDraw();
    void Openpic_1();
    void Openpic_2();
    void undo();
    void possion();
private:

    typedef struct edge {  //扫描线算法相关
        double current;
        double delta_x;
        int ymax;
    }edge;
    vector<vector<edge>> NET;
    vector<edge> AET; //活动边表
    edge eg;
    short mapP[405][405];  //哪些是顶点，哪些会被计算

    QVector<QVector<short>> regionBool;  //区域的边界值为2，内部值为1，其余为0
    QVector<QVector<long>> region;  //沿着从上到下，从左到右对点进行标号，方便快速访问
    SimplicialLDLT<SparseMatrix<double>> solver;  //矩阵求解

    Ui::possioneditClass ui;

    void possion_pre();
    void scan();
    QPoint origin_pos;  //一开始拖动时鼠标在区域内点击的位置 ， 保持相对距离不变
    QPoint left_up;  //区域边框左上角
    QPoint right_down;  //区域边框右下角

    bool pictureopened1 = 0;  //图片是否被打开
    bool pictureopened2 = 0;

    bool OnMoving = 0;  //是否在拖拽
    long num = 0;
    short mapF[405][405];  //存目标区域的特征函数（区域内存1，区域外存0）
    QVector<PointAndLine*> graphics;//操作点线
    QPen pen_pub;//公用笔
    QImage myImage1;//当前图像
    QImage myImage1_origin;//原始图像
    QImage image_part2;  //在目标图像中截取的一部分
    QImage image_part1;  //在原始图像中截取的一部分
    QImage myImage2;//当前图像
    QImage myImage2_origin;//原始图像
    QPixmap map;
};
