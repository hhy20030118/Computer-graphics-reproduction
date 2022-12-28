#include "possionedit.h"
#include "QFileDialog"
#include <QMessageBox>
#include <QImage>
#include <QPainter>
#include <QMouseEvent>
#include <QDesktopWidget>
#include <QApplication>
#include"PointAndLine.h"
#include"modeleditor.h"
#include"pictureeditor.h"
using namespace Eigen;
#define TYPE int
int cmp(const void* p1, const void* p2) {//比较函数
    if (*(TYPE*)p1 > *(TYPE*)p2)return 1;
    else if (*(TYPE*)p1 == *(TYPE*)p2)return 0;
    else return -1;
}
possionedit::possionedit(QWidget *parent)
    : QMainWindow(parent)
{
    ui.setupUi(this);
    ui.label_1->setStyleSheet("QLabel{border:2px solid rgb(255, 0, 0);}");
    ui.label_2->setStyleSheet("QLabel{border:2px solid rgb(255, 0, 0);}");
    connect(ui.actiondraw, SIGNAL(triggered()), this, SLOT(OnDraw()));
    connect(ui.pushButton_1, SIGNAL(clicked()), this, SLOT(Openpic_1()));
    connect(ui.pushButton_2, SIGNAL(clicked()), this, SLOT(Openpic_2()));
    connect(ui.actionundo, SIGNAL(triggered()), this, SLOT(undo()));
    memset(mapF, 0, sizeof(mapF));  //初始化
    memset(mapP, 0, sizeof(mapP));
    for (register int i = 0; i < 405; ++i) {
        vector<edge> n;
        NET.push_back(n);
    }
}
void possionedit::Openpic_1()//更新完毕
{
    QString fileName;
    QFileDialog* fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setFileMode(QFileDialog::AnyFile);
    fileDialog->setViewMode(QFileDialog::Detail);
    fileDialog->setGeometry(600, 600, 1000, 500);
    fileDialog->setDirectory("D:\\picture for test");
    fileDialog->setNameFilter("Image Files(*.jpg)");
    fileDialog->show();//打开文件对话框
    
    if (fileDialog->exec())
    {
        fileName = fileDialog->selectedFiles()[0];
    }
    if (fileName.isEmpty()) {//如果没有成功打开
        QMessageBox::warning(this, "Warning!", "Fail to open!");
        return;
    }
    myImage1.load(fileName);
    pictureopened1 = 1;
    map = QPixmap::fromImage(myImage1);

    ui.label_1->setAlignment(Qt::AlignCenter);
    ui.label_1->setPixmap(QPixmap::fromImage(myImage1).scaled(ui.label_1->size()));
    myImage1 = myImage1.scaled(ui.label_1->size());//根据画布调整大小
    myImage1_origin = myImage1;//保存一下原始图片

    undo();  //清除之前的操作

}
void possionedit::undo() {//更新完毕
    if (graphics.size() == 0)return;
    memset(mapF, 0, sizeof(mapF));
    memset(mapP, 0, sizeof(mapP));
    myImage1 = myImage1_origin;
    myImage2 = myImage2_origin;
    regionBool.clear();
    region.clear();
    NET.clear();
    for (register int i = 0; i < 405; ++i) {
        vector<edge> n;
        NET.push_back(n);
    }
    AET.clear();
    while (graphics.size() != 0) {
        PointAndLine* u = graphics.last();
        graphics.pop_back();
        free(u);
    } 
    ui.label_1->setPixmap(QPixmap::fromImage(myImage1).scaled(ui.label_1->size())); 
    ui.label_2->setPixmap(QPixmap::fromImage(myImage2).scaled(ui.label_2->size()));
    update();
    
}
void possionedit::Openpic_2()
{   //更新完毕
    QString fileName;
    QFileDialog* fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setFileMode(QFileDialog::AnyFile);
    fileDialog->setViewMode(QFileDialog::Detail);
    fileDialog->setGeometry(600, 600, 1000, 500);
    fileDialog->setDirectory("D:\\picture for test\\");
    fileDialog->setNameFilter("Image Files(*.jpg)");
    fileDialog->show();//打开文件对话框

    if (fileDialog->exec())
    {
        fileName = fileDialog->selectedFiles()[0];
    }
    if (fileName.isEmpty()) {//如果没有成功打开
        QMessageBox::warning(this, "Warning!", "Fail to open!");
        return;
    }
    myImage2.load(fileName);
    pictureopened2 = 1;
    map = QPixmap::fromImage(myImage2);

    ui.label_2->setAlignment(Qt::AlignCenter);
    ui.label_2->setPixmap(QPixmap::fromImage(myImage2).scaled(ui.label_2->size()));
    myImage2 = myImage2.scaled(ui.label_2->size());//根据画布调整大小
    myImage2_origin = myImage2;//保存一下原始图片
    undo();  //清除之前的操作
}
void possionedit::OnDraw() {//更改状态
    if(graphics.isEmpty())
    type = 1;
    else {
        QMessageBox box(QMessageBox::Warning, "警告", "已经画过了");
        box.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
        box.exec();//警告
    }
}
void possionedit::paintEvent(QPaintEvent*)
{
    if (!pictureopened1|| !pictureopened2 || graphics.isEmpty())return;
    //图像没有打开，停止绘画   //注意另一个图像有没有打开
    map = QPixmap::fromImage(myImage1).scaled(ui.label_1->size());
    QPainter p1(&map);
    PointAndLine& line = *graphics.at(0);
    for (int j = 0; j < line.path.size() - 1; j++) {
        p1.setPen(line.pen);
        p1.drawLine(line.path.at(j), line.path.at(j + 1));
    }
    ui.label_1->setPixmap(map); //更新图像
    p1.end();
    return;//更新图像
}
void possionedit::mousePressEvent(QMouseEvent* ev)
{
    if (graphics.size() >= 2)return;  //已经尘埃落定了
    if (!pictureopened1 || !pictureopened1) {
        QMessageBox box(QMessageBox::Warning, "警告", "图象还没有打开");
        box.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
        box.exec();//警告
    }
    if (ui.label_1->geometry().contains(ui.label_1->mapFromParent(ui.centralWidget->mapFromParent(ev->pos()))) && pictureopened1&& pictureopened2)
    {
        QPoint s = ui.label_1->mapFromParent(ui.centralWidget->mapFromParent(ev->pos())); //在画布中的位置
        if (mapF[s.x()][s.y()] != 0) {
            OnMoving = 1;  //对目标区域进行拖拽
            origin_pos = s;
            PointAndLine* r = new PointAndLine; //初始化
            r->type = 2;
            pen_pub.setColor(QColor(0, 255, 0));
            pen_pub.setWidth(1);
            r->pen = pen_pub;
            r->path.append(s);   //在image2中对应的外方框
            r->path.append(s);
            r->path.append(s);
            r->path.append(s);
            r->path.append(s);
            graphics.append(r);
        }

        if (type == 0)return;//没选择
        if (type == 1) {
            PointAndLine* r = new PointAndLine; //初始化
            r->type = 1;
            pen_pub.setColor(QColor(0, 0, 255));
            pen_pub.setWidth(3);
            r->pen = pen_pub;
            r->path.append(s);
            graphics.append(r);
        }
    }
    update();
}
void possionedit::mouseMoveEvent(QMouseEvent* ev)//绘制操作点和操作线
{
    if (pictureopened1 && pictureopened2 && !graphics.isEmpty()) {
        if (OnMoving && ui.label_2->geometry().contains(this->mapFromGlobal(QCursor::pos()))) {//正在拖拽目标区域
            PointAndLine* lastLine = graphics.last();
            int type_1 = (*graphics.last()).type;
            if (type_1 == 2) {
                lastLine->draw(ui.label_2->mapFromParent(ui.centralWidget->mapFromParent(ev->pos())) - origin_pos + left_up, ui.label_2->mapFromParent(ui.centralWidget->mapFromParent(ev->pos())) - origin_pos + right_down);
            }
            possion();  //实时计算
            update();
        }
        if (ui.label_1->geometry().contains(ui.label_1->mapFromParent(ui.centralWidget->mapFromParent(ev->pos()))) && !OnMoving)//正在绘制目标区域
        {
            if (type == 0)return;
            PointAndLine* lastLine = graphics.last();
            int type_1 = (*graphics.last()).type;
            if (type_1 == 1) {//自由绘图
                QPoint s = lastLine->path.last();
                QPoint t = ui.label_1->mapFromParent(ui.centralWidget->mapFromParent(ev->pos()));
                QPoint v;
                if (abs(s.x() - t.x()) + abs(s.y() - t.y()) > 6) {//避免重复加入
                    lastLine->draw(t);
                    if (s.y() != t.y()) { //不能与扫描线平行
                        if (lastLine->path.size() > 1) {//更新顶点
                            QPoint r = lastLine->path.at(lastLine->path.size() - 2);
                            int flag = (s.y() - r.y()) * (s.y() - t.y());//判断是不是在同一侧
                            if (flag > 0) {
                                mapP[s.x()][s.y()] = 1;   //同一侧 ， 虚值
                            }
                            else if (flag < 0) {
                                mapP[s.x()][s.y()] = 0;    //异侧 ， 实值
                            }
                            else {
                                mapP[s.x()][s.y()] = 0;   //平行 ， 实值
                            }
                        }
                        if (s.y() > t.y()) {
                            v = s;
                            s = t;
                            t = v;
                        }
                        eg.current = s.x();
                        eg.delta_x = (double)(t.x() - s.x()) / (double)(t.y() - s.y());
                        eg.ymax = t.y();
                        NET.at(s.y()).push_back(eg); //添加NET
                    }
                    else {//对顶点进行更新
                        mapP[s.x()][s.y()] = 0; //平行时不记录边，顶点记为实值
                    }

                }
            }
            update();
        }
    }
}

void possionedit::scan() {//扫描线算法
    bool isempty = 1;
    short isinterior = -1;
    int x = 0;
    int y = 0;
    int k = 0;
    int linecut[405];
    for (register int i = 0; i < 405; ++i) {//扫描每一行
        isinterior = -1;
        k = 0;
        memset(linecut, 0, sizeof(linecut));
        for (register int j = 0; j < NET.at(i).size(); ++j) {
            AET.push_back(NET.at(i).at(j));   //将NET中从i出发的那些添加到AET中
        }
        for (register int j = 0; j < AET.size(); ++j) {
            if (AET.at(j).ymax == i) {//如果边到了最值
                AET.erase(AET.begin() + j); //删除
                --j;
            }
            else {
                int b = (int)(AET.at(j).current + 0.5);//四舍五入
                if (mapP[b][i] == 0) {
                    linecut[k] = b;
                    k++;
                }
                AET.at(j).current += AET.at(j).delta_x;  //每一条相交的边加delta
            }
        }
        qsort(linecut, k, sizeof(int), cmp);  //对交点进行排序
        int u = linecut[0];
        for (register int j = 0; j < k-1; ++j) {//填充
            isinterior = -isinterior;
            while (u != linecut[j + 1]) {
                if(isinterior == 1)mapF[u][i] = 1;
                u++;
            }
            
        }
    }
    //下面是确认边框
    while (isempty) {
        for (register int j = 0; j < 405; ++j) {
            if (mapF[x][j] != 0)isempty = 0;
        }
        x++;
    }
    left_up.setX(x-1);
    isempty = 1;
    while (isempty) {
        for (register int j = 0; j < 405; ++j) {
            if (mapF[j][y] != 0)isempty = 0;
        }
        y++;
    }
    left_up.setY(y - 1);
    isempty = 1;
    x = 404;
    while (isempty) {
        for (register int j = 0; j < 405; ++j) {
            if (mapF[x][j] != 0)isempty = 0;
        }
        x--;
    }
    right_down.setX(x + 1);
    isempty = 1;
    y = 404;
    while (isempty) {
        for (register int j = 0; j < 405; ++j) {
            if (mapF[j][y] != 0)isempty = 0;
        }
        y--;
    }
    right_down.setY(y + 1);

    for (register int i = 0; i < 405; ++i) { //绘制边界
        for (register int j = 0; j < 405; ++j) {
            if (mapF[i][j] == 1) {
                if ((i > 0 && mapF[i - 1][j] == 0) || (i<404&& mapF[i+1][j]==0) || (j > 0 && mapF[i][j-1] == 0) || (j < 404 && mapF[i ][j+1] == 0)) {
                    mapF[i][j] = 2;                   
                }
            }
        }
    }

}
void possionedit::mouseReleaseEvent(QMouseEvent* ev) {
    if (graphics.size() == 1 && type == 1) {
        PointAndLine* lastLine = graphics.last();
        QPoint s = lastLine->path.last();
        QPoint t = lastLine->path.at(0);
        lastLine->draw(t);
        QPoint v;
        if (s.y() != t.y()) { //不能与扫描线平行
            if (lastLine->path.size() > 1) {//更新顶点
                QPoint r = lastLine->path.at(lastLine->path.size() - 2);
                int flag = (s.y() - r.y()) * (s.y() - t.y());//判断是不是在同一侧
                if (flag > 0) {
                    mapP[s.x()][s.y()] = 1;   //同一侧 ， 虚值
                }
                else if (flag < 0) {
                    mapP[s.x()][s.y()] = 0;    //异侧 ， 实值
                }
                else {
                    mapP[s.x()][s.y()] = 0;   //平行 ， 实值
                }
            }
            if (s.y() > t.y()) {
                v = s;
                s = t;
                t = v;
            }
            eg.current = s.x();
            eg.delta_x = (double)(t.x() - s.x()) / (double)(t.y() - s.y());
            eg.ymax = t.y();
            NET.at(s.y()).push_back(eg); //添加NET
        }
        else {//对顶点进行更新
            mapP[s.x()][s.y()] = 0; //平行时不记录边，顶点记为实值
            //if (lastLine->path.size() > 1)QPoint r = lastLine->path.at(lastLine->path.size() - 2);
        }
       

        s = lastLine->path.last();  //对最初顶点的判断
        t = lastLine->path.at(1);
        QPoint r = lastLine->path.at(lastLine->path.size() - 2);
        int flag = (s.y() - r.y()) * (s.y() - t.y());//判断是不是在同一侧
        if (flag > 0) {
            mapP[s.x()][s.y()] = 1;   //同一侧 ， 虚值
        }
        else if (flag < 0) {
            mapP[s.x()][s.y()] = 0;    //异侧 ， 实值
        }
        else {
            mapP[s.x()][s.y()] = 0;   //平行 ， 实值
        }

        scan();
        possion_pre();
        graphics.at(0)->pen.setColor(qRgb(0, 255, 0));
    }
    type = 0;
    OnMoving = 0;
    
}
void possionedit::possion_pre() {

    QVector<short> ra, ra1 ;
    QVector<long>  rb, rb1;
    long k = 0;
    //初始化region 和 regionBool
    for (register int i = 0; i <= right_down.x()-left_up.x(); ++i) {
        ra1.append(0);
        rb1.append(0);
    }
    for (register int i = 0; i <= right_down.y() - left_up.y(); ++i) {
        ra = ra1;
        rb = rb1;
        for (register int j = 0; j <= right_down.x() - left_up.x(); ++j) {
            if (mapF[left_up.x() + j][left_up.y() + i] != 0) {
                ra.replace(j, mapF[left_up.x() + j][left_up.y() + i]);
                rb.replace(j, k);
                k++;
            }
           
        }
        regionBool.append(ra);
        region.append(rb);
    }
    image_part1 = myImage1_origin.copy(left_up.x(), left_up.y(), right_down.x() - left_up.x()+1, right_down.y() - left_up.y()+1);
    e.setregionBool(regionBool);
    e.setregion(region);
    e.setPosnum(k);
    e.initial();
    num = k;
    p.setimage1(image_part1);
    p.setMatrix(e.getmatrix());
    p.setregionBool(regionBool);
    p.setnum(num);
    p.initial();
}
void possionedit::possion() {
    image_part2 = myImage2_origin.copy((*(graphics.last())).path.at(0).x(), (*(graphics.last())).path.at(0).y(), abs((*(graphics.last())).path.at(2).x() - (*(graphics.last())).path.at(0).x())+1, abs((*(graphics.last())).path.at(2).y() - (*(graphics.last())).path.at(0).y())+1);
    p.setimage2(image_part2);
    p.analysis(); //计算
    myImage2 = myImage2_origin;  //获取计算出来的图像
    QPainter imagepainter(&myImage2);  //新建画板
    imagepainter.setCompositionMode(QPainter::CompositionMode_SourceOver);   //设置重叠效果
    imagepainter.drawImage((*(graphics.last())).path.at(0).x(), (*(graphics.last())).path.at(0).y(), p.getimage());
    imagepainter.end();
    ui.label_2->setPixmap(QPixmap::fromImage(myImage2).scaled(ui.label_2->size()));//更新

}

