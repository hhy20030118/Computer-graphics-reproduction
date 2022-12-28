#include "ImageWrapping.h"
#include "qdialog.h"
#include "QFileDialog"
#include <QMessageBox>
#include <QImage>
#include "paintwidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <QDesktopWidget>
#include <QApplication>
#include "ImageOperator.h"
#include "IDWarg.h"
#include "RBFarg.h"
#include<iostream>
ImageWrapping::ImageWrapping(QWidget *parent): QMainWindow(parent)
{
    ui.setupUi(this);
    QMenu* MenuA = menuBar()->addMenu(QString::fromLocal8Bit("File"));
    MenuA->addAction(ui.actionopen);
    //ui.mainToolBar->addAction(new_file);
    //connect(new_file, SIGNAL(triggered()), this, SLOT(OpenFile()));
    connect(ui.actionopen, SIGNAL(triggered()), this, SLOT(OpenFile()));  //设置触发函数
    connect(ui.actiondrawpoint, SIGNAL(triggered()), this, SLOT(ChangeTypeToPoint()));
    connect(ui.actiondrawline, SIGNAL(triggered()), this, SLOT(ChangeTypeToLine()));
    connect(ui.actionundo, SIGNAL(triggered()), this, SLOT(undo()));
    connect(ui.actionIDW, SIGNAL(triggered()), this, SLOT(IDWoperate()));
    connect(ui.actionRBF, SIGNAL(triggered()), this, SLOT(RBFoperate()));
    connect(ui.actionrenew, SIGNAL(triggered()), this, SLOT(renew()));
    connect(ui.actionsave, SIGNAL(triggered()), this, SLOT(save()));
    //ui.actionopen->setIcon(QIcon("\\resources\\1.png"));
    this->resize(QSize(1300, 1100));  //设置窗口大小
    pen_pub.setWidth(10);           //初始化设置
    pen_pub.setColor(Qt::blue);
    ui.Basicanvas->installEventFilter(this);
    ui.mainToolBar->setMovable(false);
    //SystemParametersInfo(SPI_GETWORKAREA, 0, &rtWorkArea, 0);
    //paintwidget pw;
    //pw.graphics = &this->graphics;
    //pw.show();
}
QPoint delta(0,75);   //修正，因为工具栏引起的误差
void ImageWrapping::paintEvent(QPaintEvent *)
{
    if (!pictureopened)return;   //图像没有打开，停止绘画
    map = QPixmap::fromImage(myImage).scaled(ui.Basicanvas->size());
    QPainter p(&map);
    for (int i = 0; i < graphics.size(); i++)//画操作点和操作线
    {
        const PointAndLine& line = *graphics.at(i);     
            for (int j = 0; j < line.path.size() - 1; j++) {
                p.setPen(line.pen);
                p.drawLine(ui.Basicanvas->mapFromParent(line.path.at(j))- delta,  ui.Basicanvas->mapFromParent(line.path.at(j+1))- delta);
            }
        
    }
     ui.Basicanvas->setPixmap(map); //更新图像
}
/*bool ImageWrapping::eventFilter(QObject* watched, QEvent* event)   //用过滤器eventFilter（）拦截QLabel中的QEvent::Paint事件
{
    if (watched == ui.Basicanvas && event->type() == QEvent::Paint)
        labelPaint();

    return QWidget::eventFilter(watched, event);
}

void ImageWrapping::labelPaint()     //绘图
{
    QPainter p(this);
    p.setPen(Qt::blue);
    // painter.drawLine(100,100,200,200);
    p.drawEllipse(30, 15, 50, 65);
    p.drawLine(0, 100, 111, 100);
    for (int i = 0; i < graphics.size(); i++)
    {
        const PointAndLine& line = *graphics.at(i);
        for (int j = 0; j < line.path.size() - 1; j++) {
            p.setPen(line.pen);
            p.drawLine(line.path.at(j), line.path.at(j + 1));
        }

    }

}*/
double max(double a , double b) {
    return a > b ? a : b;
}  //最大值
QPoint po;
short invovled[1000][1000];
inline int correct(int a , int b , int c) {
    if (c < a)return a;
    if (c > b)return b;
    else return c;
}//防止越界
void ImageWrapping::IDWoperate() {//IDW算法
    IDWarg IDW;
    if (!pictureopened)return;
    QImage img(myImage.height(), myImage.width(), QImage::Format_RGB32);
    img.fill(QColor(Qt::white));  //初始化，传入
    QRgb* pixs = (QRgb*)img.bits();
    img = img.scaled(ui.Basicanvas->size());
    IDW.ChangedPos = o.ChangedPos;
    int W = ui.Basicanvas->width();
    int H = ui.Basicanvas->height();
    for (int i = 0; i < H; ++i) {//初始化，invovled存的是是否被涂色，方便找白缝
        for (int j = 0; j < W; ++j) {
            invovled[i][j] = 0;
        }
    }
    clock_t start, end;
    start = clock();
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            //img.setPixel(i, j, qRgb(55 , 55 , 55));

            po = QPoint((int)IDW.getpos(j, i).x(), (int)IDW.getpos(j, i).y());
            invovled[correct(0, H, po.y())][correct(0, W, po.x())] = 1;//通过计算算出来，然后对应点涂色
            img.setPixel(correct(0, W, po.x()), correct(0, H, po.y()), qRgb(qRed(myImage.pixel(j, i)), qGreen(myImage.pixel(j, i)), qBlue(myImage.pixel(j, i))));
        }

    }
    for (int i = 0; i < H; ++i) {//处理边界，以免把空白部分当成白缝处理
        int j = W - 1;
        while (invovled[i][j] != 1 && j >= 0) {

            invovled[i][j] = 3;
            j--;
        }
        j = 0;
        while (invovled[i][j] != 1 && j < W - 1) {

            invovled[i][j] = 3;
            j++;
        }
    }
    for (int i = 0; i < W; ++i) {
        int j = H - 1;
        while (invovled[j][i] != 1 && j >= 0) {
            invovled[j][i] = 3;
            j--;
        }
        j = 0;
        while (invovled[j][i] != 1 && j < H - 1) {
            invovled[j][i] = 3;
            j++;
        }
    }
    int tr = 4;
    /*int R, G, B;
    int op = 0;
    for (int s = 0 ; s <= 2 ; s++)
        for (int i = tr; i < H - tr; ++i) {
            for (int j = tr; j < W - tr; ++j) {
                if (invovled[i][j] == 0) {
                    op = 0;
                    for (int k = i - tr; k < i + tr; ++k) {
                        for (int l = j - tr; l < j + tr; ++l) {
                            if (invovled[k][l] == 1) {
                            R += qRed(img.pixel(l, k));//取平均算白缝
                            G += qGreen(img.pixel(l, k));
                            B += qBlue(img.pixel(l, k));
                            op++;
                            }

                        }
                        if(op!=0)img.setPixel(j, i, qRgb(R/op,G /op, B /op));
                    }

                }
            }
        }*/
        //找附近点的颜色去白缝
    for (int s = 0; s <= 2; s++)
        for (int i = tr; i < H - tr; ++i) {
            for (int j = tr; j < W - tr; ++j) {
                if (invovled[i][j] == 0) {
                    bool flag = 1;
                    int d = 1;
                    while (flag) {
                        if (d >= j || d >= i)break;
                        if (invovled[i + d][j + d] == 1) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i + d)), qGreen(img.pixel(j + d, i + d)), qBlue(img.pixel(j + d, i + d))));
                            flag = 0;
                        }
                        else if (invovled[i - d][j - d] == 1) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j - d, i - d)), qGreen(img.pixel(j - d, i - d)), qBlue(img.pixel(j - d, i - d))));
                            flag = 0;
                        }
                        else if (invovled[i - d][j + d] == 1) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i - d)), qGreen(img.pixel(j + d, i - d)), qBlue(img.pixel(j + d, i - d))));
                            flag = 0;
                        }
                        else if (invovled[i + d][j - d] == 1) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j - d, i + d)), qGreen(img.pixel(j - d, i + d)), qBlue(img.pixel(j - d, i + d))));
                            flag = 0;
                        }
                        else d++;
                    }

                }
            }
        }

        /*for (int d = -tr; d <= tr; ++d) {
                        if(!flag)break;
                        for (int f = -tr; f <= tr; ++f) {
                              if (!flag)break;
                              if (img.pixel(j + d, i + f) != Qt::white && (qRed(myImage.pixel(j + d, i + f)) + qGreen(myImage.pixel(j + d, i + f)) + qBlue(myImage.pixel(j + d, i + f))) < 600) {
                              img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i + f)), qGreen(img.pixel(j + d, i + f)), qBlue(img.pixel(j + d, i + f))));
                              flag = 0;
                              }
                        }


                }
        }*/ 
    //统计时间
    end = clock();
    double endtime = (double)(end - start)*1000 / CLOCKS_PER_SEC;
    QString str = QString::number(endtime, 'f', 2);
    QMessageBox m(this);
    m.setWindowTitle("MyAction1");
    m.setText(str+"ms");
    m.exec();   //提示，可去
    ui.Basicanvas->clear(); //去除之前的图片
    myImage = img;  //最后将img传入myimg
}
void ImageWrapping::RBFoperate() {
    RBFarg RBF;
    if (!pictureopened)return;
    QImage img(myImage.height(), myImage.width(), QImage::Format_RGB32);
    img.fill(QColor(Qt::white));
    QRgb* pixs = (QRgb*)img.bits();
    img = img.scaled(ui.Basicanvas->size());
    RBF.ChangedPos = o.ChangedPos;
    RBF.initial();
    int W = ui.Basicanvas->width();
    int H = ui.Basicanvas->height();
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            invovled[i][j] = 0;
        }
    }
    clock_t start, end;
    start = clock();
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            //img.setPixel(i, j, qRgb(55 , 55 , 55));

            po = QPoint((int)RBF.getpos(j, i).x(), (int)RBF.getpos(j, i).y());
            invovled[correct(0, H, po.y())][correct(0, W, po.x())] = 1;
            img.setPixel(correct(0, W, po.x()), correct(0, H, po.y()), qRgb(qRed(myImage.pixel(j, i)), qGreen(myImage.pixel(j, i)), qBlue(myImage.pixel(j, i))));
        }

    }
    for (int i = 0; i < H; ++i) {
        int j = W - 1;
        while (invovled[i][j] != 1 && j >= 0) {

            invovled[i][j] = 3;
            j--;
        }
        j = 0;
        while (invovled[i][j] != 1 && j < W - 1) {

            invovled[i][j] = 3;
            j++;
        }
    }
    for (int i = 0; i < W; ++i) {
        int j = H - 1;
        while (invovled[j][i] != 1 && j >= 0) {
            invovled[j][i] = 3;
            j--;
        }
        j = 0;
        while (invovled[j][i] != 1 && j < H - 1) {
            invovled[j][i] = 3;
            j++;
        }
    }
    int tr = 3;
    for (int s = 0; s <= 2; s++)
        for (int i = tr; i < H - tr; ++i) {
            for (int j = tr; j < W - tr; ++j) {
                if (invovled[i][j] == 0) {
                    //
                    bool flag = 1;
                    int d = 1;
                    while (flag) {
                        if (d >= j || d >= i)break;
                        if (img.pixel(j + d, i + d) != Qt::white && (qRed(img.pixel(j + d, i + d)) + qGreen(img.pixel(j + d, i + d)) + qBlue(img.pixel(j + d, i + d))) < 700) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i + d)), qGreen(img.pixel(j + d, i + d)), qBlue(img.pixel(j + d, i + d))));
                            flag = 0;
                        }
                        else if (img.pixel(j - d, i - d) != Qt::white && (qRed(img.pixel(j - d, i - d)) + qGreen(img.pixel(j - d, i - d)) + qBlue(img.pixel(j - d, i - d))) < 700) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j - d, i - d)), qGreen(img.pixel(j - d, i - d)), qBlue(img.pixel(j - d, i - d))));
                            flag = 0;
                        }
                        else if (img.pixel(j + d, i - d) != Qt::white && (qRed(img.pixel(j + d, i - d)) + qGreen(img.pixel(j + d, i - d)) + qBlue(img.pixel(j + d, i - d))) < 700) {
                            img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i - d)), qGreen(img.pixel(j + d, i - d)), qBlue(img.pixel(j + d, i - d))));
                            flag = 0;
                        }
                        else d++;
                    }
                    /*for (int d = -tr; d <= tr; ++d) {
                            if(!flag)break;
                            for (int f = -tr; f <= tr; ++f) {
                                  if (!flag)break;
                                  if (img.pixel(j + d, i + f) != Qt::white && (qRed(myImage.pixel(j + d, i + f)) + qGreen(myImage.pixel(j + d, i + f)) + qBlue(myImage.pixel(j + d, i + f))) < 600) {
                                  img.setPixel(j, i, qRgb(qRed(img.pixel(j + d, i + f)), qGreen(img.pixel(j + d, i + f)), qBlue(img.pixel(j + d, i + f))));
                                  flag = 0;
                                  }
                            }


                    }*/
                    // if(flag)img.setPixel(j, i, qRgb(55 , 66 , 230));  
                     //else img.setPixel(j, i, qRgb(255, 66, 47));
                }
            }
        }
    end = clock();
    double endtime = (double)(end - start) * 1000 / CLOCKS_PER_SEC;
    QString str = QString::number(endtime, 'f', 2);
    QMessageBox m(this);
    m.setWindowTitle("MyAction1");
    m.setText(str + "ms");
    m.exec();
    ui.Basicanvas->clear();
    myImage = img;

    //QPixmap u = QPixmap::fromImage(img);
    //ui.Basicanvas->setPixmap(u);
}
void ImageWrapping::ChangeTypeToPoint() {
    type = 1;//调整操作模式
}
void ImageWrapping::ChangeTypeToLine() {
    type = 2;//调整操作模式
}
void ImageWrapping::undo() {
    if (graphics.size() == 0)return;
    if (o.ChangedPos.size() == 0)return;//防止越界
    PointAndLine* u = graphics.last();
    pairpos* v = &o.ChangedPos.last();//两个同时退，而且因为操作线包括点和线，要去两次
    graphics.pop_back();
    if(u->type != 3) o.ChangedPos.pop_back();
    free(u);
    update();
}
void ImageWrapping::renew() {
    myImage = myImage_origin;//重置图片
    while (graphics.size() != 0) {
        undo();
    }
}
void ImageWrapping::OpenFile() {
    QFileDialog* fileDialog = new QFileDialog(this);
    fileDialog->setWindowTitle(tr("Open File"));
    fileDialog->setFileMode(QFileDialog::AnyFile);
    fileDialog->setViewMode(QFileDialog::Detail);
    fileDialog->setGeometry(600, 600, 1000, 500);
    fileDialog->setDirectory("D:");
    fileDialog->setNameFilter("Image Files(*.jpg)");
    fileDialog->show();//打开文件对话框
    QString fileName;
    if (fileDialog->exec())
    {
        fileName = fileDialog->selectedFiles()[0];
    }
    if (fileName.isEmpty()) {//如果没有成功打开
        QMessageBox::warning(this, "Warning!", "Fail to open!");
        return;
    }  
    myImage.load(fileName);
    pictureopened = 1;
    map = QPixmap::fromImage(myImage);
    
    ui.Basicanvas->resize(myImage.size()/max(myImage.size().height()/900.0, myImage.size().width()/900.0));
    ui.Basicanvas->setAlignment(Qt::AlignCenter);
    ui.Basicanvas->setPixmap(QPixmap::fromImage(myImage).scaled(ui.Basicanvas->size()));
    myImage = myImage.scaled(ui.Basicanvas->size());//根据画布调整大小
    myImage_origin = myImage;//保存一下原始图片
    ImageOperator basic(ui.Basicanvas->height(), ui.Basicanvas->width());
    o = basic;
    while (graphics.size() != 0) {
        undo();
    }//清空数据
}
void ImageWrapping::save() {//保存图片
    if (pictureopened) {
        QString filename1 = QFileDialog::getSaveFileName(this, tr("Save Image"), "", tr("Images (*.png *.bmp *.jpg)")); //选择路径
        myImage.save(filename1);
    }   
}
void ImageWrapping::mousePressEvent(QMouseEvent* ev)
{

    if (ui.Basicanvas->geometry().contains(this->mapFromGlobal(QCursor::pos()-delta))&&pictureopened)
    {
        if (type == 0)return;//没选择
        if (type == 1) {//点，直接加
            PointAndLine* r = new PointAndLine;
            r->type = 1;
            pen_pub.setColor(QColor(0, 0, 255));
            r->pen = pen_pub;
            r->path.append(ev->pos());
            graphics.append(r);
            pairpos l;
            l.bpos = ui.Basicanvas->mapFromParent(ev->pos()) - QPoint(0, 78);
            o.ChangedPos.append(l);
        }
        else if (type == 2) {//线，但是要加上一个红点，所以有两个
            PointAndLine* r = new PointAndLine;
            r->type = 2;
            pen_pub.setColor(QColor(0, 0, 255));
            r->pen = pen_pub;
            r->path.append(ev->pos());
            pairpos l;
            l.bpos = ui.Basicanvas->mapFromParent(ev->pos()) - QPoint(0, 78);
            o.ChangedPos.append(l);
            PointAndLine* i = new PointAndLine;
            i->type = 3;
            pen_pub.setColor(QColor(255 , 0 , 0));
            i->pen = pen_pub;
            i->path.append(ev->pos());
            graphics.append(r);
            graphics.append(i);

        }
        QVector<QPoint>& lastLine = (*graphics.last()).path;
        lastLine.append(ev->pos());
        if(o.ChangedPos.size()!=0)o.ChangedPos.last().epos = ui.Basicanvas->mapFromParent(ev->pos()) - QPoint(0, 80);
    }
        update();
}
void ImageWrapping::mouseMoveEvent(QMouseEvent* ev)//绘制操作点和操作线
{
    if (ui.Basicanvas->geometry().contains(this->mapFromGlobal(QCursor::pos() - QPoint(0, 78))) && pictureopened)
    {
        if (type == 0)return;
        QVector<QPoint>& lastLine = (*graphics.last()).path;
        int type_1 = (*graphics.last()).type;
       
       if (type_1 == 2) {
            while (lastLine.size() < 2 ) {
                QPoint p;
                lastLine.append(p);
            }
            lastLine.replace(1, ev->pos());

       }
       else if (type_1 == 3) {
           QVector<QPoint>& lastLine1 = graphics.at(graphics.size()-2)->path;
           while (lastLine1.size() < 2) {
               QPoint p;
               lastLine1.append(p);
           }
           lastLine1.replace(1, ev->pos());
        }
       if(type_1 != 1)if (o.ChangedPos.size() != 0)o.ChangedPos.last().epos = ui.Basicanvas->mapFromParent(ev->pos()) - QPoint(0, 78);
        update();
    }
}
