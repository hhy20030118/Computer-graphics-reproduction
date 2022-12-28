#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPainter>
#include <QMouseEvent>
#include <mainwindow1.h>
#include "math.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->centralwidget->setMouseTracking(true);
    ui->toolBar->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
    /*ui->actionrectangle->setIcon(QIcon("/resources/1.png"));
    ui->actionelliptical->setIcon(QIcon("/resources/2.webp"));
    ui->actionline->setIcon(QIcon(":/resources/3.webp"));
    ui->actionfree_line->setIcon(QIcon(":resources/4.webp"));
    ui->actionpolygon->setIcon(QIcon(":resources/5.webp"));
    ui->actionundo->setIcon(QIcon(":resources/6.webp"));
    ui->actioncolor->setIcon(QIcon(":resources/7.webp"));
    ui->toolBar->addAction(ui->actionrectangle);
    ui->toolBar->addAction(ui->actionelliptical);
    ui->toolBar->addAction(ui->actionline);
    ui->toolBar->addAction(ui->actionfree_line);
    ui->toolBar->addAction(ui->actionpolygon);
    ui->toolBar->addAction(ui->actionundo);
    ui->toolBar->addAction(ui->actioncolor);*/
}

void MainWindow::paintEvent(QPaintEvent *)
{
    QPainter p(this);   //在主窗体上画图
    for(int i=0; i<graphics.size(); i++)
    {
        CFigure &line = *graphics.at(i); //其中第i个
        for(int j=line.gettype()==2?1:0; j<line.path.size()-1; j++){
            p.setPen(line.getpen());   //设置笔刷
            p.drawLine(line.path.at(j), line.path.at(j+1)); //连线
        }

    }
}

void MainWindow::mouseMoveEvent(QMouseEvent *ev)
{
    if(type==0)return;  //还没选
    (*graphics.last()).draw(ev->pos());  //调用对应的绘画函数
    update();
}

void MainWindow::mousePressEvent(QMouseEvent *ev)
{
    if(type==0)return; //还没选
    if (ev->button() == Qt::RightButton)  //右键
    {
        if(type==4&&graphics.size()!=0){
            setMouseTracking(false);  //禁止不按时追踪
            ((CPolypon*)(graphics.last()))->end();  //调用结束函数 ， 画完了
            update();
        }
    }else{
    if(type==1){
        CRect* r=new CRect;  //添加新图形
        r->settype(1);
        r->setpen(pen_pub);
        graphics.append(r);

    }
    else if(type==2){
        CEllipse* r=new CEllipse;
        r->settype(2);
        r->setpen(pen_pub);
        graphics.append(r);
    }else if(type==3){
        CLine* r=new CLine;
        r->settype(3);
        r->setpen(pen_pub);
        graphics.append(r);

    }else if(type==4){
        CPolypon* r=new CPolypon;
        r->settype(4);
        r->setpen(pen_pub);
        r->setoprate(true); //正在操作
        if(graphics.size()==0||((*graphics.last()).gettype()!=4)||!((CPolypon*)(graphics.last()))->getoprate())
        graphics.append(r);

        setMouseTracking(true); //启用实时跟踪
        QVector<QPoint> &lastLine = (*graphics.last()).path;
        lastLine.append(ev->pos()); //先加一个点
    }else if(type==5){
        CFreehand* r=new CFreehand;
        r->settype(5);
        r->setpen(pen_pub);
        graphics.append(r);
        }
    QVector<QPoint> &lastLine = (*graphics.last()).path;
    lastLine.append(ev->pos());

    }
}


MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionrectangle_triggered()
{
    if(graphics.size() == 0 || !(graphics.last()->gettype() == 4 && ((CPolypon*)(graphics.last()))->getoprate()==true))  //防止多边形中途换类
    type=1;       
}

void MainWindow::on_actionelliptical_triggered()
{
    if(graphics.size() == 0 || !(graphics.last()->gettype() == 4 && ((CPolypon*)(graphics.last()))->getoprate()==true))
    type=2;
}

void MainWindow::on_actionline_triggered()
{
    if(graphics.size()==0||!(graphics.last()->gettype() == 4 && ((CPolypon*)(graphics.last()))->getoprate()==true))
    type=3;
}

void MainWindow::on_actionpolygon_triggered()
{

    type=4;
}

void MainWindow::on_actionfree_line_triggered()
{
    if(graphics.size() == 0 ||!( graphics.last()->gettype() == 4 && ((CPolypon*)(graphics.last()))->getoprate()==true))
    type=5;
}
void MainWindow::on_actionundo_triggered()
{
    if(graphics.size()==0)return;
    CFigure* u=graphics.last();//撤销 ， 去除最后一个
    graphics.pop_back();
    free(u);
    update();
}


void MainWindow::on_actioncolor_triggered()
{
    MainWindow1* w = new MainWindow1();
    w->win = this;    //弹出窗口设置笔刷
    w->R1=pen_pub.color().red();
    w->G1=pen_pub.color().green();
    w->B1=pen_pub.color().blue();
    w->width1=pen_pub.width();
    w->change();
    w->show();
}

