#include "mainwindow1.h"
#include "ui_mainwindow1.h"
#include "mainwindow.h"
MainWindow1::MainWindow1(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow1)
{
    ui->setupUi(this);
    this->pen=QPen(QColor(0,0,0,255));

}
void MainWindow1::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    this->pen.setWidth(20);   //颜色样例 ， 提供参考
    p.setPen(this->pen);
    p.drawPoint(25,25);
    update();
}
MainWindow1::~MainWindow1()
{
    delete ui;
}
void MainWindow1::change()
{
    ui->spinBox->setValue(R1);   //初始化 ， 导入当前数据
    ui->spinBox_2->setValue(G1);
    ui->spinBox_3->setValue(B1);
    ui->spinBox_4->setValue(width1);
    pen.setColor(QColor(R1 , G1 , B1));
}


void MainWindow1::on_pushButton_2_clicked()
{
   this->pen.setWidth(ui->spinBox_4->value());   //确认后调整主窗口笔刷
   this->pen.setColor(QColor(ui->spinBox->value(),ui->spinBox_2->value(),ui->spinBox_3->value(),255));
   win->pen_pub=this->pen;
}







