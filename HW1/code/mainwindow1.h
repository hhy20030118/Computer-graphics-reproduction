#ifndef MAINWINDOW1_H
#define MAINWINDOW1_H

#include <QMainWindow>
#include"mainwindow.h"
namespace Ui {
class MainWindow1;
}

class MainWindow1 : public QMainWindow
{
    Q_OBJECT
protected:
    void paintEvent(QPaintEvent *ev);
public:
    explicit MainWindow1(QWidget *parent = nullptr);
    ~MainWindow1();
    int R1 , G1 , B1 ;  //笔刷颜色
    MainWindow* win;  //存储母窗体
    int width1=10;   //笔刷粗细
    QPen pen;
    void change();
private slots:
    void on_pushButton_2_clicked();

private:
    Ui::MainWindow1 *ui;
};

#endif // MAINWINDOW1_H
