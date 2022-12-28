#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vector>
#include <QPainter>
#include <QMouseEvent>
#include <figure.h>
using namespace std;
QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
protected:
    void paintLine();
    void paintEvent(QPaintEvent *ev);
    void mouseMoveEvent(QMouseEvent *ev);
    void mousePressEvent(QMouseEvent *ev);
    //void mouseReleaseEvent(QMouseEvent *ev);

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QVector<CFigure*> graphics;  //存储的各类图形
    int type=0;  // 当前图形的类型
    QPen pen_pub;  //当前笔刷
private slots:   //对应的UI事件
    void on_actionrectangle_triggered();

    void on_actionelliptical_triggered();

    void on_actionline_triggered();

    void on_actionfree_line_triggered();

    void on_actionpolygon_triggered();

    void on_actionundo_triggered();

    void on_actioncolor_triggered();

private:
    Ui::MainWindow *ui;
};

#endif

