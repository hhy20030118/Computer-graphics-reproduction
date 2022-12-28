#include "mainwindow.h"

#include <QApplication>
#include <vector>
#include <mainwindow1.h>
#define ADD_NUM 100
using namespace std;

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    MainWindow1 r;

    return a.exec();
}
