#include "ImageWrapping.h"
#include <QtWidgets/QApplication>
#include <Eigen/Dense>  
#include<iostream>
     using namespace std;
     using namespace Eigen;
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ImageWrapping w;
    w.show();
    return a.exec();
}
