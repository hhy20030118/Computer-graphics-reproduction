#include "possionedit.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    possionedit w;
    w.show();
    return a.exec();
}
