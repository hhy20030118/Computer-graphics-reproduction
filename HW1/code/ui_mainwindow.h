/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.12.12
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QIcon>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionrectangle;
    QAction *actionelliptical;
    QAction *actionline;
    QAction *actionfree_line;
    QAction *actionpolygon;
    QAction *actionundo;
    QAction *actioncolor;
    QWidget *centralwidget;
    QStatusBar *statusbar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(800, 600);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/new/prefix1/resouses/7.webp"), QSize(), QIcon::Normal, QIcon::Off);
        MainWindow->setWindowIcon(icon);
        actionrectangle = new QAction(MainWindow);
        actionrectangle->setObjectName(QString::fromUtf8("actionrectangle"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/new/prefix1/resouses/1.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionrectangle->setIcon(icon1);
        actionelliptical = new QAction(MainWindow);
        actionelliptical->setObjectName(QString::fromUtf8("actionelliptical"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/new/prefix1/resouses/2.webp"), QSize(), QIcon::Normal, QIcon::Off);
        actionelliptical->setIcon(icon2);
        actionline = new QAction(MainWindow);
        actionline->setObjectName(QString::fromUtf8("actionline"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/new/prefix1/resouses/3.webp"), QSize(), QIcon::Normal, QIcon::Off);
        actionline->setIcon(icon3);
        actionfree_line = new QAction(MainWindow);
        actionfree_line->setObjectName(QString::fromUtf8("actionfree_line"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/new/prefix1/resouses/4.webp"), QSize(), QIcon::Normal, QIcon::Off);
        actionfree_line->setIcon(icon4);
        actionpolygon = new QAction(MainWindow);
        actionpolygon->setObjectName(QString::fromUtf8("actionpolygon"));
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/new/prefix1/resouses/5.webp"), QSize(), QIcon::Normal, QIcon::Off);
        actionpolygon->setIcon(icon5);
        actionundo = new QAction(MainWindow);
        actionundo->setObjectName(QString::fromUtf8("actionundo"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/new/prefix1/resouses/6.webp"), QSize(), QIcon::Normal, QIcon::Off);
        actionundo->setIcon(icon6);
        actioncolor = new QAction(MainWindow);
        actioncolor->setObjectName(QString::fromUtf8("actioncolor"));
        actioncolor->setIcon(icon);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        MainWindow->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);

        toolBar->addAction(actionrectangle);
        toolBar->addAction(actionelliptical);
        toolBar->addAction(actionline);
        toolBar->addAction(actionfree_line);
        toolBar->addAction(actionpolygon);
        toolBar->addAction(actionundo);
        toolBar->addAction(actioncolor);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MiniDraw", nullptr));
        actionrectangle->setText(QApplication::translate("MainWindow", "rectangle", nullptr));
        actionelliptical->setText(QApplication::translate("MainWindow", "elliptical", nullptr));
        actionline->setText(QApplication::translate("MainWindow", "line", nullptr));
#ifndef QT_NO_TOOLTIP
        actionline->setToolTip(QApplication::translate("MainWindow", "lineqew", nullptr));
#endif // QT_NO_TOOLTIP
        actionfree_line->setText(QApplication::translate("MainWindow", "freehand", nullptr));
        actionpolygon->setText(QApplication::translate("MainWindow", "polygon", nullptr));
        actionundo->setText(QApplication::translate("MainWindow", "undo", nullptr));
        actioncolor->setText(QApplication::translate("MainWindow", "color and width", nullptr));
        toolBar->setWindowTitle(QApplication::translate("MainWindow", "MiniDraw", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
