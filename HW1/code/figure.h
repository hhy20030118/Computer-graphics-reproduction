#ifndef FIGURE_H
#define FIGURE_H

#include <QMainWindow>
#include <vector>
#include <QPainter>
#include <QMouseEvent>
#include <math.h>
using namespace std;

class CFigure{       //关于图形的类 ， 主要存储对应图像的相关数据
protected :

    bool isfilled;  //关于填充颜色 ， 目前还没实现
    int type ;      //图像的类型
    bool enddrawing;
    QPen pen;   //笔刷
public:
    QVector<QPoint>path;

    CFigure(){

    }
    ~CFigure(){

    }
    int gettype(){
        return this->type;
    }
    void settype(int t){
        this->type = t;
    }
    bool get_enddrawing(){
        return this->enddrawing;
    }
    void set_enddrawing(bool b){
        this->enddrawing = b;
    }
    QPen getpen(){
        return pen;
    }
    void setpen(QPen p){
        pen = p;
    }
    virtual void draw(QPoint pos) {};//虚函数 ， 每一个类再实现

};
class CLine:public CFigure{
    void draw(QPoint pos){
        while(path.count()<2){
            QPoint p;
            path.append(p);
        }
        path.replace(1,pos); //只需要更新线段尾
    }
};
class CEllipse:public CFigure{
    void draw(QPoint pos){
        int h=200;    //采取的密度 ， 可以调整
        //QVector<QPoint> &lastLine = path;
        while(path.size()<2*h){
            QPoint p;
            path.append(p);//防止越界
        }
        double a=(pos.x()-path.at(0).x())/2;    //长轴
        double b=(pos.y()-path.at(0).y())/2;    //短轴
        path.replace(1,QPoint(path.at(0).x(),(pos.y()+path.at(0).y())/2)); //起点
        double d;

        for(long i=2 ; i<h ; ++i){
            d=(abs(2*a*i/h-a)/a)*(abs(2*a*i/h-a)/a); //中间数 ， 方便计算
            path.replace(i,QPoint(path.at(0).x()+2*a*i/h,path.at(1).y()+sqrt(1-d)*b)); //对称的设置点
            path.replace(2*h-i,QPoint(path.at(0).x()+2*a*i/h,path.at(1).y()-sqrt(1-d)*b));
            //lastLine.replace(i,QPoint(lastLine.at(0).x()+4*i,(ev->pos().y()+lastLine.at(0).y())/2+sqrt((1-abs(4*i-abs(ev->pos().x()-lastLine.at(0).x())/2)*abs(4*i-abs(ev->pos().x()-lastLine.at(0).x())/2)/((abs(ev->pos().x()-lastLine.at(0).x())/2)*(abs(ev->pos().x()-lastLine.at(0).x())/2)))*(abs(ev->pos().y()-lastLine.at(0).y())/2)*(abs(ev->pos().y()-lastLine.at(0).y())/2))));
        }
        path.replace(h,QPoint(pos.x(),path.at(1).y()));//把边界点补上
        path.replace(2*h-1,QPoint(path.at(0).x(),(pos.y()+path.at(0).y())/2));
    }

};
class CRect:public CFigure{
    void draw(QPoint pos){
        //QVector<QPoint> &lastLine = path;
        while(path.count()<5){
             QPoint p;
             path.append(p);//防止越界
        }
        QPoint p(path.at(0).x(),pos.y());
        path.replace(1,p);//对应长方体的四个点
        QPoint q(pos.x(),path.at(0).y());
        path.replace(3,q);
        path.replace(2,pos);
        path.replace(4,path.at(0));
    }
};
class CPolypon:public CFigure{
protected:
    bool isoperating = false;
    bool onending = false;
public:
    void draw(QPoint pos){
        if(!onending)    //保证不在结束的时候乱动造成干扰
        path.replace(path.size()-1,pos);
    }
    void end(){
        isoperating = false;
        onending = true;  //设置结束模式
        path.append(path.at(0));  //添加回去第一个点
    }
    void setoprate(bool b){
        isoperating = b;
    }
    bool getoprate(){
        return isoperating;
    }
};
class CFreehand:public CFigure{
    void draw(QPoint pos){
        path.append(pos);  //直接加点
    }
};
#endif // FIGURE_H
