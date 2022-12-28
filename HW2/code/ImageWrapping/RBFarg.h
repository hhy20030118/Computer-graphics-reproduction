#pragma once
#include "ImageOperator.h"
#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
class RBFarg :
    public ImageOperator
{
    void changepos() {

    }
    private:
        QPointF x, y;
        QVector<QPointF> alph;
        MatrixXd m;
        VectorXd l , f;

        double inline d(QPointF a, QPointF b) {//�������
            return sqrt((((double)a.x() - (double)b.x()) * ((double)a.x() - (double)b.x()) + ((double)a.y() - (double)b.y()) * ((double)a.y() - (double)b.y())));
        }

        double inline getr(int i) {
            double dis = 10000.0;
            for (int j = 0; j < ChangedPos.size(); j++) {
                if (d(ChangedPos.at(j).bpos, ChangedPos.at(i).bpos) < dis)dis = d(ChangedPos.at(j).bpos, ChangedPos.at(i).bpos);
            }
            return dis;
        }
        double inline limitrange(double r) {//ȡ������
            return r > 0 ? r : 0;
        }
        double inline getf(int i, QPointF p) {//����Ȩ��
            double dis = d((QPointF)ChangedPos.at(i).bpos, p);
            double disr = getr(i);
            //double range = limitrange(1.0 - pow((double)dis / 900.0, 2));
            double P = sqrt(dis * dis + disr * disr);
            return P;//������ƺ���֮��
        }
        QPointF getF(int i, QPointF p) {
            return ChangedPos.at(i).epos + p - ChangedPos.at(i).bpos;
        }
    public:
        void initial() {//��ʼ��
            m = Matrix<double, Dynamic, Dynamic>();
            m.resize((double)2 * (double)ChangedPos.size(), (double)2.0 * (double)ChangedPos.size());
            for (int i = 0; i <  ChangedPos.size(); i++) {
                for (int j = 0; j <  ChangedPos.size(); j++) {//��Ϊ�Ƕ�ά������Ҫ�����Ƿ�2�����������ͬ
                    m(2 * i, (2 * j)) = getf(j, (QPointF)ChangedPos.at(i).bpos);
                    m((2 * i+1), ((2 * j) + 1)) = getf(j, (QPointF)ChangedPos.at(i).bpos);
                    m((2 * i), (2 * j + 1)) = 0.0;
                    m((2 * i + 1), 2 * j) = 0.0;
                }
            }
            l = Vector<double, Dynamic>();//AX=B��B
            l.resize(2 * ChangedPos.size());//
            f = Vector<double, Dynamic>();//AX=B��X
            f.resize(2 * ChangedPos.size());
            for (int i = 0; i < ChangedPos.size(); i++) {
                l((2 * i)) = (double)(ChangedPos.at(i).epos - ChangedPos.at(i).bpos).x();
                l((2 * i+1)) = (double)(ChangedPos.at(i).epos - ChangedPos.at(i).bpos).y();
            }
            f = m.householderQr().solve(l);//��ⷽ����
            for (int i = 0; i < ChangedPos.size(); i++) {
                QPointF aa;
                aa.setX(f(2*i));
                aa.setY(f(2 * i+1));
                alph.append(aa);//ת���ɵ�
            }
        }
        QPointF inline getpos(int i, int j) {
            y = QPointF(i, j);
            x = QPointF(i, j);
            
            for (register int a = 0; a < ChangedPos.size(); ++a) {
                x += alph.at(a) * getf(a, y);
            }
            return x;
        }


};

