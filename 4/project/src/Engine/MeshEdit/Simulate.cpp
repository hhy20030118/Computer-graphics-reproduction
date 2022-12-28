#include <Engine/MeshEdit/Simulate.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace Ubpa;

using namespace std;
using namespace Eigen;


void Simulate::Clear() {
	this->positions.clear();
	this->velocity.clear();
}
inline void Simulate::Matrix_value(Matrix3d s, MatrixXd* t, int i, int j) {//给矩阵赋值
	(*t)(3 * i, 3 * j) += s(0, 0);
	(*t)(3 * i + 1, 3 * j) += s(1, 0);
	(*t)(3 * i + 2, 3 * j) += s(2, 0);
	(*t)(3 * i, 3 * j + 1) += s(0, 1);
	(*t)(3 * i + 1, 3 * j + 1) += s(1, 1);
	(*t)(3 * i + 2, 3 * j + 1) += s(2, 1);
	(*t)(3 * i, 3 * j + 2) += s(0, 2);
	(*t)(3 * i + 1, 3 * j + 2) += s(1, 2);
	(*t)(3 * i + 2, 3 * j + 2) += s(2, 2);
}
bool Simulate::Init() {
	//Clear();
	isinit = false;
	this->velocity .resize(positions.size());//速度初始化
	for (int i = 0; i < positions.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->velocity[i][j] = 0;
		}
	}
	this->edge_origin.resize(edgelist.size());
	for (int k = 0; k < edgelist.size()/2; k++)//原长初始化
	{
		int i = edgelist.at(2 * k);
		int j = edgelist.at(2 * k + 1);
		this->edge_origin[k] = sqrt((this->positions[j][0]- this->positions[i][0])*(this->positions[j][0] - this->positions[i][0])+(this->positions[j][1] - this->positions[i][1])*(this->positions[j][1] - this->positions[i][1])+(this->positions[j][2] - this->positions[i][2])*(this->positions[j][2] - this->positions[i][2]));
	}
	this->force.resize(positions.size());
	this->positions_next.resize(positions.size());
	for (int i = 0; i < positions.size(); i++)
	{
		this->positions_next[i] = this->positions[i];
	}
	xn.resize(3 * positions.size());
	vn.resize(3 * positions.size());
	for (int i = 0; i < positions.size(); i++)  //上一帧的位置
	{
		xn(3 * i) = this->positions[i][0];
		xn(3 * i+1) = this->positions[i][1];
		xn(3 * i+2) = this->positions[i][2];
		vn(3 * i) = this->velocity[i][0];
		vn(3 * i + 1) = this->velocity[i][1];
		vn(3 * i + 2) = this->velocity[i][2];
	}
	int l = 2;
	for (int i = 0; i < 400; ++i) { //自己写
		if (l % 21 == 1) {
			l++;
			i--;
			continue;
		}
		facelist.push_back(l-1);
		facelist.push_back(l-2);
		facelist.push_back(l +19);
		facelist.push_back(l-1);
		facelist.push_back(l+19);
		facelist.push_back(l + 20);
		l++;
	}
	return true;
}
void Simulate::Matrix_init() {
	//矩阵初始化
	size_t n = 3 * positions.size();
	size_t s = 3 * edge_origin.size();
	L.resize(n, n);
	J.resize(n, s);
	L = MatrixXd::Zero(n, n);
	J = MatrixXd::Zero(n, s);
	for (register int i = 0; i < n; ++i) {
		for (register int j = 0; j < n; ++j) {
			L(i, j) = 0;
		}
		for (register int j = 0; j < s; ++j) {
			J(i, j) = 0;
		}
	}
	for (register int c = 0; c < edgelist.size() / 2; ++c) {
		int i = edgelist.at(2 * c);
		int j = edgelist.at(2 * c + 1);
		double k = ks / edge_origin[c]*0.1;
		Matrix_value(k * Matrix3d::Identity(), &L, i, i);
		Matrix_value(-k * Matrix3d::Identity(), &L, j, i);
		Matrix_value(k * Matrix3d::Identity(), &L, j, j);
		Matrix_value(-k * Matrix3d::Identity(), &L, i, j);
		Matrix_value(-k * Matrix3d::Identity(), &J, i, c);
		Matrix_value(k * Matrix3d::Identity(), &J, j, c);
	}
	A = (mess * MatrixXd::Identity(n, n) + h * h * L).sparseView();
	B = (h * h * J).sparseView();
	A.makeCompressed();
	solverA.compute(A);
	if (solverA.info() != Eigen::Success)
	{
		throw std::exception("Compute Matrix is error");
		return ;
	}
	B.makeCompressed();
}
bool Simulate::Run() {
	//type 决定运行哪个
	if (type == 1) {
		SimulateOnce();
		divide = 30;
	}
	else if (type == 2) {
		SimulateOnce_implict();
		divide = 10;
	}
	else if (type == 4) {
		h = 0.001f;
		SimulateOnce_PBD();
	}
	else {	
		if (!isinit) {
			cout << "On Initing" << endl;
			Matrix_init();
			cout << "done" << endl;
			isinit = true;
		}
		SimulateOnce_implict_fast();	
	}
	return true;
}
void Ubpa::Simulate::SetLeftFix()
{
	//固定网格x坐标最小点
	fixed_id.clear();
	double x = 100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][0] < x)
		{
			x = positions[i][0];
		}
	}

	for (int i = 0; i < positions.size(); i++)
	{
		if (abs(positions[i][0] - x) < 1e-5)
		{
			fixed_id.push_back(i);
		}
	}

	Init();
}
void Ubpa::Simulate::SetpointFix()
{
	//固定网格x坐标最小点
	fixed_id.clear();
	double y = -100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][1] > y)
		{
			y = positions[i][1];
		}
	}
	double x1 = -100000, x2 = 100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][0] < x2)
		{
			x2 = positions[i][0];
		}
		if (positions[i][0] > x1)
		{
			x1 = positions[i][0];
		}
	}
	for (int i = 0; i < positions.size(); i++)
	{
		if (abs(positions[i][0] - x1) < 0.2 && abs(positions[i][1] - y) < 1e-5)
		{
			fixed_id.push_back(i);
		}
		if (abs(positions[i][0] - x2) < 0.2 && abs(positions[i][1] - y) < 1e-5)
		{
			fixed_id.push_back(i);
		}
	}

	Init();
}
void Ubpa::Simulate::SetUponFix()
{
	fixed_id.clear();
	double y = -100000;
	for (int i = 0; i < positions.size(); i++)
	{
		if (positions[i][1] > y)
		{
			y = positions[i][1];
		}
	}

	for (int i = 0; i < positions.size(); i++)
	{
		if (abs(positions[i][1] - y) < 1e-5)
		{
			fixed_id.push_back(i);
		}
	}

	Init();
}
void Ubpa::Simulate::SetMode() {
	if (mode == 1)mode = 2;
	else mode = 1;
}
void Simulate::SimulateOnce() {
	//半隐式（作了一些修改）
	for (register int r = 0; r < divide; r++) {
		for (int i = 0; i < positions.size(); i++)
		{	
			this->positions_next[i][0] += h /divide/2 * this->velocity[i][0];
			this->positions_next[i][1] += h / divide/2 * this->velocity[i][1];
			this->positions_next[i][2] += h / divide/2 * this->velocity[i][2];

		}
		for (int i = 0; i < positions.size(); i++)
		{
			pointf3 v;
			v[0] = this->velocity[i][0];
			v[1] = this->velocity[i][1];
			v[2] = this->velocity[i][2];
			auto l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
			//double f_1 = 0.1*(1 - exp(-sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])));
			//if (f_1 < 0)cout << "wrong!" << endl;
			double f_1 = 0.0;
			this->force[i][0] = -f_1 * v[0] ;
			this->force[i][1] = -f_1 * v[1] ;
			this->force[i][2] = -f_1 * v[2]  - 9.8 * mess;
			//this->force[i][0] = 0;
			//this->force[i][1] = 0;
			//this->force[i][2] = - 9.8 * mess;
		}

		//		<< "\t" << "not implemented" << endl;
		for (register int a = 0; a < edgelist.size() / 2; ++a) {
			int i = edgelist.at(2 * a);
			int j = edgelist.at(2 * a + 1);
			pointf3 d, v, f , d_last;
			d[0] = this->positions_next[j][0] - this->positions_next[i][0];
			d[1] = this->positions_next[j][1] - this->positions_next[i][1];
			d[2] = this->positions_next[j][2] - this->positions_next[i][2];
			d_last[0] = this->positions[j][0] - this->positions[i][0];
			d_last[1] = this->positions[j][1] - this->positions[i][1];
			d_last[2] = this->positions[j][2] - this->positions[i][2];
			v[0] = this->velocity[j][0]-this->velocity[i][0];
			v[1] = this->velocity[j][1] - this->velocity[i][1];
			v[2] = this->velocity[j][2] - this->velocity[i][2];
			auto l_0 = edge_origin[a];
			float l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			float l_last = sqrt(d_last[0] * d_last[0] + d_last[1] * d_last[1] + d_last[2] * d_last[2]);
			auto v_re = abs(v[0] * d[0] + v[1] * d[1] + v[2] * d[2]) / l;
			auto v_pe = abs(v[0] * d[0] + v[1] * d[1] + v[2] * d[2]);
			double f_0 = ks * (l / l_0 - 1);
			/*if ((l / l_0 < 1.25 && l / l_0 > 0.5)) {
				auto x = (l_0 - l) / l_0;
				auto y = x * x * (-16 * x + 12);
				f_0 = f_0 * y;
			}
			if (l / l_0 < 1)f_0 = f_0 * 0.01;
			if (l / l_0 < 0.75)f_0 = 0.0;
			if (l / l_0 > 1.25)f_0 = f_0*exp((l/l_0-1.25));*/  //this is some other method

			//if (l / l_0 < 1)f_0 = f_0 * 0.01;

			//f_0 *= (1 - exp(-abs(l/l_0-1))*exp(-0.1*abs(l_last-l)));

			//double f_0 = ks * (l / l0 - 1.0);
			//double f_0 = ks * (l / l0 - 1.0) + kd * (v[0] * d[0] + v[1] * d[1] + v[2] * d[2]) / (l0 * l);
			f[0] = f_0 * d[0] / l;
			f[1] = f_0 * d[1] / l;
			f[2] = f_0 * d[2] / l;
			force[i][0] += f[0];
			force[i][1] += f[1];
			force[i][2] += f[2];
			force[j][0] -= f[0];
			force[j][1] -= f[1];
			force[j][2] -= f[2];   //换成vector操作更简洁一些


		}
		int p = 0;
		for (int j = 0; j < fixed_id.size(); ++j) {
			if (mode == 2)fixed_id.pop_back();
		}
		for (int i = 0; i < positions.size(); i++)
		{

			if (fixed_id.size()>0&&fixed_id.at(p) == i) {
				p++;
				if (p > fixed_id.size() - 1)p = fixed_id.size() - 1;
				continue;
			}
			this->velocity[i][0] += h / divide * force[i][0] / mess;
			this->velocity[i][1] += h / divide * force[i][1] / mess;
			this->velocity[i][2] += h / divide * force[i][2] / mess;
			this->positions[i][0] += h / divide * this->velocity[i][0];
			this->positions[i][1] += h / divide * this->velocity[i][1];
			this->positions[i][2] += h / divide * this->velocity[i][2];
			/*if (this->positions[i][2] < -0.3) {
				this->positions[i][2] = -0.3;
				this->velocity[i][2] = 0.0;
			}*/
			this->positions_next[i][0] = this->positions[i][0];
			this->positions_next[i][1] = this->positions[i][1];
			this->positions_next[i][2] = this->positions[i][2];

		}


	}
	

}
void Simulate::SimulateOnce_PBD() {
	int p = 0;
	for (int i = 0; i < positions.size(); i++)
	{
		if (fixed_id.size() > 0 && fixed_id.at(p) == i) {
			p++;
			if (p > fixed_id.size() - 1)p = fixed_id.size() - 1;
			continue;
		}
		this->positions_next[i][0] += h  * this->velocity[i][0];
		this->positions_next[i][1] += h  * this->velocity[i][1];
		this->positions_next[i][2] += h  * this->velocity[i][2];

	}
	for (register int u = 0; u < 5; ++u) {
		for (register int a = 0; a < facelist.size() / 3; ++a) {
			int i = facelist.at(3 * a);
			int j = facelist.at(3 * a + 1);
			int k = facelist.at(3 * a + 2);
			int l = 0;
			if (j == i - 1)l = k + 1;
			else l = i - 1;
			Vector3d s , d1 , d2 , d3;
			s[0] = (this->positions_next[i][0] + this->positions_next[j][0] + this->positions_next[k][0]) / 3;
			s[1] = (this->positions_next[i][1] + this->positions_next[j][1] + this->positions_next[k][1]) / 3;
			s[2] = (this->positions_next[i][2] + this->positions_next[j][2] + this->positions_next[k][2]) / 3;
			/*d1[0] = s[0] - this->positions_next[i][0];
			d1[1] = s[1] - this->positions_next[i][1];
			d1[2] = s[2] - this->positions_next[i][2];
			d2[0] = s[0] - this->positions_next[j][0];
			d2[1] = s[1] - this->positions_next[j][1];
			d2[2] = s[2] - this->positions_next[j][2];
			d3[0] = s[0] - this->positions_next[k][0];
			d3[1] = s[1] - this->positions_next[k][1];
			d3[2] = s[2] - this->positions_next[k][2];
			double r1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
			double r2 = sqrt(d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2]);
			double r3 = sqrt(d3[0] * d3[0] + d3[1] * d3[1] + d3[2] * d3[2]);*/
			if (j == i - 1) {
				d2[0] = s[0] - this->positions_next[l][0];
				d2[1] = s[1] - this->positions_next[l][1];
				d2[2] = s[2] - this->positions_next[l][2];
				double r2 = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
				double kk = -((s[0] - this->positions_next[i][0]) * (this->positions_next[l][0] - s[0]) + (s[1] - this->positions_next[i][1]) * (this->positions_next[l][1] - s[1]) + (s[2] - this->positions_next[i][2]) * (this->positions_next[l][2] - s[2])) / (r2);
				r2 = sqrt(r2);
				d1[0] = kk * (-d2[0]) + s[0] - this->positions_next[i][0];
				d1[1] = kk * (-d2[1]) + s[1] - this->positions_next[i][1];
				d1[2] = kk * (-d2[2]) + s[2] - this->positions_next[i][2];
				double r1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
				d3 = d1;
				d1 = d1 / r1 * 0.7071 * 0.05 + d2 / r2 * 0.2357 * 0.05;
				d3 = d3 / r1 * 0.7071 * 0.05 - d2 / r2 * 0.2357 * 0.05;
				d2 = d2 / r2 * 0.471405 * 0.05;
				d1 = -d1;
				/*double x1 = r1 - 0.745356 * 0.05;
				double x2 = r2 - 0.471405*0.05;
				double x3 = r3 - 0.745356*0.05;
				d1 = d1 / r1 * x1;
				d2 = d2 / r2 * x2;
				d3 = d3 / r3 * x3;*/
			}
			else {
				d3[0] = s[0] - this->positions_next[l][0];
				d3[1] = s[1] - this->positions_next[l][1];
				d3[2] = s[2] - this->positions_next[l][2];
				double r3 = d3[0] * d3[0] + d3[1] * d3[1] + d3[2] * d3[2];
				double kk = -((s[0] - this->positions_next[i][0]) * (this->positions_next[l][0] - s[0]) + (s[1] - this->positions_next[i][1]) * (this->positions_next[l][1] - s[1]) + (s[2] - this->positions_next[i][2]) * (this->positions_next[l][2] - s[2])) / (r3);
				r3 = sqrt(r3);
				d1[0] = kk * (-d3[0]) + s[0] - this->positions_next[i][0];
				d1[1] = kk * (-d3[1]) + s[1] - this->positions_next[i][1];
				d1[2] = kk * (-d3[2]) + s[2] - this->positions_next[i][2];
				double r1 = sqrt(d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2]);
				d2 = d1;
				d1 = d1 / r1 * 0.7071 * 0.05 + d3 / r3 * 0.2357 * 0.05;
				d2 = d2 / r1 * 0.7071 * 0.05 - d3 / r3 * 0.2357 * 0.05;
				d3 = d3 / r3 * 0.471405 * 0.05;
				d1 = -d1;
			}
			this->positions_next[i][0] = s[0] + d1[0];
			this->positions_next[i][1] = s[1] + d1[1];
			this->positions_next[i][2] = s[2] + d1[2];
			this->positions_next[j][0] = s[0] + d2[0];
			this->positions_next[j][1] = s[1] + d2[1];
			this->positions_next[j][2] = s[2] + d2[2];
			this->positions_next[k][0] = s[0] + d3[0];
			this->positions_next[k][1] = s[1] + d3[1];
			this->positions_next[k][2] = s[2] + d3[2];
		}
		/*for (register int a = 0; a < edgelist.size() / 2; ++a) { //基于边的
			int i = edgelist.at(2 * a);
			int j = edgelist.at(2 * a + 1);
			Vector3d d, c1, c2;
			d[0] = this->positions_next[j][0] - this->positions_next[i][0];
			d[1] = this->positions_next[j][1] - this->positions_next[i][1];
			d[2] = this->positions_next[j][2] - this->positions_next[i][2];
			double l = edge_origin[a];
			double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
			double x = r - l;
			d = d / r * x*0.5;
			this->positions_next[i][0] += d[0];
			this->positions_next[i][1] += d[1];
			this->positions_next[i][2] += d[2];
			this->positions_next[j][0] -= d[0];
			this->positions_next[j][1] -= d[1];
			this->positions_next[j][2] -= d[2];
		}*/
	}
	p = 0;
	for (int i = 0; i < positions.size(); i++)
	{
		if (fixed_id.size() > 0 && fixed_id.at(p) == i) {
			p++;
			if (p > fixed_id.size() - 1)p = fixed_id.size() - 1;
			continue;
		}
		this->velocity[i][0] = (this->positions_next[i][0] - this->positions[i][0]) / h;
		this->velocity[i][1] = (this->positions_next[i][1] - this->positions[i][1]) / h;
		this->velocity[i][2] = (this->positions_next[i][2] - this->positions[i][2]) / h;
		this->velocity[i][2] -= h * 0.1;
		this->positions[i][0] = this->positions_next[i][0];
		this->positions[i][1] = this->positions_next[i][1];
		this->positions[i][2] = this->positions_next[i][2];
	}
	for (register int i = 0; i < fixed_id.size(); ++i) {
		this->positions_next[fixed_id[i]][0] = this->positions[fixed_id[i]][0];
		this->positions_next[fixed_id[i]][1] = this->positions[fixed_id[i]][1];
		this->positions_next[fixed_id[i]][2] = this->positions[fixed_id[i]][2];
	}
}
void Simulate::SimulateOnce_implict() {
	VectorXd x = Vector<double, Dynamic>(), b = Vector<double, Dynamic>();
	VectorXd y = Vector<double, Dynamic>();
	VectorXd f_g = Vector<double, Dynamic>();
	VectorXd v = Vector<double, Dynamic>();
	MatrixXd a;

	size_t n = 3 * positions.size();
	a.resize(n, n);
	x.resize(n);
	y.resize(n);
	f_g.resize(n);
	b.resize(n);
	v.resize(n);
	for (register int r = 0; r < divide; r++) {
		a = MatrixXd::Zero(n, n);
		//外力矩阵
		for (register size_t i = 0; i < positions.size(); ++i) {
			f_g(3 * i) = 0.0;
			f_g(3 * i + 1) = 0.0;
			f_g(3 * i + 2) = -1.0 * mess;
		}
		//y
		y = xn + h /divide * vn + h / divide * h / divide / mess * f_g;
		x = y;
		for (register int k = 0; k < 2; ++k) {
			b = mess * (x - y);

			//求解梯度矩阵
			for (register int c = 0; c < edgelist.size() / 2; ++c) {
				int i = edgelist.at(2 * c);
				int j = edgelist.at(2 * c + 1);
				pointf3 v, f;
				Matrix3d gradient;
				Vector3d d;
				d(0) = x(3 * j) - x(3 * i);
				d(1) = x(3 * j + 1) - x(3 * i + 1);
				d(2) = x(3 * j + 2) - x(3 * i + 2);
				double r = sqrt(d(0) * d(0) + d(1) * d(1) + d(2) * d(2)), l = edge_origin[c];
				gradient = ks * (1 / r - 1 / l) * Matrix3d::Identity() - ks / r / r / r * (d * d.transpose());
				double f1 = h / divide * h / divide * ks * (r / l - 1) / r;
				if (r < l) {
					gradient *= 0.01;
					f1 *= 0.01;
				}
				else {
					//gradient *= 10.0;
					//f1 *= 10.0;
				}
				b(3 * i) -= f1 * d(0);
				b(3 * i + 1) -= f1 * d(1);
				b(3 * i + 2) -= f1 * d(2);
				b(3 * j) -= -f1 * d(0);
				b(3 * j + 1) -= -f1 * d(1);
				b(3 * j + 2) -= -f1 * d(2);
				Matrix_value(gradient, &a, i, i);  //赋值
				Matrix_value(gradient, &a, j, j);
				Matrix_value(-gradient, &a, i, j);
				Matrix_value(-gradient, &a, j, i);
			}
			a = mess * MatrixXd::Identity(n, n) - h / divide * h / divide * a;
			x = x - a.reverse() * b;
			for (register int i = 0; i < fixed_id.size(); ++i) {//处理固定点
				x(3 * fixed_id[i]) = xn(3 * fixed_id[i]);
				x(3 * fixed_id[i] + 1) = xn(3 * fixed_id[i] + 1);
				x(3 * fixed_id[i] + 2) = xn(3 * fixed_id[i] + 2);
			}
		}
		vn = (x - xn) / (h / divide);
		xn = x;
		for (int i = 0; i < positions.size(); i++)
		{
			this->positions[i][0] = xn(3 * i);
			this->positions[i][1] = xn(3 * i + 1);
			this->positions[i][2] = xn(3 * i + 2);
			this->velocity[i][0] = vn(3 * i);
			this->velocity[i][1] = vn(3 * i + 1);
			this->velocity[i][2] = vn(3 * i + 2);
		}
	}
}
void Simulate::SimulateOnce_implict_fast() {
	VectorXd c = Vector<double, Dynamic>(), e = Vector<double, Dynamic>();
	VectorXd x = Vector<double, Dynamic>(), d = Vector<double, Dynamic>();
	VectorXd y = Vector<double, Dynamic>();
	VectorXd f_g = Vector<double, Dynamic>();

	size_t n = 3 * positions.size();
	size_t s = 3 * edge_origin.size();
	vector<double> e1;
	//i.resize(n, n);
	x.resize(n);
	y.resize(n);
	f_g.resize(n);
	d.resize(s);
	e.resize(n);

	//x_next.resize(n);
	vector<double> y_debug;
	vector<double> x_debug;

	//测试
	//矩阵初始化
	/*L = MatrixXd::Zero(n, n);
	J = MatrixXd::Zero(n, s);
	for (register int i = 0; i < n; ++i) {
		for (register int j = 0; j < n; ++j) {
			L(i, j) = 0.0;
		}
		for (register int j = 0; j < s; ++j) {
			J(i, j) = 0.0;
		}
	}
	for (register int c = 0; c < edgelist.size() / 2; ++c) {
		int i = edgelist.at(2 * c);
		int j = edgelist.at(2 * c + 1);
		double k = ks / edge_origin[c];
		Vector3d d1;
		d1(0) = xn(3 * j) - xn(3 * i);
		d1(1) = xn(3 * j + 1) - xn(3 * i + 1);
		d1(2) = xn(3 * j + 2) - xn(3 * i + 2);
		double r = sqrt(d1(0) * d1(0) + d1(1) * d1(1) + d1(2) * d1(2)), l = edge_origin[c];
		if ((r / l < 1.25 && r / l > 0.5)) {
			auto x = (l-r) / l;
			auto y = x * x * (-16 * x + 12);
			k = k * y;
		}
		if (r/l < 1)k = k * 0.01;
		//if (r / l < 0.75)k = 0.0;
		if (r / l > 1.25)k = k * exp((r / l - 1.25));
		//if (r < l)k *= 0.1;
		//else k *= 1.0*exp(r/l-1);
		Matrix_value(k * Matrix3d::Identity(), &L, i, i);
		Matrix_value(-k * Matrix3d::Identity(), &L, j, i);
		Matrix_value(k * Matrix3d::Identity(), &L, j, j);
		Matrix_value(-k * Matrix3d::Identity(), &L, i, j);
		Matrix_value(-k * Matrix3d::Identity(), &J, i, c);
		Matrix_value(k * Matrix3d::Identity(), &J, j, c);
	}
	A = (mess * MatrixXd::Identity(n, n) + h * h * L).sparseView();
	B = (h * h * J).sparseView();
	A.makeCompressed();
	solverA.compute(A);
	if (solverA.info() != Eigen::Success)
	{
		throw std::exception("Compute Matrix is error");
		return;
	}
	B.makeCompressed();*/

	//外力矩阵
	for (register size_t i = 0; i < positions.size(); ++i) {
		f_g(3 * i) = 0.0;
		f_g(3 * i + 1) = 0.0;
		f_g(3 * i + 2) = -9.8 * mess;
	}
	//y
	y = xn + h * vn + h * h / mess * f_g;
	y = mess * y;
	x = y;
	for (register int k = 0; k < times; ++k) {
		//迭代d
		for (register int c = 0; c < edgelist.size() / 2; ++c) {
			int i = edgelist.at(2 * c);
			int j = edgelist.at(2 * c + 1);
			Vector3d d1;
			d1(0) = x(3 * j) - x(3 * i);
			d1(1) = x(3 * j + 1) - x(3 * i + 1);
			d1(2) = x(3 * j + 2) - x(3 * i + 2);
			double r = sqrt(d1(0) * d1(0) + d1(1) * d1(1) + d1(2) * d1(2)), l = edge_origin[c];
			d1 = d1 * l / r;
			d(3 * c) = d1(0);
			d(3 * c + 1) = d1(1);
			d(3 * c + 2) = d1(2);
		}
		//迭代x76
		e = B * d +y;
		for (register int i = 0; i < n; ++i) {
			e1.push_back(e(i));
		}
		x = solverA.solve(e);
		for (register int i = 0; i < fixed_id.size(); ++i) {
			x(3 * fixed_id[i]) = xn(3 * fixed_id[i]);
			x(3 * fixed_id[i] + 1) = xn(3 * fixed_id[i] + 1);
			x(3 * fixed_id[i] + 2) = xn(3 * fixed_id[i] + 2);
		}
	}
	vn = (x - xn) / h;
	xn = x;
	for (int i = 0; i < positions.size(); i++)
	{
		this->positions[i][0] = xn(3 * i);
		this->positions[i][1] = xn(3 * i + 1);
		this->positions[i][2] = xn(3 * i + 2);
		this->velocity[i][0] = vn(3 * i);
		this->velocity[i][1] = vn(3 * i + 1);
		this->velocity[i][2] = vn(3 * i + 2);
	}



}
