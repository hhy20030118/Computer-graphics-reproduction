#include <Engine/MeshEdit/MinSurf_normal.h>
#include <Engine/Primitive/TriMesh.h>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
using namespace Ubpa;
using namespace Eigen;
using namespace std;
MinSurf_normal::MinSurf_normal(Ptr<TriMesh> triMesh)
	: MinSurf(triMesh)
{
	MinSurf::Init(triMesh);
	method = triMesh->method;
	type = triMesh->boundary_type;
}
inline double MinSurf_normal::distance(int a, int b) {
	return sqrt((heMesh->Vertices()[a]->pos.at(0) - heMesh->Vertices()[b]->pos.at(0)) * (heMesh->Vertices()[a]->pos.at(0) - heMesh->Vertices()[b]->pos.at(0)) + (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1)) * (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1)) + (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)) * (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)));
}
void MinSurf_normal::Minimize() {
	//初始化
	size_t nV = heMesh->NumVertices();  
	SparseMatrix<double> a;
	SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	VectorXi sizes;
	VectorXd x1 = Vector<double, Dynamic>(), x2 = Vector<double, Dynamic>(), x3 = Vector<double, Dynamic>(), b1 = Vector<double, Dynamic>(), b2 = Vector<double, Dynamic>(), b3 = Vector<double, Dynamic>();
	a.resize(nV, nV);
	x1.resize(nV);
	x2.resize(nV);
	x3.resize(nV);
	b1.resize(nV);
	b2.resize(nV);
	b3.resize(nV);
	sizes.resize(nV);
	for (register long i = 0; i < nV; ++i) {
		sizes(i) = 100;
	}
	a.reserve(sizes);
	vector<long> bdy;
	long current = -1;
	bool color[100000];//注意点的数量不能超
	for (register int i = 0; i < 100000; ++i) {
		color[i] = 0;
	}
	for (auto v : heMesh->Vertices()) {  //找到第一个边界
		if (v->IsBoundary()) {
			bdy.push_back(heMesh->Index(v));
			color[heMesh->Index(v)] = 1;
			for (auto u : v->AdjVertices()) {
				if (u->IsBoundary()) {
					bdy.push_back(heMesh->Index(u));
					current = heMesh->Index(u);
					color[heMesh->Index(u)] = 1;
				}
			}
			break;
		}
	}
	if (current == -1) {
		printf("ERROR::MinSurf::Run\n"
			"\t""has only one boundary\n");
		return;
	}
	bool finish = 0;
	while (!finish) {
		auto v = heMesh->Vertices()[current];
		finish = 1;
		for (auto u : v->AdjVertices()) {
			if (u->IsBoundary() && !color[heMesh->Index(u)]) {
				bdy.push_back(heMesh->Index(u));
				current = heMesh->Index(u);
				color[heMesh->Index(u)] = 1;
				finish = 0;
			}
		}
	}
	if (type == 1 && (bdy.size() + 1) >= 4) { //边界是正方形
		long q1 = bdy.size() / 4;
		long q2 = bdy.size() / 2;
		long q3 = bdy.size() - q1;
		auto v = heMesh->Vertices()[0];//遍历一遍边界
		v->pos.at(0) = 0.0;
		v->pos.at(1) = 0.0;
		v->pos.at(2) = 0.0;
		for (register int i = 1; i < q1; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = (1.0 / (double)q1) * i;
			v->pos.at(1) = 0.0;
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q1)];//遍历一遍边界
		v->pos.at(0) = 1.0;
		v->pos.at(1) = 0.0;
		v->pos.at(2) = 0.0;
		for (register int i = q1 + 1; i < q2; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = 1.0;
			v->pos.at(1) = (1.0 / (double)q1) * (i - q1);
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q2)];//遍历一遍边界
		v->pos.at(0) = 1.0;
		v->pos.at(1) = 1.0;
		v->pos.at(2) = 0.0;
		for (register int i = q2 + 1; i < q3; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = 1.0 - (1.0 / (double)q1) * (i - q2);
			v->pos.at(1) = 1.0;
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q3)];//遍历一遍边界
		v->pos.at(0) = 0.0;
		v->pos.at(1) = 1.0;
		v->pos.at(2) = 0.0;
		for (register int i = q3 + 1; i < bdy.size(); ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = 0.0;
			v->pos.at(1) = 1.0 - (1.0 / (double)q1) * (i - q3);
			v->pos.at(2) = 0.0;
		}

	}
	else if ((bdy.size() + 1) < 4) {
		printf("ERROR::MinSurf::Run\n""\t""point not enough\n");
		return;
	}

	if (type == 2 && (bdy.size() + 1) >= 4) { //边界是圆形
		auto v = heMesh->Vertices()[0];//遍历一遍边界
		v->pos.at(0) = 0.5;
		v->pos.at(1) = 1.0;
		v->pos.at(2) = 0.0;
		for (register int i = 1; i < bdy.size(); ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			double theta = 2 * 3.1415926535 / (double)(bdy.size()) * i;
			v->pos.at(0) = 0.5 + 0.5 * sin(theta);
			v->pos.at(1) = 0.5 + 0.5 * cos(theta);
			v->pos.at(2) = 0.0;
		}
	}
	else if ((bdy.size() + 1) < 4) {
		printf("ERROR::MinSurf::Run\n""\t""point not enough\n");
	}
	vector<vector<double>> w;
	vector<double> w1;
	for (register int i = 0; i < nV; ++i) {
		w1.push_back(0.0);
	}
	for (register int i = 0; i < nV; ++i) {
		w.push_back(w1);
	}
	if (method == 1) {
		for (auto v : heMesh->Vertices()) {
			int k = 0;
			for (auto u : v->AdjVertices()) {
				k++;
			}
			for (auto u : v->AdjVertices()) {
				w.at(heMesh->Index(v)).at(heMesh->Index(u)) = 1.0 / (double)k;
			}

		}
	}
	else {

		for (auto v : heMesh->Vertices()) {
			double sum = 0;
			for (auto u : v->AdjVertices()) {
				vector<int> target;
				for (auto v1 : v->AdjVertices()) {
					for (auto v2 : v1->AdjVertices()) {
						if (heMesh->Index(v2) == heMesh->Index(u)) {
							target.push_back(heMesh->Index(v1));
						}
					}
				}
				if (target.size() != 2) {
					//printf("ERROR::MinSurf::Run\n""\t""fault pointnum\n");
					continue;
				}
				else {
					double a = distance(heMesh->Index(v), heMesh->Index(u));
					double b = distance(heMesh->Index(v), target.at(0));
					double c = distance(target.at(0), heMesh->Index(u));
					double cos1 = (b * b + c * c - a * a) / (2 * b * c);
					double cot1 = cos1 / sqrt(1.0 - cos1 * cos1);
					if (a > b && a > c)cot1 = -cot1;
					b = distance(heMesh->Index(v), target.at(1));
					c = distance(target.at(1), heMesh->Index(u));
					double cos2 = (b * b + c * c - a * a) / (2 * b * c);

					double cot2 = cos2 / sqrt(1.0 - cos2 * cos2);
					if (a > b && a > c)cot2 = -cot2;
					w.at(heMesh->Index(v)).at(heMesh->Index(u)) = cot1 + cot2;
				}
			}
		}

	}
	for (auto v : heMesh->Vertices()) {
		if (!v->IsBoundary()) {
			int num = 0;
			double sum = 0;
			for (auto u : v->AdjVertices()) {
				if (!u->IsBoundary()) {
					if (method == 1)a.insert(heMesh->Index(v), heMesh->Index(u)) = -1.0;
					else {

						a.insert(heMesh->Index(v), heMesh->Index(u)) = -w.at(heMesh->Index(v)).at(heMesh->Index(u));

					}
				}
				num++;
				if (w.at(heMesh->Index(v)).at(heMesh->Index(u)) == 0) {
					printf("ERROR::MinSurf::Run\n""\t""point not enough\n");
				}
				sum += w.at(heMesh->Index(v)).at(heMesh->Index(u));
			}

			if (method == 1)a.insert(heMesh->Index(v), heMesh->Index(v)) = num;
			else { a.insert(heMesh->Index(v), heMesh->Index(v)) = sum; }
			b1(heMesh->Index(v)) = 0;
			b2(heMesh->Index(v)) = 0;
			b3(heMesh->Index(v)) = 0;
			for (auto u : v->AdjVertices()) {
				if (u->IsBoundary()) {
					//if (type == 0) {
					if (method == 1) {   //普通
						b1(heMesh->Index(v)) += u->pos.at(0);
						b2(heMesh->Index(v)) += u->pos.at(1);
						b3(heMesh->Index(v)) += u->pos.at(2);
					}
					else {  //cot权重
						b1(heMesh->Index(v)) += w.at(heMesh->Index(v)).at(heMesh->Index(u)) * (u->pos.at(0));
						b2(heMesh->Index(v)) += w.at(heMesh->Index(v)).at(heMesh->Index(u)) * (u->pos.at(1));
						b3(heMesh->Index(v)) += w.at(heMesh->Index(v)).at(heMesh->Index(u)) * (u->pos.at(2));
					}

					//}

				}
			}
		}
		else {
			a.insert(heMesh->Index(v), heMesh->Index(v)) = 1.0;
			b1(heMesh->Index(v)) = v->pos.at(0);
			b2(heMesh->Index(v)) = v->pos.at(1);
			b3(heMesh->Index(v)) = v->pos.at(2);
		}

	}
	a.makeCompressed();
	solver.compute(a);
	if (solver.info() != Eigen::Success)
	{
		throw std::exception("Compute Matrix is error");
		return;
	}
	x1 = solver.solve(b1);
	x2 = solver.solve(b2);
	x3 = solver.solve(b3);
	for (auto v : heMesh->Vertices()) {
		v->pos.at(0) = x1(heMesh->Index(v));
		v->pos.at(1) = x2(heMesh->Index(v));
		v->pos.at(2) = x3(heMesh->Index(v));
	}

}
bool MinSurf_normal::Run() {
	//std::cout<<"test"<<std::endl;
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}
	Minimize();

	// half-edge structure -> triangle mesh
	size_t nV = heMesh->NumVertices();
	size_t nF = heMesh->NumPolygons();
	vector<pointf3> positions;
	vector<unsigned> indice;
	positions.reserve(nV);
	indice.reserve(3 * nF);
	for (auto v : heMesh->Vertices())
		positions.push_back(v->pos.cast_to<pointf3>());
	for (auto f : heMesh->Polygons()) { // f is triangle
		for (auto v : f->BoundaryVertice()) // vertices of the triangle
			indice.push_back(static_cast<unsigned>(heMesh->Index(v)));
	}

	triMesh->Init(indice, positions);

	return true;
}

