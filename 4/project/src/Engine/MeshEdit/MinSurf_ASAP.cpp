#include <Engine/MeshEdit/MinSurf_ASAP.h>
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
MinSurf_ASAP::MinSurf_ASAP(Ptr<TriMesh> triMesh)
	: MinSurf(triMesh)
{
	MinSurf::Init(triMesh);
}
inline double MinSurf_ASAP::distance(int a, int b) {
	return sqrt((heMesh->Vertices()[a]->pos.at(0) - heMesh->Vertices()[b]->pos.at(0)) * (heMesh->Vertices()[a]->pos.at(0) - heMesh->Vertices()[b]->pos.at(0)) + (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1)) * (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1)) + (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)) * (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)));
}
inline double MinSurf_ASAP::acos(int a, int b, int c) {//���ء�A��cosֵ
	auto a1 = distance(a, b);
	auto b1 = distance(a, c);
	auto c1 = distance(b, c);
	return (a1 * a1 + b1 * b1 - c1 * c1) / (2 * a1 * b1);
}
void MinSurf_ASAP::ASAP() {
	size_t nV = heMesh->NumVertices();
	SparseMatrix<double> a;
	SimplicialLDLT<SparseMatrix<double>> solver;
	//HouseholderQR<SparseMatrix<double>> solver;;
	VectorXi sizes;
	VectorXd x1 = Vector<double, Dynamic>(), x2 = Vector<double, Dynamic>(), x3 = Vector<double, Dynamic>(), b1 = Vector<double, Dynamic>(), b2 = Vector<double, Dynamic>(), b3 = Vector<double, Dynamic>();
	size_t nT = heMesh->NumPolygons();  //������Էŵ�ͷ�ļ���
	size_t n = 2 * (nV + nT);
	//vector<vector<double>> x;  //�涥��Ⱦ��������ԭʼ����xij
	QVector<float> x_sub;//�ݴ�
	//vector<vector<double>> cot;//����������cotֵ
	QVector<float> cot_sub;//�ݴ�
	vector<vector<float>> a_pre;//�ݴ����
	vector<float> a_pre_sub;
	for (register size_t j = 0; j < 6; ++j) {
		x_sub.push_back(0.0);
	}
	for (register size_t i = 0; i < nT; ++i) {
		x.push_back(x_sub);
	}
	for (register size_t j = 0; j < 3; ++j) {
		cot_sub.push_back(0.0);
	}
	for (register size_t i = 0; i < nT; ++i) {
		cot.push_back(cot_sub);
	}

	a.resize(n, n);
	x1.resize(n);
	b1.resize(n);
	sizes.resize(n);
	for (register size_t i = 0; i < n; ++i) {
		sizes(i) = 50;
	}
	a.reserve(sizes);
	for (register size_t i = 0; i < n; ++i) {
		a_pre_sub.push_back(0.0);
	}
	for (register size_t j = 0; j < n; ++j) {
		a_pre.push_back(a_pre_sub);
	}
	for (register size_t j = 0; j < n; ++j) {
		b1(j) = 0.0;
	}
	for (register size_t j = 0; j < n; ++j) {
		x1(j) = 0.0;
	}
	//�����ǶԾ���ĳ�ʼ��
	size_t k = 0;
	for (auto m : heMesh->Polygons()) {  //��ʼ��x��cot
		auto v1 = m->BoundaryVertice()[0];
		auto v2 = m->BoundaryVertice()[1];
		auto v3 = m->BoundaryVertice()[2];
		auto cos = acos(heMesh->Index(v1), heMesh->Index(v2), heMesh->Index(v3));
		auto sin = sqrt(1.0 - cos * cos);
		auto cos2 = acos(heMesh->Index(v2), heMesh->Index(v1), heMesh->Index(v3));
		auto sin2 = sqrt(1.0 - cos2 * cos2);
		auto cos3 = acos(heMesh->Index(v3), heMesh->Index(v2), heMesh->Index(v1));
		auto sin3 = sqrt(1.0 - cos3 * cos3);
		x[k][2] = distance(heMesh->Index(v1), heMesh->Index(v2));
		x[k][4] = distance(heMesh->Index(v1), heMesh->Index(v3)) * cos;
		x[k][5] = distance(heMesh->Index(v1), heMesh->Index(v3)) * sin;
		cot[k][0] = cos / sin;
		cot[k][1] = cos2 / sin2;
		cot[k][2] = cos3 / sin3;
		k++;
	}

	auto fixed_point1 = heMesh->Index(heMesh->Boundaries()[0][0]->Origin());
	auto fixed_point2 = heMesh->Index(heMesh->Boundaries()[0][1]->Origin());
	for (auto v : heMesh->Vertices()) {
		if (v->IsBoundary() && distance(fixed_point1, heMesh->Index(v)) > distance(fixed_point1, fixed_point2)) {
			fixed_point2 = heMesh->Index(v);
		}
	}

	//���ê��
	a_pre.at(fixed_point1).at(fixed_point1) = 1.0;
	a_pre.at(fixed_point2).at(fixed_point2) = 1.0;
	a_pre.at(nV + fixed_point1).at(nV + fixed_point1) = 1.0;
	a_pre.at(nV + fixed_point2).at(nV + fixed_point2) = 1.0;
	b1(fixed_point1) = 0.5;
	b1(fixed_point2) = 0.5;
	b1(nV + fixed_point1) = 0.0;
	b1(nV + fixed_point2) = 1.0;
	//���ԸĽ�
	for (auto h : heMesh->HalfEdges()) {//�Ծ���ֵ
		auto v1 = h->Origin();
		auto v2 = h->End();
		auto m = h->Polygon();
		if (m == NULL)continue;
		auto v0_index = heMesh->Index(m->BoundaryVertice()[0]);
		int order1 = 0;//���v1���������λ��
		int order2 = 0;//���v2���������λ��
		int order3 = 0;//���theta���������λ��
		if (v0_index == heMesh->Index(v1)) {
			order1 = 0;
			order2 = 1;
			order3 = 2;
		}
		else if (v0_index == heMesh->Index(v2)) {
			order1 = 2;
			order2 = 0;
			order3 = 1;
		}
		else {
			order1 = 1;
			order2 = 2;
			order3 = 0;
		}
		auto cotij = cot.at(heMesh->Index(m)).at(order3);
		auto delta1 = x.at(heMesh->Index(m)).at(2 * order1) - x.at(heMesh->Index(m)).at(2 * order2); //��һ������Ĳ�ֵ
		auto delta2 = x.at(heMesh->Index(m)).at(2 * order1 + 1) - x.at(heMesh->Index(m)).at(2 * order2 + 1);//�ڶ�������Ĳ�ֵ
		if (heMesh->Index(v1) != fixed_point1 && heMesh->Index(v1) != fixed_point2) {//����ê��			
			a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2)b1(heMesh->Index(v1)) -= -cotij * b1(heMesh->Index(v2));
			else a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v2)) += -cotij;
			a_pre.at(heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta1;
			a_pre.at(heMesh->Index(v1)).at(2 * nV + nT + heMesh->Index(m)) += -cotij * delta2;
			//���Ƕ�ui1��
			a_pre.at(nV + heMesh->Index(v1)).at(nV + heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2)b1(nV + heMesh->Index(v1)) -= -cotij * b1(nV + heMesh->Index(v2));
			else a_pre.at(nV + heMesh->Index(v1)).at(nV + heMesh->Index(v2)) += -cotij;
			a_pre.at(nV + heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta2;
			a_pre.at(nV + heMesh->Index(v1)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta1;
			//��ui2��
		}
		if (heMesh->Index(v2) != fixed_point1 && heMesh->Index(v2) != fixed_point2) {//����ê��			

			if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2)b1(heMesh->Index(v2)) -= -cotij * b1(heMesh->Index(v1));
			else a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v1)) += -cotij;

			a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v2)) += cotij;
			a_pre.at(heMesh->Index(v2)).at(2 * nV + heMesh->Index(m)) += cotij * delta1;
			a_pre.at(heMesh->Index(v2)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta2;
			//���Ƕ�ui1��

			if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2)b1(nV + heMesh->Index(v2)) -= -cotij * b1(nV + heMesh->Index(v1));
			else a_pre.at(nV + heMesh->Index(v2)).at(nV + heMesh->Index(v1)) += -cotij;
			a_pre.at(nV + heMesh->Index(v2)).at(nV + heMesh->Index(v2)) += cotij;
			a_pre.at(nV + heMesh->Index(v2)).at(2 * nV + heMesh->Index(m)) += cotij * delta2;
			a_pre.at(nV + heMesh->Index(v2)).at(2 * nV + nT + heMesh->Index(m)) += -cotij * delta1;
			//��ui2��
		}
		if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2) {
			b1(2 * nV + heMesh->Index(m)) -= -cotij * delta1 * b1(heMesh->Index(v1));
			b1(2 * nV + heMesh->Index(m)) -= -cotij * delta2 * b1(nV + heMesh->Index(v1));
		}
		else {
			a_pre.at(2 * nV + heMesh->Index(m)).at(heMesh->Index(v1)) += -cotij * delta1;
			a_pre.at(2 * nV + heMesh->Index(m)).at(nV + heMesh->Index(v1)) += -cotij * delta2;
		}
		if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2) {
			b1(2 * nV + heMesh->Index(m)) -= cotij * delta1 * b1(heMesh->Index(v2));
			b1(2 * nV + heMesh->Index(m)) -= cotij * delta2 * b1(nV + heMesh->Index(v2));
		}
		else {
			a_pre.at(2 * nV + heMesh->Index(m)).at(heMesh->Index(v2)) += cotij * delta1;
			a_pre.at(2 * nV + heMesh->Index(m)).at(nV + heMesh->Index(v2)) += cotij * delta2;
		}
		a_pre.at(2 * nV + heMesh->Index(m)).at(2 * nV + heMesh->Index(m)) += cotij * (delta1 * delta1 + delta2 * delta2);
		a_pre.at(2 * nV + heMesh->Index(m)).at(2 * nV + nT + heMesh->Index(m)) += 0.0;
		//��at��
		if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2) {
			b1(2 * nV + nT + heMesh->Index(m)) -= -cotij * delta2 * b1(heMesh->Index(v1));
			b1(2 * nV + nT + heMesh->Index(m)) -= cotij * delta1 * b1(nV + heMesh->Index(v1));
		}
		else {
			a_pre.at(2 * nV + nT + heMesh->Index(m)).at(heMesh->Index(v1)) += -cotij * delta2;
			a_pre.at(2 * nV + nT + heMesh->Index(m)).at(nV + heMesh->Index(v1)) += cotij * delta1;
		}
		if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2) {
			b1(2 * nV + nT + heMesh->Index(m)) -= cotij * delta2 * b1(heMesh->Index(v2));
			b1(2 * nV + nT + heMesh->Index(m)) -= -cotij * delta1 * b1(nV + heMesh->Index(v2));
		}
		else {
			a_pre.at(2 * nV + nT + heMesh->Index(m)).at(heMesh->Index(v2)) += cotij * delta2;
			a_pre.at(2 * nV + nT + heMesh->Index(m)).at(nV + heMesh->Index(v2)) += -cotij * delta1;
		}
		a_pre.at(2 * nV + nT + heMesh->Index(m)).at(2 * nV + heMesh->Index(m)) += 0.0;
		a_pre.at(2 * nV + nT + heMesh->Index(m)).at(2 * nV + nT + heMesh->Index(m)) += cotij * (delta1 * delta1 + delta2 * delta2);
	}
	for (register size_t i = 0; i < n; ++i) {
		for (register size_t j = 0; j < n; ++j) {
			if (a_pre.at(i).at(j) != 0) {
				a.insert(i, j) = a_pre.at(i).at(j);
			}
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
	vector<double> l;
	for (register size_t i = 0; i < nV; ++i) {
		l.push_back(x1(i));
		l.push_back(x1(i + nV));
	}
	for (auto v : heMesh->Vertices()) {
		v->pos.at(0) = x1(heMesh->Index(v));
		v->pos.at(1) = x1(nV + heMesh->Index(v));
		v->pos.at(2) = 0;
	}
}
bool MinSurf_ASAP::Run() {
	//std::cout<<"test"<<std::endl;
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}
	ASAP();

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