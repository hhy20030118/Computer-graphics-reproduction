#include <Engine/MeshEdit/MinSurf.h>

#include <Engine/Primitive/TriMesh.h>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
using namespace Ubpa;

using namespace std;
using namespace Eigen;

MinSurf::MinSurf(Ptr<TriMesh> triMesh)
	: heMesh(make_shared<HEMesh<V>>())
{
	Init(triMesh);
}

void MinSurf::Clear() {
	heMesh->Clear();
	triMesh = nullptr;
}

bool MinSurf::Init(Ptr<TriMesh> triMesh) {
	Clear();

	if (triMesh == nullptr)
		return true;

	if (triMesh->GetType() == TriMesh::INVALID) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is invalid\n");
		return false;
	}
	
	// init half-edge structure
	size_t nV = triMesh->GetPositions().size();
	vector<vector<size_t>> triangles;
	triangles.reserve(triMesh->GetTriangles().size());
	for (auto triangle : triMesh->GetTriangles())
		triangles.push_back({ triangle->idx[0], triangle->idx[1], triangle->idx[2] });
	heMesh->Reserve(nV);
	heMesh->Init(triangles);

	if (!heMesh->IsTriMesh() || !heMesh->HaveBoundary()) {
		printf("ERROR::MinSurf::Init:\n"
			"\t""trimesh is not a triangle mesh or hasn't a boundaries\n");
		heMesh->Clear();
		return false;
	}

	// triangle mesh's positions ->  half-edge structure's positions
	for (int i = 0; i < nV; i++) {
		auto v = heMesh->Vertices().at(i);
		v->pos = triMesh->GetPositions()[i].cast_to<vecf3>();
	}

	this->triMesh = triMesh;
	return true;
}
//using namespace vector;
bool MinSurf::Run() {  //ʵ���ϲ�������
	if (heMesh->IsEmpty() || !triMesh) {
		printf("ERROR::MinSurf::Run\n"
			"\t""heMesh->IsEmpty() || !triMesh\n");
		return false;
	}
	
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
/*inline double MinSurf::distance(int a , int b) {
	return sqrt((heMesh->Vertices()[a]->pos.at(0)- heMesh->Vertices()[b]->pos.at(0))* (heMesh->Vertices()[a]->pos.at(0) - heMesh->Vertices()[b]->pos.at(0))+ (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1)) * (heMesh->Vertices()[a]->pos.at(1) - heMesh->Vertices()[b]->pos.at(1))+ (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)) * (heMesh->Vertices()[a]->pos.at(2) - heMesh->Vertices()[b]->pos.at(2)));
}
inline double MinSurf::acos(int a, int b, int c) {//���ء�A��cosֵ
	auto a1 = distance(a,b);
	auto b1 = distance(a,c);
	auto c1 = distance(b,c);
	return (a1 * a1 + b1 * b1 - c1 * c1) / (2 * a1 * b1);
}*/
/*void MinSurf::one_way_ARAP() {
	//ÿһ�ζ�������
	size_t nV = heMesh->NumVertices();
	SparseMatrix<double> a;
	SimplicialLDLT<SparseMatrix<double>> solver;
	//HouseholderQR<SparseMatrix<double>> solver;;
	VectorXi sizes;
	VectorXd x1 = Vector<double, Dynamic>(), b1 = Vector<double, Dynamic>();
	size_t nT = heMesh->NumPolygons();  //������Էŵ�ͷ�ļ���
	size_t n = 2 * nV;
	//QVector<QVector<double>> x = e;  //�涥��Ⱦ��������ԭʼ����xij //�ȴ�����
	//QVector<QVector<double>> cot = b;//����������cotֵ
	//�Ѿ�����h�ļ�����
	vector<vector<double>> a_pre;//�ݴ����
	vector<double> a_pre_sub;//�ݴ����
	vector<double> at, bt;  //��Lt����
	a.resize(n, n);
	x1.resize(n);
	b1.resize(n);
	sizes.resize(n);
	for (register size_t i = 0; i < n; ++i) {
		sizes(i) = 20;
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

	/*int k = 0;//��һ���ƺ�Ӧ�÷���one����way���Ϊÿ��Ҫ����
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
	*//*
	auto fixed_point1 = heMesh->Index(heMesh->Boundaries()[0][0]->Origin());
	//a_pre[0][1] = 1;
	//���ê��
	a_pre.at(fixed_point1).at(fixed_point1) = 1.0;
	a_pre.at(nV + fixed_point1).at(nV + fixed_point1) = 1.0;
	b1(fixed_point1) = 0.5;
	b1(nV + fixed_point1) = 0.0;
	//���ԸĽ�

	//��Ҫ������󣬶�ÿ��������
	for (auto m : heMesh->Polygons()) {
		Matrix2d l;
		l(0, 0) = 0.0;
		l(1, 0) = 0.0;
		l(0, 1) = 0.0;
		l(1, 1) = 0.0;
		for (register int i = 0; i < 3; ++i) {
			auto order1 = i;
			auto order2 = (i + 1) % 3;
			auto order3 = (i + 2) % 3;
			auto c = cot[heMesh->Index(m)][order3];
			auto delta1_x = x.at(heMesh->Index(m)).at(2 * order1) - x.at(heMesh->Index(m)).at(2 * order2); //��һ������Ĳ�ֵ
			auto delta2_x = x.at(heMesh->Index(m)).at(2 * order1 + 1) - x.at(heMesh->Index(m)).at(2 * order2 + 1);//�ڶ�������Ĳ�ֵ
			auto delta1_u = heMesh->Vertices()[heMesh->Index(m->BoundaryVertice()[order1])]->pos.at(0) - heMesh->Vertices()[heMesh->Index(m->BoundaryVertice()[order2])]->pos.at(0);
			auto delta2_u = heMesh->Vertices()[heMesh->Index(m->BoundaryVertice()[order1])]->pos.at(1) - heMesh->Vertices()[heMesh->Index(m->BoundaryVertice()[order2])]->pos.at(1);
			l(0, 0) += c * delta1_u * delta1_x;
			l(0, 1) += c * delta1_u * delta2_x;
			l(1, 0) += c * delta2_u * delta1_x;
			l(1, 1) += c * delta2_u * delta2_x;
		}
		JacobiSVD<Eigen::MatrixXd> svd(l, ComputeThinU | ComputeThinV);
		Matrix2d U = svd.matrixU();
		Matrix2d V = svd.matrixV();
		Matrix2d S = U * V.transpose();
		at.push_back(S(0, 0));
		bt.push_back(S(0, 1));
		//cout << S(0, 0) << "  " << S(0, 1) << endl;
		//cout << S(1, 0) << "  " << S(1, 1) << endl<<endl;
	}
	//�����Ѿ��õ���Lt����
	for (auto h : heMesh->HalfEdges()) {//�Ծ���ֵ //��asapһ��
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
		if (heMesh->Index(v1) != fixed_point1) {//����ê��			
			a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1)b1(heMesh->Index(v1)) -= -cotij * b1(heMesh->Index(v2));
			else a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v2)) += -cotij;
			//a_pre.at(heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta1;  ����
			b1(heMesh->Index(v1)) -= -cotij * delta1 * at.at(heMesh->Index(m)); //����at
			//a_pre.at(heMesh->Index(v1)).at(2 * nV + nT + heMesh->Index(m)) += -cotij * delta2; ����
			b1(heMesh->Index(v1)) -= -cotij * delta2 * bt.at(heMesh->Index(m));
			//���Ƕ�ui1��
			a_pre.at(nV + heMesh->Index(v1)).at(nV + heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1)b1(nV + heMesh->Index(v1)) -= -cotij * b1(nV + heMesh->Index(v2));
			else a_pre.at(nV + heMesh->Index(v1)).at(nV + heMesh->Index(v2)) += -cotij;
			//a_pre.at(nV + heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta2;
			b1(nV + heMesh->Index(v1)) -= -cotij * delta2 * at.at(heMesh->Index(m));
			//a_pre.at(nV + heMesh->Index(v1)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta1;
			b1(nV + heMesh->Index(v1)) -= cotij * delta1 * bt.at(heMesh->Index(m));
			//��ui2��
		}
		if (heMesh->Index(v2) != fixed_point1) {//����ê��			

			if (heMesh->Index(v1) == fixed_point1)b1(heMesh->Index(v2)) -= -cotij * b1(heMesh->Index(v1));
			else a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v1)) += -cotij;

			a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v2)) += cotij;
			//a_pre.at(heMesh->Index(v2)).at(2 * nV + heMesh->Index(m)) += cotij * delta1;
			b1(heMesh->Index(v2)) -= cotij * delta1 * at.at(heMesh->Index(m));
			//a_pre.at(heMesh->Index(v2)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta2;
			b1(heMesh->Index(v2)) -= cotij * delta2 * bt.at(heMesh->Index(m));
			//���Ƕ�ui1��

			if (heMesh->Index(v1) == fixed_point1)b1(nV + heMesh->Index(v2)) -= -cotij * b1(nV + heMesh->Index(v1));
			else a_pre.at(nV + heMesh->Index(v2)).at(nV + heMesh->Index(v1)) += -cotij;
			a_pre.at(nV + heMesh->Index(v2)).at(nV + heMesh->Index(v2)) += cotij;
			//a_pre.at(nV + heMesh->Index(v2)).at(2 * nV + heMesh->Index(m)) += cotij * delta2;
			b1(nV + heMesh->Index(v2)) -= cotij * delta2 * at.at(heMesh->Index(m));
			//a_pre.at(nV + heMesh->Index(v2)).at(2 * nV + nT + heMesh->Index(m)) += -cotij * delta1;
			b1(nV + heMesh->Index(v2)) -= -cotij * delta1 * bt.at(heMesh->Index(m));
			//��ui2��
		}
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
		vector<double> l;//���ڲ���
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
void MinSurf::ARAP() {
	size_t nV = heMesh->NumVertices();
	size_t nT = heMesh->NumPolygons();  //������Էŵ�ͷ�ļ���
	size_t n = 2 * nV;
	//vector<vector<double>> x;  //�涥��Ⱦ��������ԭʼ����xij
	QVector<float> x_sub;//�ݴ�
	//vector<vector<double>> cot;//����������cotֵ
	QVector<float> cot_sub;//�ݴ�
	for (register int j = 0; j < 6; ++j) {
		x_sub.push_back(0.0);
	}
	for (register size_t i = 0; i < nT; ++i) {
		x.push_back(x_sub);
	}
	for (register int j = 0; j < 3; ++j) {
		cot_sub.push_back(0.0);
	}
	for (register size_t i = 0; i < nT; ++i) {
		cot.push_back(cot_sub);
	}	
	int k = 0;//��һ���ƺ�Ӧ�÷���one����way���Ϊÿ��Ҫ����
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
	ASAP();//��ʼ������
	//Minimize();
	
	int times = 20;
	vector<pointf2> he_pre;
	pointf2 p_pre;
	p_pre.at(0) = 0.0;
	p_pre.at(1) = 0.0;
	for (register size_t i = 0; i < nV; ++i) {
		he_pre.push_back(p_pre);
	}
	fixed_point1 = heMesh->Index(heMesh->Boundaries()[0][0]->Origin());
	fixed_point2 = heMesh->Index(heMesh->Boundaries()[0][1]->Origin());
	for (auto v : heMesh->Vertices()) {
		if (v->IsBoundary() && distance(fixed_point1, heMesh->Index(v)) > distance(fixed_point1, fixed_point2)) {
			fixed_point2 = heMesh->Index(v);
		}
	}
	for (register int j = 0; j < times; ++j) {
		one_way_ARAP();
		bool isfinished = 1;
		for (register size_t i = 0; i < nV; ++i) {
			if (abs(he_pre.at(i).at(0) - heMesh->Vertices()[i]->pos.at(0)) > 0.01)isfinished = 0;
			if (abs(he_pre.at(i).at(1) - heMesh->Vertices()[i]->pos.at(1)) > 0.01)isfinished = 0;
		}
		if (isfinished)break;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
		for (register size_t i = 0; i < nV; ++i) {
			he_pre.at(i).at(0) = heMesh->Vertices()[i]->pos.at(0);
			he_pre.at(i).at(1) = heMesh->Vertices()[i]->pos.at(1);
		}
	}

}*/
/*
void MinSurf::ASAP() {
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
		x[k][4] = distance(heMesh->Index(v1), heMesh->Index(v3))*cos;
		x[k][5] = distance(heMesh->Index(v1), heMesh->Index(v3))*sin;
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
		else{
			order1 = 1;
			order2 = 2;
			order3 = 0;
		}
		auto cotij = cot.at(heMesh->Index(m)).at(order3);
		auto delta1 = x.at(heMesh->Index(m)).at(2 * order1) - x.at(heMesh->Index(m)).at(2 * order2); //��һ������Ĳ�ֵ
		auto delta2 = x.at(heMesh->Index(m)).at(2 * order1+1) - x.at(heMesh->Index(m)).at(2 * order2+1);//�ڶ�������Ĳ�ֵ
		if (heMesh->Index(v1) != fixed_point1&& heMesh->Index(v1) != fixed_point2) {//����ê��			
			a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2)b1(heMesh->Index(v1)) -= -cotij * b1(heMesh->Index(v2));
			else a_pre.at(heMesh->Index(v1)).at(heMesh->Index(v2)) += -cotij;
			a_pre.at(heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta1;
			a_pre.at(heMesh->Index(v1)).at(2 * nV + nT + heMesh->Index(m)) += -cotij * delta2;
			//���Ƕ�ui1��
			a_pre.at(nV+heMesh->Index(v1)).at(nV+heMesh->Index(v1)) += cotij;
			if (heMesh->Index(v2) == fixed_point1 || heMesh->Index(v2) == fixed_point2)b1(nV + heMesh->Index(v1)) -= -cotij * b1(nV+heMesh->Index(v2));
			else a_pre.at(nV + heMesh->Index(v1)).at(nV+heMesh->Index(v2)) += -cotij;
			a_pre.at(nV + heMesh->Index(v1)).at(2 * nV + heMesh->Index(m)) += -cotij * delta2;
			a_pre.at(nV +heMesh->Index( v1)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta1;
			//��ui2��
		}
		if (heMesh->Index(v2) != fixed_point1&& heMesh->Index(v2) != fixed_point2) {//����ê��			
			
			if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2)b1(heMesh->Index(v2)) -= -cotij * b1(heMesh->Index(v1));
			else a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v1)) += -cotij;

			a_pre.at(heMesh->Index(v2)).at(heMesh->Index(v2)) += cotij;
			a_pre.at(heMesh->Index(v2)).at(2 * nV + heMesh->Index(m)) += cotij * delta1;
			a_pre.at(heMesh->Index(v2)).at(2 * nV + nT + heMesh->Index(m)) += cotij * delta2;
			//���Ƕ�ui1��

			if (heMesh->Index(v1) == fixed_point1 || heMesh->Index(v1) == fixed_point2)b1(nV + heMesh->Index(v2)) -= -cotij * b1(nV + heMesh->Index(v1));
			else a_pre.at(nV +heMesh->Index(v2)).at(nV + heMesh->Index(v1)) += -cotij;
			a_pre.at(nV +heMesh->Index( v2)).at(nV + heMesh->Index(v2)) += cotij;
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
			b1(2 * nV + heMesh->Index(m)) -= cotij * delta2 * b1(nV+heMesh->Index(v2));
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
		l.push_back(x1(i+nV));
	}
	for (auto v : heMesh->Vertices()) {
		v->pos.at(0) = x1(heMesh->Index(v));
		v->pos.at(1) = x1(nV+heMesh->Index(v));
		v->pos.at(2) = 0;
	}
}*/
/*void MinSurf::Minimize() {

	size_t nV = heMesh->NumVertices();
	SparseMatrix<double> a;
	SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	VectorXi sizes;
	VectorXd x1 = Vector<double, Dynamic>(), x2 = Vector<double, Dynamic>(), x3 = Vector<double, Dynamic>(), b1 = Vector<double, Dynamic>(), b2 = Vector<double, Dynamic>(), b3 = Vector<double, Dynamic>();
	a.resize(nV,nV);
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
	bool color[100000];//ע�����������ܳ�
	for (register int i = 0; i < 100000; ++i) {
		color[i] = 0;
	}
	for (auto v : heMesh->Vertices()) {  //�ҵ���һ���߽�
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
			if (u->IsBoundary()&&!color[heMesh->Index(u)]) {
				bdy.push_back(heMesh->Index(u));
				current = heMesh->Index(u);
				color[heMesh->Index(u)] = 1;
				finish = 0;
			}
		}
	}
	if (type == 1 && (bdy.size() + 1) >= 4) { //�߽���������
		long q1 = bdy.size() / 4;
		long q2 = bdy.size() / 2;
		long q3 = bdy.size() - q1;
		auto v = heMesh->Vertices()[0];//����һ��߽�
		v->pos.at(0) = 0.0;
		v->pos.at(1) = 0.0;
		v->pos.at(2) = 0.0;
		for (register int i = 1; i < q1; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = (1.0/(double)q1)*i;
			v->pos.at(1) = 0.0;
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q1)];//����һ��߽�
		v->pos.at(0) = 1.0;
		v->pos.at(1) = 0.0;
		v->pos.at(2) = 0.0;
		for (register int i = q1+1; i < q2; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = 1.0;
			v->pos.at(1) = (1.0 / (double)q1) * (i-q1);
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q2)];//����һ��߽�
		v->pos.at(0) = 1.0;
		v->pos.at(1) = 1.0;
		v->pos.at(2) = 0.0;
		for (register int i = q2 + 1; i < q3; ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			v->pos.at(0) = 1.0- (1.0 / (double)q1) * (i - q2);
			v->pos.at(1) = 1.0;
			v->pos.at(2) = 0.0;
		}
		v = heMesh->Vertices()[bdy.at(q3)];//����һ��߽�
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
	
	if (type == 2 && (bdy.size() + 1) >= 4) { //�߽���Բ��
		auto v = heMesh->Vertices()[0];//����һ��߽�
		v->pos.at(0) = 0.5;
		v->pos.at(1) = 1.0;
		v->pos.at(2) = 0.0;
		for (register int i = 1; i < bdy.size(); ++i) {
			v = heMesh->Vertices()[bdy.at(i)];
			double theta = 2 * 3.1415926535 / (double)(bdy.size()) * i;
			v->pos.at(0) = 0.5+0.5*sin(theta);
			v->pos.at(1) = 0.5+0.5*cos(theta);
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
	if (method == 1){
		for (auto v : heMesh->Vertices()) {
			int k = 0;
			for (auto u : v->AdjVertices()) {
				k++;
			}
			for (auto u : v->AdjVertices()) {
				w.at(heMesh->Index(v)).at(heMesh->Index(u)) = 1.0/(double)k;
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
					printf("ERROR::MinSurf::Run\n""\t""fault pointnum\n");
					continue;
				}
				else {
					double a = distance(heMesh->Index(v), heMesh->Index(u));
					double b = distance(heMesh->Index(v), target.at(0));
					double c = distance(target.at(0), heMesh->Index(u));
					double cos1 = (b * b + c * c - a * a) / (2 * b * c);
					double cot1 = cos1 / sqrt(1.0 - cos1 * cos1);
					//if (a > b && a > c)cot1 = -cot1;
					b = distance(heMesh->Index(v), target.at(1));
					c = distance(target.at(1), heMesh->Index(u));
					double cos2 = (b * b + c * c - a * a) / (2 * b * c);
					
					double cot2 = cos2 / sqrt(1.0 - cos2 * cos2);
					//if (a > b && a > c)cot2 = -cot2;
					w.at(heMesh->Index(v)).at(heMesh->Index(u)) = cot1+cot2 ;
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
			
			if(method == 1)a.insert(heMesh->Index(v), heMesh->Index(v)) = num;
			else {a.insert(heMesh->Index(v), heMesh->Index(v)) = sum;}
			b1(heMesh->Index(v)) = 0;
			b2(heMesh->Index(v)) = 0;
			b3(heMesh->Index(v)) = 0;
			for (auto u : v->AdjVertices()) {
				if (u->IsBoundary()) {
					//if (type == 0) {
					if (method == 1) {
						b1(heMesh->Index(v)) += u->pos.at(0);
						b2(heMesh->Index(v)) += u->pos.at(1);
						b3(heMesh->Index(v)) += u->pos.at(2);
					}
					else {
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
	/*cout << "WARNING::MinSurf::Minimize:" << endl
		<< "\t" << "not implemented" << endl;*//*
}*/