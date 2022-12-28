#include <Engine/MeshEdit/Paramaterize.h>

#include <Engine/MeshEdit/MinSurf.h>
#include <Engine/MeshEdit/MinSurf_normal.h>
#include <Engine/MeshEdit/MinSurf_ASAP.h>
#include <Engine/MeshEdit/MinSurf_ARAP.h>
#include <Engine/Primitive/TriMesh.h>

using namespace Ubpa;

using namespace std;

Paramaterize::Paramaterize(Ptr<TriMesh> triMesh) {
	type = triMesh->boundary_type;
	method = triMesh->method;
	Init(triMesh);
	PtriMesh = triMesh;
}

void Paramaterize::Clear() {
	PtriMesh = nullptr;
}

bool Paramaterize::Init(Ptr<TriMesh> triMesh) {
	if (type < 3) {
		surf1 = MinSurf_normal::New(triMesh);
		surf1->Init(triMesh);
	}
	else if (type == 3) {
		surf2 = MinSurf_ASAP::New(triMesh);
		surf2->Init(triMesh);
	}
	else {
		surf3 = MinSurf_ARAP::New(triMesh);
		surf3->Init(triMesh);
	}
	return true;
}

bool Paramaterize::Run() {
	vector<pointf3> pos;
	int k = 0;
	pointf3 aa;
	aa.at(0) = 0;
	aa.at(1) = 0;
	aa.at(2) = 0;
	if (type < 3) {  //略微繁琐
		for (auto v : surf1->heMesh->Vertices()) {
			pos.push_back(aa);//预先存位置，否则参数化完不知道位置了
			pos.at(k).at(0) = v->pos.at(0);
			pos.at(k).at(1) = v->pos.at(1);
			pos.at(k).at(2) = v->pos.at(2);
			k++;
		}
		surf1->Run();
	}
	else if (type == 3) {
		for (auto v : surf2->heMesh->Vertices()) {
			pos.push_back(aa);
			pos.at(k).at(0) = v->pos.at(0);
			pos.at(k).at(1) = v->pos.at(1);
			pos.at(k).at(2) = v->pos.at(2);
			k++;
		}
		surf2->Run();
	}
	else {
		for (auto v : surf3->heMesh->Vertices()) {
			pos.push_back(aa);
			pos.at(k).at(0) = v->pos.at(0);
			pos.at(k).at(1) = v->pos.at(1);
			pos.at(k).at(2) = v->pos.at(2);
			k++;
		}
		surf3->Run();
	}
	vector<pointf3>& posq = pos;
	PtriMesh->Update(posq);
	vector<pointf2> textpos;
	pointf2 a;
	if (type < 3) {
		for (auto v : surf1->heMesh->Vertices()) {
			a.at(0) = v->pos.at(0);
			a.at(1) = v->pos.at(1);
			textpos.push_back(a);
		}
	}
	else if (type == 3) {
		for (auto v : surf2->heMesh->Vertices()) {
			a.at(0) = v->pos.at(0);
			a.at(1) = v->pos.at(1);
			textpos.push_back(a);
		}
	}
	else {
		for (auto v : surf3->heMesh->Vertices()) {
			a.at(0) = v->pos.at(0);
			a.at(1) = v->pos.at(1);
			textpos.push_back(a);
		}
	}
	vector<pointf2>& texposq = textpos;
	PtriMesh->Update(texposq);
	return true;
}
