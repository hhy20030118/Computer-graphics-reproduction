#pragma once

#include <Basic/HeapObj.h>
#include <Engine/MeshEdit/MinSurf.h>
#include <Engine/MeshEdit/MinSurf_normal.h>
#include <Engine/MeshEdit/MinSurf_ASAP.h>
#include <Engine/MeshEdit/MinSurf_ARAP.h>
#include<vector>
namespace Ubpa {
	class TriMesh;
	class MinSurf;

	class Paramaterize : public HeapObj {
	public:
		Paramaterize(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<Paramaterize> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<Paramaterize>(triMesh);
		}
	private:
		class V;
		class E;
		class P;
		class V : public TVertex<V, E, P> {
		public:
			vecf3 pos;
		};
		class E : public TEdge<V, E, P> { };
		class P :public TPolygon<V, E, P> { };
	private:
		Ptr<MinSurf_normal> surf1;  //三种不同方法
		Ptr<MinSurf_ASAP> surf2;
		Ptr<MinSurf_ARAP> surf3;
		Ptr<TriMesh> triMesh;
		int type = 0;
		int method = 1;
	public:
		void Clear();
		bool Init(Ptr<TriMesh> triMesh);
		bool Run();
		Ptr<TriMesh> PtriMesh;
		//Ptr<HEMesh<V>> PheMesh; // vertice order is same with triMesh
	};
}
