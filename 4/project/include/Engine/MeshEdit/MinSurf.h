#pragma once

#include <Basic/HeapObj.h>
#include <UHEMesh/HEMesh.h>
#include <UGM/UGM>
#include<qvector.h>
#include<vector>
	namespace Ubpa {

		class TriMesh;
		class Paramaterize;

		class MinSurf : public HeapObj {
		public:
			MinSurf(Ptr<TriMesh> triMesh);
		public:
			static const Ptr<MinSurf> New(Ptr<TriMesh> triMesh) {
				return Ubpa::New<MinSurf>(triMesh);
			}
		public:
			// clear cache data
			void Clear();

			// init cache data (eg. half-edge structure) for Run()
			bool Init(Ptr<TriMesh> triMesh);

			// call it after Init()
			bool Run();
			// kernel part of the algorithm
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
			
		public:
			friend class Paramaterize;
			Ptr<TriMesh> triMesh;
			const Ptr<HEMesh<V>> heMesh; // vertice order is same with triMesh
		};
	}

