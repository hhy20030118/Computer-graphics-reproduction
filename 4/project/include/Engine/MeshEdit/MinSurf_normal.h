#pragma once
#include "MinSurf.h"
namespace Ubpa {
	class TriMesh;
	class Paramaterize;


	class MinSurf_normal :
		public MinSurf
	{
	public:
		MinSurf_normal(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<MinSurf_normal> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<MinSurf_normal>(triMesh);
		}
		void changetype(int t) { type = t; };   //‘› ±√ª”√
		void changemethod(int t) { method = t; };
		bool Run();

	private:
		inline double distance(int a, int b);
		void Minimize();
		int type = 0;
		int method = 1;
	};
}



