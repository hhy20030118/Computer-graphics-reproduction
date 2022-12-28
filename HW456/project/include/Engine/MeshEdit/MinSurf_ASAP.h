#pragma once
#include "MinSurf.h"
namespace Ubpa {
	class TriMesh;
	class Paramaterize;


    class MinSurf_ASAP :
        public MinSurf
    {
	public:
		MinSurf_ASAP(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<MinSurf_ASAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<MinSurf_ASAP>(triMesh);
		}

		bool Run();

	private:
		void ASAP();
		inline double distance(int a, int b);
		inline double acos(int a, int b, int c);
		int type = 3;
		size_t fixed_point1, fixed_point2;
		QVector<QVector<float>> x;
		QVector<QVector<float>> cot;
    };
}

