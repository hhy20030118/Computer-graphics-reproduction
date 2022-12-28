#pragma once
#include "MinSurf.h"
namespace Ubpa {
	class TriMesh;
	class Paramaterize;


	class MinSurf_ARAP :
		public MinSurf
	{
	public:
		MinSurf_ARAP(Ptr<TriMesh> triMesh);
	public:
		static const Ptr<MinSurf_ARAP> New(Ptr<TriMesh> triMesh) {
			return Ubpa::New<MinSurf_ARAP>(triMesh);
		}

		bool Run();

	private:
		inline double distance(int a, int b);
		inline double acos(int a, int b, int c);
		void ASAP();
		void ARAP();
		void one_way_ARAP();
		int type = 4;
		size_t fixed_point1, fixed_point2;  //锚点
		QVector<QVector<float>> x;   //存几何信息
		QVector<QVector<float>> cot;
	};
}

