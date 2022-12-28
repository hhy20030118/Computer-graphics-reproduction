#pragma once

#include <Basic/HeapObj.h>
//#include <Engine/Primitive/MassSpring.h>
#include <UGM/UGM>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/Dense>
#include <Eigen/Core>
namespace Ubpa {
	class Simulate : public HeapObj {
	public:
		Simulate(const std::vector<pointf3>& plist,
			const std::vector<unsigned>& elist) {
			edgelist = elist;
			this->positions.resize(plist.size());
			for (int i = 0; i < plist.size(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					this->positions[i][j] = plist[i][j];
				}
			}
		};
	public:
		static const Ptr<Simulate> New(const std::vector<pointf3>& plist,
			const std::vector<unsigned> &elist) {
			return Ubpa::New<Simulate>(plist, elist);
		}
	public:
		// clear cache data
		void Clear();
		int type = 1;
		// init cache data (eg. half-edge structure) for Run()
		bool Init();
		//bool Init();

		// call it after Init()
		bool Run();
		
		const std::vector<pointf3>& GetPositions() const { return positions; };

		const float GetStiff() { return stiff; };
		const float Getks() { return ks; };
		const float Getkd() { return kd; };
		const float Getmess() { return mess; };
		void SetStiff(float k) { stiff = k; };
		void Setks(float k) { ks = k;  };
		void Setkd(float k) { kd = k;  };   //可能会导致发散，暂时不计入
		void Setmess(float k) { mess = k; };
		const float GetTimeStep() { return h; };
		void SetTimeStep(float k) { h = k; Init();};
		std::vector<unsigned>& GetFix() { return this->fixed_id; };
		void SetFix(const std::vector<unsigned>& f) { this->fixed_id = f; Init();};
		const std::vector<pointf3>& GetVelocity() { return velocity; };
		//void SetVelocity(const std::vector<pointf3>& v) { velocity = v; };
		void SetType(int t) { type = t; };
		void SetLeftFix();
		void SetUponFix();
		void SetpointFix();  //设置固定点
		void SetMode();
	private:
		// kernel part of the algorithm
		void SimulateOnce();
		void SimulateOnce_implict();
		void SimulateOnce_implict_fast();
		void SimulateOnce_PBD();
		void Matrix_init();  //因为比较耗时间，有的时候不需要模拟，所以单独放一个位置，开始模拟的时候再调用
		void Matrix_value(Eigen::Matrix3d s, Eigen::MatrixXd* t, int i,int j); //对矩阵赋值（九宫格1赋值，方便起见）
	private:
		int mode = 1;
		Eigen::VectorXd xn= Eigen::Vector<double, Eigen::Dynamic>(), vn= Eigen::Vector<double, Eigen::Dynamic>();  //上一帧的位移速度
		float h = 0.01f;  //步长
		float stiff;
		float ks = 1000.0f, kd = 0.01f;//kd暂时保留
		float mess = 0.1f;
		int times = 6;  //迭代次数
		int divide = 60;  //对于发散，加细步长
		bool isinit = false;
		Eigen::MatrixXd L, J;
		Eigen::SparseMatrix<double> A, B;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverA;
		std::vector<unsigned> fixed_id;  //fixed point id


		//mesh data
		std::vector<unsigned> edgelist;
		std::vector<unsigned> facelist;
		std::vector<pointf3> force;  //受力（半隐式）

		//simulation data
		std::vector<pointf3> positions;
		std::vector<pointf3> positions_next; //下一帧（半隐式）
		std::vector<double> edge_origin;  //弹簧原长
		std::vector<pointf3> velocity;
		
	};
}
