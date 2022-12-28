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
		void Setkd(float k) { kd = k;  };   //���ܻᵼ�·�ɢ����ʱ������
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
		void SetpointFix();  //���ù̶���
		void SetMode();
	private:
		// kernel part of the algorithm
		void SimulateOnce();
		void SimulateOnce_implict();
		void SimulateOnce_implict_fast();
		void SimulateOnce_PBD();
		void Matrix_init();  //��Ϊ�ȽϺ�ʱ�䣬�е�ʱ����Ҫģ�⣬���Ե�����һ��λ�ã���ʼģ���ʱ���ٵ���
		void Matrix_value(Eigen::Matrix3d s, Eigen::MatrixXd* t, int i,int j); //�Ծ���ֵ���Ź���1��ֵ�����������
	private:
		int mode = 1;
		Eigen::VectorXd xn= Eigen::Vector<double, Eigen::Dynamic>(), vn= Eigen::Vector<double, Eigen::Dynamic>();  //��һ֡��λ���ٶ�
		float h = 0.01f;  //����
		float stiff;
		float ks = 1000.0f, kd = 0.01f;//kd��ʱ����
		float mess = 0.1f;
		int times = 6;  //��������
		int divide = 60;  //���ڷ�ɢ����ϸ����
		bool isinit = false;
		Eigen::MatrixXd L, J;
		Eigen::SparseMatrix<double> A, B;
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solverA;
		std::vector<unsigned> fixed_id;  //fixed point id


		//mesh data
		std::vector<unsigned> edgelist;
		std::vector<unsigned> facelist;
		std::vector<pointf3> force;  //����������ʽ��

		//simulation data
		std::vector<pointf3> positions;
		std::vector<pointf3> positions_next; //��һ֡������ʽ��
		std::vector<double> edge_origin;  //����ԭ��
		std::vector<pointf3> velocity;
		
	};
}
