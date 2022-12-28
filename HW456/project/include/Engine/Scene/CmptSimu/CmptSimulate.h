#pragma once

#include <Engine/Scene/Component.h>
#include <Engine/Scene/CmptSimu/MassSpring.h>
#include <vector>


namespace Ubpa {
	class CmptSimulate : public Component {
	public:
		/*CmptSimulate(
			Ptr<SObj> sobj, Ptr<Simulate> simu,float stiff = 1000.0f)
			: Component(sobj), simulate(simu),stiffness(stiff) { }*/
		CmptSimulate(
			Ptr<SObj> sobj = nullptr, Ptr<Primitive> primitive = nullptr)
			: Component(sobj), primitive(primitive) { }


	public:
		/*static const Ptr<CmptSimulate> New(Ptr<SObj> sobj, Ptr<Simulate> simu, float stiff = 1000.0f) {
			return Ubpa::New<CmptSimulate>(sobj, simu, stiff);
		}*/
		static const  Ptr<CmptSimulate> New(Ptr<SObj> sobj, Ptr<Primitive> primitive) {
			return Ubpa::New<CmptSimulate>(sobj, primitive);
		}
	protected:
		virtual ~CmptSimulate() = default;

	public:
		void Init(std::vector<unsigned>& fix) { SetFix(fix); };
		// theta is in radian

	public:
		const std::vector<unsigned>& GetFix() const { return fix_id; };
		float GetStiff() const { return CastTo<MassSpring>(primitive)->GetSimu()->GetStiff(); };
		float Getks() const { return CastTo<MassSpring>(primitive)->GetSimu()->Getks(); };
		float Getkd() const { return CastTo<MassSpring>(primitive)->GetSimu()->Getkd(); };
		float Getmess() const { return CastTo<MassSpring>(primitive)->GetSimu()->Getmess(); };
		Ptr<Primitive> GetMesh() {
			return primitive;
		}// return CastTo<MassSpring>(primitive); }

	public:
		void SetStiff(float stiff) { CastTo<MassSpring>(primitive)->GetSimu()->SetStiff(stiff); };
		void Setks(float stiff) { CastTo<MassSpring>(primitive)->GetSimu()->Setks(stiff); };
		void Setkd(float stiff) { CastTo<MassSpring>(primitive)->GetSimu()->Setkd(stiff); };
		void Setmess(float stiff) { CastTo<MassSpring>(primitive)->GetSimu()->Setmess(stiff); };
		void SetFix(std::vector<unsigned>& fix) { fix_id = fix; };
		void SetLeftFix() {CastTo<MassSpring>(primitive)->GetSimu()->SetLeftFix();};
		void SetMode() { CastTo<MassSpring>(primitive)->GetSimu()->SetMode(); };
		void SetUponFix() { CastTo<MassSpring>(primitive)->GetSimu()->SetUponFix(); };
		void SetpointFix() { CastTo<MassSpring>(primitive)->GetSimu()->SetpointFix(); };
		void SetType(int t) { CastTo<MassSpring>(primitive)->GetSimu()->SetType(t); };

	private:
		std::vector<unsigned> fix_id;
		float stiffness = 1000;
		Ptr<Primitive> primitive;
	};
}
