#include "PathTracer.h"

#include <UBL/Image.h>

#include <iostream>

#include <thread>
#include <queue>
using namespace Ubpa;
using namespace std;

#define ALIAS 1
PathTracer::PathTracer(const Scene* scene, const SObj* cam_obj, Image* img, size_t spp)
	: scene{ scene },
	bvh{ const_cast<Scene*>(scene) },
	img{ img },
	cam{ cam_obj->Get<Cmpt::Camera>() },
	ccs{ cam->GenCoordinateSystem(cam_obj->Get<Cmpt::L2W>()->value) },
	spp{ spp }
{
	IntersectorVisibility::Instance();
	IntersectorClosest::Instance();

	scene->Each([this](const Cmpt::Light* light) ->bool {
		if (!vtable_is<EnvLight>(light->light.get()))
			return true; // continue

		env_light = static_cast<const EnvLight*>(light->light.get());
		return false; // stop
	});
	cout << "begin initial" << endl;
	// TODO: preprocess env_light here
	double sum = 0;
	auto w = env_light->texture->img->width;
	auto h = env_light->texture->img->height;
	queue<int> more , less;
	for (register int i = 0; i < w; ++i) {
		for (register int j = 0; j < h; ++j) {
			sum += env_light->texture->img->At(i, j).to_rgb().illumination();
		}
	}

	Ui.resize(w * h);
	Ki.resize(w * h);
	for (register int i = 0; i < w; ++i) {
		for (register int j = 0; j < h; ++j) {
			Ui[i*h+j] = env_light->texture->img->At(i, j).to_rgb().illumination()*w*h/sum;
		}
	}
	for (register int i = 0; i < w*h; ++i) {
		if (abs(Ui[i] - 1.0) < 1e-8) {
			Ki[i] = i;
		}
		else if (Ui[i] > 1.0) {
			more.push(i);
		}
		else {
			less.push(i);
		}
	}
	while (!more.empty() && !less.empty()) {
		auto m = more.front();
		auto l = less.front();
		Ki[l] = m;
		Ui[m] = Ui[m] + Ui[l] - 1;
		less.pop();
		if (abs(Ui[m] - 1.0) < 1e-8) {
			Ki[m] = m;
		}
		else if (Ui[m] < 1.0) {
			more.pop();
			less.push(m);
		}
	}
	while (!more.empty()) {
		auto m = more.front();
		more.pop();
		Ui[m] = 1.0;
	}
	while (!less.empty()) {
		auto l = less.front();
		less.pop();
		Ui[l] = 1.0;
	}
}

void PathTracer::Run() {
	img->SetAll(0.f);
	//const size_t spp = 100; // samples per pixel

#ifdef NDEBUG
	const size_t core_num = std::thread::hardware_concurrency();
	auto work = [this, core_num](size_t id) {
		for (size_t j = id; j < img->height; j += core_num) {
			for (size_t i = 0; i < img->width; i++) {
				for (size_t k = 0; k < spp; k++) {
					float u = (i + rand01<float>() - 0.5f) / img->width;
					float v = (j + rand01<float>() - 0.5f) / img->height;
					rayf3 r = cam->GenRay(u, v, ccs);
					rgbf Lo;
					do { Lo = Shade(IntersectorClosest::Instance().Visit(&bvh, r), -r.dir, true); }
					while (Lo.has_nan());
					img->At<rgbf>(i, j) += Lo / float(spp);
				}
			}
			float progress = (j + 1) / float(img->height);
			cout << progress << endl;
		}
	};
	vector<thread> workers;
	for (size_t i = 0; i < core_num; i++)
		workers.emplace_back(work, i);
	for (auto& worker : workers)
		worker.join();
#else
	for (size_t j = 0; j < img->height; j++) {
		for (size_t i = 0; i < img->width; i++) {
			for (size_t k = 0; k < spp; k++) {
				float u = (i + rand01<float>() - 0.5f) / img->width;
				float v = (j + rand01<float>() - 0.5f) / img->height;
				rayf3 r = cam->GenRay(u, v, ccs);
				rgbf Lo;
				do { Lo = Shade(IntersectorClosest::Instance().Visit(&bvh, r), -r.dir, true); }
				while (Lo.has_nan());
				img->At<rgbf>(i, j) += Lo / static_cast<float>(spp);
			}
		}
		float progress = (j + 1) / float(img->height);
		cout << progress << endl;
	}
#endif
}

rgbf PathTracer::Shade(const IntersectorClosest::Rst& intersection, const vecf3& wo, bool last_bounce_specular) const {
	// TODO: HW9 - Trace
	// [ Tips ]
	// - EnvLight::Radiance(<direction>), <direction> is pointing to environment light
	// - AreaLight::Radiance(<uv>)
	// - rayf3: point, dir, tmin, **tmax**
	// - IntersectorVisibility::Instance().Visit(&bvh, <rayf3>)
	//   - tmin = EPSILON<float>
	//   - tmax = distance to light - EPSILON<float>
	// - IntersectorCloest::Instance().Visit(&bvh, <rayf3>)
	//   - tmin as default (EPSILON<float>)
	//   - tmax as default (FLT_MAX)
	//
	// struct IntersectorClosest::Rst {
	//	 bool IsIntersected() const noexcept { return sobj != nullptr; }
	//	 const SObj* sobj{ nullptr }; // intersection sobj
	//	 pointf3 pos; // intersection point's position
	//	 pointf2 uv; // texcoord
	//   normalf n; // normal, normalized
	//	 vecf3 tangent; // perpendicular to normal, normalized
	// };

	constexpr rgbf error_color = rgbf{1.f,0.f,1.f};
	constexpr rgbf todo_color = rgbf{ 0.f,1.f,0.f };
	constexpr rgbf zero_color = rgbf{ 0.f,0.f,0.f };

	if (!intersection.IsIntersected()) {
		if (last_bounce_specular && env_light != nullptr) {
			// TODO: environment light
			 
			return env_light->Radiance(-wo);
		}
		else
			return zero_color;
	}
	
	if (!intersection.sobj->Get<Cmpt::Material>()) {
		auto light = intersection.sobj->Get<Cmpt::Light>();
		if(!light) return error_color;

		if (last_bounce_specular) { // avoid double-count
			auto area_light = dynamic_cast<const AreaLight*>(light->light.get());
			if (!area_light) return error_color;

			// TODO: area light

			return area_light->Radiance(intersection.uv);
			//return zero_color;
		}else
			return zero_color;
	}
	
	rgbf L_dir{ 0.f };
	rgbf L_indir{ 0.f };

	scene->Each([=, &L_dir](const Cmpt::Light* light, const Cmpt::L2W* l2w, const Cmpt::SObjPtr* ptr) {
		// TODO: L_dir += ...
		// - use PathTracer::BRDF to get BRDF value

		SampleLightResult sample_light_rst = SampleLight(intersection, wo, light, l2w, ptr);
		if (sample_light_rst.pd <= 0)
			return;
		if (sample_light_rst.is_infinity) {
			// TODO: L_dir of environment light
			// - only use SampleLightResult::L, n, pd
			// - SampleLightResult::x is useless
			auto wi = -sample_light_rst.n.cast_to<vecf3>();
			auto l = sqrt(wi[0] * wi[0] + wi[1] * wi[1]+ wi[2] * wi[2]);
			//wi = wi.normalize();
			rayf3 r(intersection.pos,wi);
			bool aviliable = IntersectorVisibility::Instance().Visit(&bvh, r);
			if (aviliable) {
				L_dir += sample_light_rst.L * BRDF(intersection, wi, wo) * abs(wi.cos_theta(intersection.n.cast_to<vecf3>())) * abs(wi.cos_theta(sample_light_rst.n.cast_to<vecf3>())) / l / l / sample_light_rst.pd;
				//L_dir += 10.0;
			}
			/*auto wi = -sample_light_rst.n.cast_to<vecf3>();
			auto r = rayf3(intersection.pos, wi);
			auto v = IntersectorVisibility::Instance().Visit(&bvh, r);
			//auto v = intersectors.visibility.Visit(&bvh, r);
			if (v)
				L_dir += BRDF(intersection, wi, wo) * sample_light_rst.L
				* abs(wi.cos_theta(intersection.n.cast_to<vecf3>()))
				/ sample_light_rst.pd;*/
			
		}
		else {
			// TODO: L_dir of area light
			auto wi = sample_light_rst.x-intersection.pos;
			auto l = sqrt(wi[0] * wi[0] + wi[1] * wi[1] + wi[2] * wi[2]);
			wi = wi.normalize();
			rayf3 r(intersection.pos, wi, EPSILON<float>, l - EPSILON<float>);
			bool aviliable = IntersectorVisibility::Instance().Visit(&bvh, r);
			if (aviliable&& wi.cos_theta(sample_light_rst.n.cast_to<vecf3>())<0.3) {
				L_dir += sample_light_rst.L * BRDF(intersection, wi, wo.normalize()) * abs(wi.cos_theta(intersection.n.cast_to<vecf3>())) * abs(wi.cos_theta(sample_light_rst.n.cast_to<vecf3>())) / l / l / sample_light_rst.pd;
				//L_dir += 10.0;
			}
		}
	});

	// TODO: Russian Roulette
	// - rand01<float>() : random in [0, 1)
	auto rou = rand01<float>();
	const float pdf_RR = 0.9;
	if (rou < pdf_RR) {
		auto sample = SampleBRDF(intersection, wo);
		auto wi = std::get<0>(sample);
		auto pd = std::get<1>(sample);
		if (abs(wi[0]) + abs(wi[1]) + abs(wi[2]) < 0.00001) {
			return L_dir;
		}
		auto wi_p = wi[0] * intersection.n[0] + wi[1] * intersection.n[1] + wi[2] * intersection.n[2];
		if (wi_p <= 0.0)wi = wi - 2.0*wi_p* intersection.n.cast_to<vecf3>();
		wi = wi.normalize();
		rayf3 r(intersection.pos , wi);
			auto y = IntersectorClosest::Instance().Visit(&bvh, r);
			if (y.IsIntersected() && y.sobj->Get<Cmpt::Material>()) {
				auto bounce = true;
				//if (y.sobj->name == "metal_ball")bounce = true;
				L_indir += Shade(y, -wi, bounce) * BRDF(intersection, wi, wo) * abs(wi.cos_theta(intersection.n.cast_to<vecf3>())) / pd / pdf_RR;
			}
		
	}
	// TODO: recursion
	// - use PathTracer::SampleBRDF to get wi and pd (probability density)
	// wi may be **under** the surface
	// - use PathTracer::BRDF to get BRDF value

	// TODO: combine L_dir and L_indir
	return L_dir+L_indir; // you should commemt this line

}

PathTracer::SampleLightResult PathTracer::SampleLight(const IntersectorClosest::Rst& intersection, const vecf3& wo, const Cmpt::Light* light, const Cmpt::L2W* l2w, const Cmpt::SObjPtr* ptr) const {
	PathTracer::SampleLightResult rst;

	auto mat = intersection.sobj->Get<Cmpt::Material>();
	if (!mat) return rst; // invalid
	auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
	if (!brdf) return rst; // not support

	if (wo.dot(intersection.n.cast_to<vecf3>()) < 0)
		return rst;

	rgbf albedo = brdf->Albedo(intersection.uv);
	float metalness = brdf->Metalness(intersection.uv);
	float roughness = brdf->Roughness(intersection.uv);
	//          roughness    0     0.5     1
	// metalness----------------------------
	//     0    |           0.5    0.38    0
	//    0.5   |           0.75   0.56    0
	//     1    |            1     0.75    0
	float p_mat = (1 + metalness) / 2 * (1 - stdBRDF::Alpha(roughness)); // 0 - 1

	auto w2l = l2w->value->inverse();

	float pd_mat, pd_light; // dwi / dA
	vecf3 wi;
	vecf3 light_wi; // wi in light space

	// multi-importance sampling, MIS

	if (vtable_is<AreaLight>(light->light.get())) {
		// [1] area light

		auto area_light = static_cast<const AreaLight*>(light->light.get());
		auto geo = ptr->value->Get<Cmpt::Geometry>();
		if (!geo) return rst; // invalid
		if (!vtable_is<Square>(geo->primitive.get())) return rst; // not support

		rst.n = (l2w->value * normalf{ 0,1,0 }).normalize();
		auto light_p = w2l * intersection.pos; // intersection point's position in light space
		scalef3 world_s = l2w->WorldScale();
		float area = world_s[0] * world_s[1] * Square::area;

		if (rand01<float>() < p_mat) {
			// [1.1] sample material

			// pd_mat : dwi
			tie(wi, pd_mat) = SampleBRDF(intersection, wo);
			light_wi = (w2l * wi).normalize(); // wi in light space

			auto light_r = rayf3{ light_p, light_wi, -std::numeric_limits<float>::max() }; // ray in light space
			auto [isIntersected, t, xz] = light_r.intersect_std_square();
			if (isIntersected) {

				pointf3 p_on_light = pointf3{ xz[0], 0.f, xz[1] };

				pd_light = 1 / area;

				rst.x = l2w->value * p_on_light;
				rst.L = area_light->Radiance({ (xz[0] + 1) / 2, (1 - xz[1]) / 2 });

				// pd_mat : dw -> dA
				float dist2 = light_p.distance2(p_on_light);
				float cos_theta_l = (-light_wi)[1];
				pd_mat *= std::abs(cos_theta_l) / dist2;
			}
			else {
				pd_light = 0.f;
				rst.L = 0.f;
				rst.x = 0.f;
			}
		}
		else {
			// [1.2] sample area light

			auto Xi = uniform_in_square<float>(); // [0, 1] x [0, 1]
			pointf3 p_on_light{ 2 * Xi[0] - 1, 0, 2 * Xi[1] - 1 }; // light space
			vecf3 diff = p_on_light - light_p;
			float dist2 = diff.norm2();
			light_wi = diff / std::sqrt(dist2);
			wi = (l2w->value * light_wi).normalize();

			pd_light = 1.f / area;

			rst.L = area_light->Radiance(Xi.cast_to<pointf2>());
			rst.x = l2w->value * p_on_light;
			rst.n = l2w->UpInWorld().cast_to<normalf>();

			// pd_mat : dw
			matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
			matf3 world_to_surface = surface_to_world.inverse();
			svecf s_wo = (world_to_surface * wo).cast_to<svecf>();
			svecf s_wi = (world_to_surface * wi).cast_to<svecf>();
			pd_mat = brdf->PDF(albedo, metalness, roughness, s_wi, s_wo);

			// pd_mat : dw -> dA
			float cos_theta_l = (-light_wi)[1];
			pd_mat *= std::abs(cos_theta_l) / dist2;
		}
	}
	else if (vtable_is<EnvLight>(light->light.get())) {
		// [2] env light
		auto env_light = static_cast<const EnvLight*>(light->light.get());
		auto light_n = (w2l * intersection.n).normalize(); // intersetion point's normal in light space

		rst.is_infinity = true;
		rst.x = std::numeric_limits<float>::max();
		auto w = env_light->texture->img->width;
		auto h = env_light->texture->img->height;
		auto pi = PI<double>;
		if (rand01<float>() < p_mat) {
			tie(wi, pd_mat) = SampleBRDF(intersection, wo);
			light_wi = (w2l * wi).normalize();
			rst.L = env_light->Radiance(light_wi);
			// pd_light : dwi
			if (!ALIAS) {
				pd_light = env_light->PDF(light_wi, light_n); // TODO: use your PDF			
			}
			else {
				auto texcoord = light_wi.cast_to<normalf>().to_sphere_texcoord();

				auto theta = (1 - texcoord[1]) * pi;
				int i = texcoord[0] * w;
				if (i >= w)i = w - 1;
				int j = texcoord[1] * h;
				if (j >= h)j = h - 1;
				auto pd_img = Ui[i * h + j];
				pd_light = w * h / (2 * pi * pi * sin(theta)) * pd_img;
			}
			
		}
		else {
			// pd_light : dwi
			if (!ALIAS) {
tie(rst.L, light_wi, pd_light) = env_light->Sample(light_n);

			}else{
			// // TODO: use your sampling method
			float x = rand01<float>();
			int i = w * h * x;
			
			float y = w * h * x - i;
			if (abs(x - 1.0) < 1e-5)i--;
			if (i >= h * w) {
				cout << "i";
			}
			int r = 0;
			if (y < Ui[i])r = i;
			else r = Ki[i];
			auto qw = i / h;
			auto qh = i - h * qw;

			float theta = pi *(1 - (float)qh / h);
			float phi = i * 2 * pi / w;

			light_wi = vecf3(sin(theta) * sin(phi), cos(theta), sin(theta) * cos(phi));
			rst.L = env_light->Radiance(light_wi);
			pd_light = Ui[i] * w * h / (2 * pi *pi *sin(theta));

}
			wi = (l2w->value * light_wi).normalize();

			matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
			matf3 world_to_surface = surface_to_world.inverse();
			svecf s_wo = (world_to_surface * wo).cast_to<svecf>();
			svecf s_wi = (world_to_surface * wi).cast_to<svecf>();
			pd_mat = brdf->PDF(albedo, metalness, roughness, s_wi, s_wo);
		}

		rst.n = -wi.cast_to<normalf>();
	}
	else
		return rst; // not support

			
	rst.pd = p_mat * pd_mat + (1 - p_mat) * pd_light;
	
	return rst;
}

std::tuple<vecf3, float> PathTracer::SampleBRDF(const IntersectorClosest::Rst& intersection, const vecf3& wo) {
	auto mat = intersection.sobj->Get<Cmpt::Material>();
	if (!mat) return { vecf3{0.f}, 0.f };
	auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
	if (!brdf) return { vecf3{0.f}, 0.f };

	matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
	matf3 world_to_surface = surface_to_world.inverse();
	svecf s_wo = (world_to_surface * wo).cast_to<svecf>();

	rgbf albedo = brdf->Albedo(intersection.uv);
	float metalness = brdf->Metalness(intersection.uv);
	float roughness = brdf->Roughness(intersection.uv);

	auto [s_wi, pdf] = brdf->Sample(albedo, metalness, roughness, s_wo);
	if (pdf == 0.f)
		return { vecf3{0.f}, 0.f };

	vecf3 wi = surface_to_world * s_wi;

	return { wi,pdf };
}

rgbf PathTracer::BRDF(IntersectorClosest::Rst intersection, const vecf3& wi, const vecf3& wo) {
	auto mat = intersection.sobj->Get<Cmpt::Material>();
	if (!mat) return rgbf{ 1.f,0.f,1.f };
	auto brdf = dynamic_cast<const stdBRDF*>(mat->material.get());
	if (!brdf) return rgbf{ 1.f,0.f,1.f };

	matf3 surface_to_world = svecf::TBN(intersection.n.cast_to<vecf3>(), intersection.tangent);
	matf3 world_to_surface = surface_to_world.inverse();
	svecf s_wi = (world_to_surface * wi).cast_to<svecf>();
	svecf s_wo = (world_to_surface * wo).cast_to<svecf>();

	rgbf albedo = brdf->Albedo(intersection.uv);
	float metalness = brdf->Metalness(intersection.uv);
	float roughness = brdf->Roughness(intersection.uv);

	return brdf->BRDF(albedo, metalness, roughness, s_wi, s_wo);
}
