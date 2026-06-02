//std
#include <cmath>

//math
#include "Math/inc/Linear/Vector.hpp"
#include "Math/inc/Solvers/Bisection.hpp"
#include "Math/inc/Miscellaneous/bits.hpp"
#include "Math/inc/Miscellaneous/util.hpp"

//Materials
#include "Materials/inc/Mechanic/Stress.hpp"
#include "Materials/inc/Mechanic/States.hpp"
#include "Materials/inc/Mechanic/Mechanic.hpp"

//static
static double fun_stress(double v, const void** args)
{
	//data
	const uint32_t st = *(uint32_t*) args[0];
	const math::Vector& um = *(math::Vector*) args[1];
	const materials::Mechanic* material = (materials::Mechanic*) args[2];
	//return
	const math::Vector sm = v * um;
	return material->yield_function(sm.data(), st) - material->plastic_stress(0);
}
static void fun_yield_function(double* f, const double* s, const void** args)
{
	//data
	const uint32_t st = *(uint32_t*) args[0];
	const materials::Mechanic* material = (materials::Mechanic*) args[1];
	//yield function
	f[0] = material->yield_function(s, st);
}
static void fun_yield_gradient(double* df, const double* s, const void** args)
{
	//data
	const uint32_t st = *(uint32_t*) args[0];
	const materials::Mechanic* material = (materials::Mechanic*) args[1];
	//yield function
	material->yield_gradient(df, s, st);
}
static void fun_plastic_stress(double* sy, const double* p, const void** args)
{
	//data
	const materials::Mechanic* material = (materials::Mechanic*) args[0];
	//plastic stress
	sy[0] = material->plastic_stress(p[0]);
}

namespace materials
{
	//constructors
	Mechanic::Mechanic(void)
	{
		return;
	}

	//destructor
	Mechanic::~Mechanic(void)
	{
		return;
	}

	//serialization
	void Mechanic::load(FILE* file)
	{
		Material::load(file);
		if(fscanf(file, "%lf", &m_break_strain) != 1)
		{
			printf("\tError: Unable to load Mechanic material!\n");
			exit(EXIT_FAILURE);
		}
	}
	void Mechanic::save(FILE* file) const
	{
		Material::save(file);
		fprintf(file, "%+.6e ", m_break_strain);
	}

	//properties
	bool Mechanic::fracture(void)
	{
		return m_break;
	}
	bool Mechanic::fracture(bool fracture)
	{
		return m_break = fracture;
	}
	
	bool Mechanic::hardening(void)
	{
		return m_hardening;
	}
	bool Mechanic::hardening(bool hardening)
	{
		return m_hardening = hardening;
	}

	bool Mechanic::inelastic(void)
	{
		return m_inelastic;
	}
	bool Mechanic::inelastic(bool inelastic)
	{
		return m_inelastic = inelastic;
	}

	double Mechanic::break_strain(void) const
	{
		return m_break_strain;
	}
	double Mechanic::break_strain(double break_strain)
	{
		return m_break_strain = break_strain;
	}

	double Mechanic::bulk_modulus(void) const
	{
		const double n = poisson_ratio();
		const double E = elastic_modulus();
		return E / (1 - 2 * n) / 3;
	}
	double Mechanic::shear_modulus(void) const
	{
		const double n = poisson_ratio();
		const double E = elastic_modulus();
		return E / (1 + n) / 2;
	}

	double* Mechanic::elastic_stiffness(double* C, uint32_t st) const
	{
		//data
		const double v = poisson_ratio();
		const double E = elastic_modulus();
		const uint32_t ns = math::bit_count(st);
		//stiffness
		if(ns == 6)
		{
			C[3 + 6 * 3] = E / (1 + v) / 2;
			C[4 + 6 * 4] = E / (1 + v) / 2;
			C[5 + 6 * 5] = E / (1 + v) / 2;
			C[3 + 6 * 0] = C[0 + 6 * 3] = 0;
			C[4 + 6 * 0] = C[0 + 6 * 4] = 0;
			C[5 + 6 * 0] = C[0 + 6 * 5] = 0;
			C[3 + 6 * 1] = C[1 + 6 * 3] = 0;
			C[4 + 6 * 1] = C[1 + 6 * 4] = 0;
			C[5 + 6 * 1] = C[1 + 6 * 5] = 0;
			C[3 + 6 * 2] = C[2 + 6 * 3] = 0;
			C[4 + 6 * 2] = C[2 + 6 * 4] = 0;
			C[5 + 6 * 2] = C[2 + 6 * 5] = 0;
			C[4 + 6 * 3] = C[3 + 6 * 4] = 0;
			C[5 + 6 * 3] = C[3 + 6 * 5] = 0;
			C[5 + 6 * 4] = C[4 + 6 * 5] = 0;
			C[0 + 6 * 0] = E * (1 - v) / (1 + v) / (1 - 2 * v);
			C[1 + 6 * 1] = E * (1 - v) / (1 + v) / (1 - 2 * v);
			C[2 + 6 * 2] = E * (1 - v) / (1 + v) / (1 - 2 * v);
			C[1 + 6 * 0] = C[0 + 6 * 1] = E * v / (1 + v) / (1 - 2 * v);
			C[2 + 6 * 0] = C[0 + 6 * 2] = E * v / (1 + v) / (1 - 2 * v);
			C[2 + 6 * 1] = C[1 + 6 * 2] = E * v / (1 + v) / (1 - 2 * v);
			return C;
		}
		else
		{
			elastic_flexibility(C, st);
			math::Matrix(C, ns, ns) = math::Matrix(C, ns, ns).inverse();
			return C;
		}
	}
	double* Mechanic::elastic_flexibility(double* D, uint32_t st) const
	{
		//properties
		const double v = poisson_ratio();
		const double E = elastic_modulus();
		//flexibility
		double Dc[36];
		Dc[0 + 6 * 0] = 1;
		Dc[1 + 6 * 1] = 1;
		Dc[2 + 6 * 2] = 1;
		Dc[3 + 6 * 3] = 2 * (1 + v);
		Dc[4 + 6 * 4] = 2 * (1 + v);
		Dc[5 + 6 * 5] = 2 * (1 + v);
		Dc[3 + 6 * 0] = Dc[0 + 6 * 3] = 0;
		Dc[4 + 6 * 0] = Dc[0 + 6 * 4] = 0;
		Dc[5 + 6 * 0] = Dc[0 + 6 * 5] = 0;
		Dc[3 + 6 * 1] = Dc[1 + 6 * 3] = 0;
		Dc[4 + 6 * 1] = Dc[1 + 6 * 4] = 0;
		Dc[5 + 6 * 1] = Dc[1 + 6 * 5] = 0;
		Dc[3 + 6 * 2] = Dc[2 + 6 * 3] = 0;
		Dc[4 + 6 * 2] = Dc[2 + 6 * 4] = 0;
		Dc[5 + 6 * 2] = Dc[2 + 6 * 5] = 0;
		Dc[4 + 6 * 3] = Dc[3 + 6 * 4] = 0;
		Dc[5 + 6 * 3] = Dc[3 + 6 * 5] = 0;
		Dc[5 + 6 * 4] = Dc[4 + 6 * 5] = 0;
		Dc[1 + 6 * 0] = Dc[0 + 6 * 1] = -v;
		Dc[2 + 6 * 0] = Dc[0 + 6 * 2] = -v;
		Dc[2 + 6 * 1] = Dc[1 + 6 * 2] = -v;
		return Stress(nullptr, st).filter_matrix(D, Dc, E);
	}

	//yield surface
	double Mechanic::yield_stress(const double* u, uint32_t st) const
	{
		//data
		math::solvers::Bisection solver;
		const double sy = plastic_stress(0);
		const uint32_t ns = math::bit_count(st);
		if(math::Vector(u, ns).norm() == 0) return 0;
		const math::Vector um = math::Vector(u, ns).unit();
		//solver
		solver.m_x1 = solver.m_x2 = sy;
		solver.m_tolerance = 1e-5 * sy;
		solver.m_system_2 = fun_stress;
		const void* args[] = {&st, &um, this};
		//setup
		while(true)
		{
			if(fun_stress(solver.m_x1, args) > 0) solver.m_x1 /= 2; else break;
		}
		while(true)
		{
			if(fun_stress(solver.m_x2, args) < 0) solver.m_x2 *= 2; else break;
		}
		//solve
		if(!solver.solve(args))
		{
			printf("Mechanic: yield stress function failed\n");
			exit(EXIT_FAILURE);
		}
		//return
		return solver.m_xs;
	}

	//return mapping
	void Mechanic::stress(double* s, const double* e, const Point& mp) const
	{
		//data
		double C[36];
		const uint32_t st = mp.m_stress_types;
		const uint32_t ns = math::bit_count(st);
		//stress
		elastic_stiffness(C, st);
		if(!m_inelastic)
		{
			math::Vector(s, ns) = math::Matrix(C, ns, ns) * math::Vector(e, ns);
		}
		else
		{
			const math::Vector ep(mp.m_plastic_strain_new, ns);
			math::Vector(s, ns) = math::Matrix(C, ns, ns) * (math::Vector(e, ns) - ep);
		}
	}
	void Mechanic::stiffness(double* C, const double* e, const Point& mp) const
	{
		//size
		const uint32_t st = mp.m_stress_types;
		const uint32_t ns = math::bit_count(st);
		//data
		double da, Kp;
		const math::Vector em(e, ns);
		math::Vector sm(ns), dfs(ns);
		math::Matrix Cm(C, ns, ns), Dm(ns, ns), dfss(ns, ns);
		//stiffness
		elastic_stiffness(C, st);
		if(m_inelastic && mp.m_state_new == state::broken)
		{
			math::Matrix(C, ns, ns).zeros();
		}
		if(m_inelastic && mp.m_state_new == state::plastic)
		{
			//state
			elastic_flexibility(Dm.data(), st);
			const double a_old = mp.m_hardening_old;
			const double a_new = mp.m_hardening_new;
			const math::Vector ep(mp.m_plastic_strain_new, ns);
			//stress
			da = a_new - a_old;
			sm = Cm * (em - ep);
			Cm = (Dm + da * dfss).inverse();
			yield_hessian(dfss.data(), sm.data(), st);
			yield_gradient(dfs.data(), sm.data(), st);
			Kp = m_hardening * plastic_gradient(a_new);
			Cm -= (Cm * dfs).outer() / (Kp + dfs.inner(Cm * dfs));
		}
	}
	void Mechanic::return_mapping(double* s, double* C, const double* e, Point& mp) const
	{
		if(!broken(s, C, e, mp) && !predictor(s, C, e, mp)) corrector(s, C, e, mp);
	}

	//test
	void Mechanic::test_yield_hessian(uint32_t st) const
	{
		//data
		const uint32_t nt = 100000;
		const void* args[] = { &st, this };
		const uint32_t ns = math::bit_count(st);
		math::Matrix s(ns, 1), Ha(ns, ns), Hn(ns, ns);
		//test
		srand(uint32_t(time(nullptr)));
		for(uint32_t i = 0; i < nt; i++)
		{
			s.randu();
			yield_hessian(Ha.data(), s.data(), st);
			math::ndiff(fun_yield_gradient, Hn.data(), s.data(), args, ns, ns, 1e-7);
			if((Ha - Hn).norm() > 1e-5 * (Ha.norm() ? Ha.norm() : 1))
			{
				s.print();
				Ha.print("Ha");
				Hn.print("Hn");
				(Ha - Hn).print("Hr", 1e-5);
				printf("Yield hessian numerical test failed!\n");
				break;
			}
			else
			{
				printf("Yield hessian test: %d ok!\n", i);
			}
		}
	}
	void Mechanic::test_yield_gradient(uint32_t st) const
	{
		//data
		const uint32_t nt = 100000;
		const void* args[] = { &st, this };
		const uint32_t ns = math::bit_count(st);
		math::Vector s(ns), ga(ns), gn(ns), gr(ns);
		//test
		srand(uint32_t(time(nullptr)));
		for(uint32_t i = 0; i < nt; i++)
		{
			s.randu();
			yield_gradient(ga.data(), s.data(), st);
			math::ndiff(fun_yield_function, gn.data(), s.data(), args, 1, ns, 1e-7);
			gr = ga - gn;
			if(gr.norm() > 1e-5)
			{
				ga.print("ga");
				gn.print("gn");
				gr.print("gr", 1e-5);
				printf("Yield gradient numerical test failed!\n");
				break;
			}
			else
			{
				printf("Yield gradient test: %d ok!\n", i);
			}
		}
	}
	void Mechanic::test_plastic_gradient(void) const
	{
		//data
		double p, kpa, kpn;
		const uint32_t nt = 100000;
		const void* args[] = { this };
		//test
		srand(uint32_t(time(nullptr)));
		for(uint32_t i = 0; i < nt; i++)
		{
			p = rand() / RAND_MAX;
			kpa = plastic_gradient(p);
			math::ndiff(fun_plastic_stress, &kpn, &p, args, 1, 1, 1e-5);
			if(fabs(kpa - kpn) > 1e-5)
			{
				printf("Plastic gradient numerical test failed!\n");
				break;
			}
			else
			{
				printf("Plastic gradient test: %d ok!\n", i);
			}
		}
	}

	//return mapping
	bool Mechanic::broken(double* s, double* C, const double* e, Point& mp) const
	{
		//size
		const uint32_t st = mp.m_stress_types;
		const uint32_t ns = math::bit_count(st);
		//data
		math::Vector sm(s, ns);
		math::Matrix Cm(C, ns, ns);
		double& a_new = mp.m_hardening_new;
		const double& a_old = mp.m_hardening_old;
		math::Vector ep_new(mp.m_plastic_strain_new, ns);
		const math::Vector ep_old(mp.m_plastic_strain_old, ns);
		//check
		if(mp.m_state_old == state::broken || (m_break && math::Vector(e, ns).norm() > break_strain()))
		{
			sm.zeros();
			Cm.zeros();
			a_new = a_old;
			ep_new = ep_old;
			mp.m_state_new = state::broken;
			return true;
		}
		//return
		return false;
	}
	bool Mechanic::predictor(double* s, double* C, const double* e, Point& mp) const
	{
		//size
		const uint32_t st = mp.m_stress_types;
		const uint32_t ns = math::bit_count(st);
		//data
		math::Vector sm(s, ns);
		math::Matrix Cm(C, ns, ns);
		const math::Vector em(e, ns);
		//new state
		elastic_stiffness(C, st);
		double& p_new = mp.m_hardening_new;
		math::Vector ep_new(mp.m_plastic_strain_new, ns);
		//old state
		const double& p_old = mp.m_hardening_old;
		const math::Vector ep_old(mp.m_plastic_strain_old, ns);
		//predictor
		p_new = p_old;
		ep_new = ep_old;
		sm = Cm * (em - ep_new);
		mp.m_state_new = state::elastic;
		//return
		return yield_function(s, st) - plastic_stress(m_hardening * p_new) < 0;
	}
	bool Mechanic::corrector(double* s, double* C, const double* e, Point& mp) const
	{
		//size
		const uint32_t st = mp.m_stress_types;
		const uint32_t ns = math::bit_count(st);
		//data
		double dp, ddp, fp, Kp;
		const math::Vector em(e, ns);
		const double sr = plastic_stress(0);
		math::Vector sm(s, ns), rm(ns), dfs(ns);
		math::Matrix Dm(ns, ns), Cm(C, ns, ns), dfss(ns, ns);
		//new state
		elastic_flexibility(Dm.data(), st);
		double& p_new = mp.m_hardening_new;
		math::Vector ep_new(mp.m_plastic_strain_new, ns);
		//old state
		const double& p_old = mp.m_hardening_old;
		const math::Vector ep_old(mp.m_plastic_strain_old, ns);
		//loop
		for(uint32_t i = 0; i < 50; i++)
		{
			//state
			dp = p_new - p_old;
			yield_hessian(dfss.data(), s, st);
			yield_gradient(dfs.data(), s, st);
			//yield
			Cm = (Dm + dp * dfss).inverse();
			rm = em - ep_old - dp * dfs - Dm * sm;
			Kp = m_hardening * plastic_gradient(p_new);
			fp = yield_function(s, st) - plastic_stress(m_hardening * p_new);
			//check
			if(rm.norm() < 1e-5 * em.norm() && fabs(fp) < 1e-5 * sr)
			{
				ep_new = ep_old + dp * dfs;
				mp.m_state_new = state::plastic;
				Cm -= (Cm * dfs).outer() / (Kp + dfs.inner(Cm * dfs));
				return true;
			}
			//increment
			ddp = (fp + dfs.inner(Cm * rm)) / (Kp + dfs.inner(Cm * dfs));
			//update
			p_new += ddp;
			sm += Cm * (rm - ddp * dfs);
		}
		//return
		return false;
	}

	//data
	bool Mechanic::m_break = false;
	bool Mechanic::m_hardening = false;
	bool Mechanic::m_inelastic = false;
}