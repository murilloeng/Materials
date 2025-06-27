//std
#include <cmath>
#include <cfloat>

//fea
#include "Galileo/fea/inc/Mesh/Materials/Types.hpp"
#include "Galileo/fea/inc/Mesh/Materials/Mechanic/Hyperelastic.hpp"

namespace fea
{
	namespace mesh
	{
		namespace materials
		{
			//constructors
			Hyperelastic2::Hyperelastic2(void) : 
				m_model(materials::model::hencky),
				m_poisson_ratio(0.2), m_elastic_modulus(1e6), m_thermal_coefficient(1e-5)
			{
				return;
			}

			//destructor
			Hyperelastic2::~Hyperelastic2(void)
			{
				return;
			}

			//serialization
			void Hyperelastic2::load(FILE* file)
			{
				uint32_t model;
				Mechanic::load(file);
				if(fscanf(file, "%d %lf %lf %lf", 
					&model, &m_poisson_ratio, &m_elastic_modulus, &m_thermal_coefficient) != 4)
				{
					printf("\tError: Unable to load Hyperelastic material!\n");
					exit(EXIT_FAILURE);
				}
				m_model = materials::model(model);
			}
			void Hyperelastic2::save(FILE* file) const
			{
				Mechanic::save(file);
				fprintf(file, "%d %+.6e %+.6e %+.6e", 
					uint32_t(m_model), m_poisson_ratio, m_elastic_modulus, m_thermal_coefficient);
			}

			//name
			const char* Hyperelastic2::model_name(void) const
			{
				return model_name(m_model);
			}
			const char* Hyperelastic2::model_name(materials::model model)
			{
				const char* names[] = {
					"Gent", "Simple", "Neo Hookean"
				};
				return names[uint32_t(model)];
			}

			//type
			materials::type Hyperelastic2::type(void) const
			{
				return materials::type::hyperelastic;
			}

			//data
			double Hyperelastic2::poisson_ratio(void) const
			{
				return m_poisson_ratio;
			}
			double Hyperelastic2::poisson_ratio(double poisson_ratio)
			{
				return m_poisson_ratio = poisson_ratio;
			}

			double Hyperelastic2::elastic_modulus(void) const
			{
				return m_elastic_modulus;
			}
			double Hyperelastic2::elastic_modulus(double elastic_modulus)
			{
				return m_elastic_modulus = elastic_modulus;
			}

			double Hyperelastic2::thermal_coefficient(void) const
			{
				return m_thermal_coefficient;
			}
			double Hyperelastic2::thermal_coefficient(double thermal_coefficient)
			{
				return m_thermal_coefficient = thermal_coefficient;
			}

			materials::model Hyperelastic2::model(void) const
			{
				return m_model;
			}
			materials::model Hyperelastic2::model(materials::model model)
			{
				return m_model = model;
			}

			//plastic stress
			double Hyperelastic2::plastic_stress(double) const
			{
				return DBL_MAX;
			}
			double Hyperelastic2::plastic_gradient(double) const
			{
				return 0;
			}

			//yield criteria
			double Hyperelastic2::yield_function(const double*, uint32_t) const
			{
				return 0;
			}
			void Hyperelastic2::yield_hessian(double*, const double*, uint32_t) const
			{
				return;
			}
			void Hyperelastic2::yield_gradient(double*, const double*, uint32_t) const
			{
				return;
			}

			//energy
			double Hyperelastic2::energy(const double* s) const
			{
				double (Hyperelastic2::*methods[])(const double*) const = {
					&Hyperelastic2::energy_gent, 
					&Hyperelastic2::energy_hencky, 
					&Hyperelastic2::energy_neo_hookean,
					&Hyperelastic2::energy_saint_venant
				};
				return (this->*methods[uint32_t(m_model)])(s);
			}
			void Hyperelastic2::energy_gradient(double* dW, const double* s) const
			{
				void (Hyperelastic2::*methods[])(double*, const double*) const = {
					&Hyperelastic2::energy_gent_gradient, 
					&Hyperelastic2::energy_hencky_gradient, 
					&Hyperelastic2::energy_neo_hookean_gradient,
					&Hyperelastic2::energy_saint_venant_gradient
				};
				(this->*methods[uint32_t(m_model)])(dW, s);
			}
			void Hyperelastic2::energy_hessian(double* d2W, const double* s) const
			{
				void (Hyperelastic2::*methods[])(double*, const double*) const = {
					&Hyperelastic2::energy_gent_hessian, 
					&Hyperelastic2::energy_hencky_hessian, 
					&Hyperelastic2::energy_neo_hookean_hessian,
					&Hyperelastic2::energy_saint_venant_hessian
				};
				(this->*methods[uint32_t(m_model)])(d2W, s);
			}

			//states
			void Hyperelastic2::plane_stress(double* s, double* dW, double* d2W) const
			{
				s[2] = 1;
				while(true)
				{
					energy_gradient(dW, s);
					energy_hessian(d2W, s);
					if(fabs(dW[2]) < 1e-5)
					{
						d2W[0 + 3 * 0] -= d2W[0 + 3 * 2] * d2W[0 + 3 * 2] / d2W[2 + 3 * 2];
						d2W[1 + 3 * 0] -= d2W[0 + 3 * 2] * d2W[1 + 3 * 2] / d2W[2 + 3 * 2];
						d2W[0 + 3 * 1] -= d2W[0 + 3 * 2] * d2W[1 + 3 * 2] / d2W[2 + 3 * 2];
						d2W[1 + 3 * 1] -= d2W[1 + 3 * 2] * d2W[1 + 3 * 2] / d2W[2 + 3 * 2];
						break;
					}
					s[2] -= dW[2] / d2W[2 + 3 * 2];
				}
			}
			void Hyperelastic2::uniaxial_stress(double* s, double* dW, double* d2W) const
			{
				//data
				s[1] = s[2] = 1;
				double d, ds21, ds31;
				//stress
				while(true)
				{
					energy_gradient(dW, s);
					energy_hessian(d2W, s);
					d = dW[1 + 3 * 1] * dW[2 + 3 * 2] - dW[1 + 3 * 2] * dW[1 + 3 * 2];
					ds21 = -(dW[0 + 3 * 1] * dW[2 + 3 * 2] - dW[0 + 3 * 2] * dW[1 + 3 * 2]) / d;
					ds31 = -(dW[0 + 3 * 2] * dW[1 + 3 * 1] - dW[0 + 3 * 1] * dW[1 + 3 * 2]) / d;
					if(sqrt(dW[1] * dW[1] + dW[2] * dW[2]) < 1e-5)
					{
						dW[0 + 3 * 0] += dW[0 + 3 * 1] * ds21 + dW[0 + 3 * 2] * ds31;
						break;
					}
					s[1] += ds21;
					s[2] += ds31;
				}
			}

			//energy
			double Hyperelastic2::energy_gent(const double*) const
			{
				return 0;
			}
			void Hyperelastic2::energy_gent_hessian(double*, const double*) const
			{
				return;
			}
			void Hyperelastic2::energy_gent_gradient(double*, const double*) const
			{
				return;
			}

			double Hyperelastic2::energy_hencky(const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double e1 = log(s[0]);
				const double e2 = log(s[1]);
				const double e3 = log(s[2]);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				return G * (e1 * e1 + e2 * e2 + e3 * e3) + l / 2 * am * am;
			}
			void Hyperelastic2::energy_hencky_hessian(double* d2W, const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = log(s1);
				const double e2 = log(s2);
				const double e3 = log(s3);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = l / s1 / s2;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = l / s1 / s3;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = l / s2 / s3;
				d2W[0 + 3 * 0] = (2 * G * (1 - e1) + l * (1 - am)) / s1 / s1;
				d2W[1 + 3 * 1] = (2 * G * (1 - e2) + l * (1 - am)) / s2 / s2;
				d2W[2 + 3 * 2] = (2 * G * (1 - e3) + l * (1 - am)) / s3 / s3;
			}
			void Hyperelastic2::energy_hencky_gradient(double* dW, const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = log(s1);
				const double e2 = log(s2);
				const double e3 = log(s3);
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				dW[0] = (2 * G * e1 + l * am) / s1;
				dW[1] = (2 * G * e2 + l * am) / s2;
				dW[2] = (2 * G * e3 + l * am) / s3;
			}

			double Hyperelastic2::energy_neo_hookean(const double*) const
			{
				return 0;
			}
			void Hyperelastic2::energy_neo_hookean_hessian(double*, const double*) const
			{
				return;
			}
			void Hyperelastic2::energy_neo_hookean_gradient(double*, const double*) const
			{
				return;
			}

			double Hyperelastic2::energy_saint_venant(const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				return G * (e1 * e1 + e2 * e2 + e3 * e3) + l / 2 * am * am;
			}
			void Hyperelastic2::energy_saint_venant_hessian(double* d2W, const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				d2W[0 + 3 * 1] = d2W[1 + 3 * 0] = l * s1 * s2;
				d2W[0 + 3 * 2] = d2W[2 + 3 * 0] = l * s1 * s3;
				d2W[1 + 3 * 2] = d2W[2 + 3 * 1] = l * s2 * s3;
				d2W[0 + 3 * 0] = (2 * G + l) * (3 * e1 + 1) + l * (e2 + e3);
				d2W[1 + 3 * 1] = (2 * G + l) * (3 * e2 + 1) + l * (e1 + e3);
				d2W[2 + 3 * 2] = (2 * G + l) * (3 * e3 + 1) + l * (e1 + e2);
			}
			void Hyperelastic2::energy_saint_venant_gradient(double* dW, const double* s) const
			{
				//data
				const double v = m_poisson_ratio;
				const double E = m_elastic_modulus;
				//stretch
				const double s1 = s[0];
				const double s2 = s[1];
				const double s3 = s[2];
				const double e1 = (s1 * s1 - 1) / 2;
				const double e2 = (s2 * s2 - 1) / 2;
				const double e3 = (s3 * s3 - 1) / 2;
				const double am = e1 + e2 + e3;
				//material
				const double G = E / 2 / (1 + v);
				const double l = E * v / (1 + v) / (1 - 2 * v);
				//energy
				dW[0] = (2 * G * e1 + l * am) * s1;
				dW[1] = (2 * G * e2 + l * am) * s2;
				dW[2] = (2 * G * e3 + l * am) * s3;
			}
		}
	}
}