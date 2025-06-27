#pragma once

//std
#include <cstdio>
#include <cstdint>

namespace materials
{
	enum class type : uint32_t;
}

namespace materials
{
	class Material
	{
	protected:
		//constructors
		Material(void);

		//destructor
		virtual ~Material(void);

		//serialization
		virtual void load(FILE*);
		virtual void save(FILE*) const;

		//create
		static Material* create(Material*&, const Material*);
		static Material* create(Material*&, materials::type, const Material* = nullptr);

	public:
		//data
		const char* label(void) const;
		const char* label(const char*);

		double specific_mass(double);
		double specific_mass(void) const;

		//type
		virtual materials::type type(void) const = 0;

		//name
		const char* name(void) const;
		static const char* name(materials::type);

		//index
		uint32_t index(void) const;

		//analysis
		virtual bool check(void) const;

	protected:
		//data
		uint32_t m_index;
		char m_label[200];
		double m_specific_mass;
	};
}