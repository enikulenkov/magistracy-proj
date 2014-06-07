#ifndef CGA_H
#define CGA_H 1

#include <boost/numeric/conversion/cast.hpp>
#include <string>

#include <pagmo/pagmo.h>
//#include <pagmo.h>

class cga : public pagmo::problem::base
{
	public:
		cga(int=1);
		pagmo::problem::base_ptr clone() const;
		std::string get_name() const;
    void objfun_impl_test(pagmo::fitness_vector &, const pagmo::decision_vector &) const;
	protected:
		void objfun_impl(pagmo::fitness_vector &, const pagmo::decision_vector &) const;
	private:
		friend class boost::serialization::access;
		template <class Archive>
		void serialize(Archive &ar, const unsigned int)
		{
			ar & boost::serialization::base_object<base>(*this);
		}
};

BOOST_CLASS_EXPORT_KEY(cga);

#endif /* CGA_H */
