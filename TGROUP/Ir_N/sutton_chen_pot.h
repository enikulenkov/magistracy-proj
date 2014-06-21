#ifndef SUTTON_CHEN_POT_H
#define SUTTON_CHEN_POT_H 1

#include <boost/numeric/conversion/cast.hpp>
#include <string>

#include <pagmo/pagmo.h>
//#include <pagmo.h>

class sutton_chen_pot : public pagmo::problem::base
{
	public:
		sutton_chen_pot(int=1);
		pagmo::problem::base_ptr clone() const;
		std::string get_name() const;
    void objfun_impl_test(pagmo::fitness_vector &, const pagmo::decision_vector &) const;
    void d_objfun(const double *vec, int vec_size, double *df) const;
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

BOOST_CLASS_EXPORT_KEY(sutton_chen_pot);

#endif /* SUTTON_CHEN_POT_H */
