#ifndef METRICS_HPP
#define METRICS_HPP

static const size_t MAX_NORM=~size_t(0);
template<class FloatType,size_t P>
struct Lnorm_impl
{
	static FloatType call(const FloatType* f1,const FloatType* f2,size_t n,const std::vector<bool>& mask)
	{
		FloatType sm=0.0;
		for(size_t di=0;di<n;di++)
		{
			if(mask.size() == n && mask[di])
			{			
				FloatType diff=f2[di]-f1[di];
				sm+=std::pow(std::fabs(diff),P);
			}
		}
		return std::pow(sm,1.0/static_cast<double>(P));
	}
};

template<class FloatType>
struct Lnorm_impl<FloatType,MAX_NORM>
{
	static FloatType call(const FloatType* f1,const FloatType* f2,size_t n,const std::vector<bool>& mask)
	{
		FloatType sm=0.0;
		for(size_t di=0;di<n;di++)
		{
			if(mask.size() == n && mask[di])
			{
				FloatType diff=f2[di]-f1[di];
				sm=std::max(sm,std::fabs(diff));
			}
		}
		return sm;
	}
};
template<class FloatType>
struct Lnorm_impl<FloatType,1>
{
	static FloatType call(const FloatType* f1,const FloatType* f2,size_t n,const std::vector<bool>& mask)
	{
		FloatType sm=0.0;
		for(size_t di=0;di<n;di++)
		{
			if(mask.size() == n && mask[di])
			{
				FloatType diff=f2[di]-f1[di];
				sm+=std::abs(diff);
			}
		}
		return sm;
	}
};
template<class FloatType>
struct Lnorm_impl<FloatType,2>
{
	static FloatType call(const FloatType* f1,const FloatType* f2,size_t n,const std::vector<bool>& mask)
	{
		FloatType sm=0.0;
		for(size_t di=0;di<n;di++)
		{
			if(mask.size() == n && mask[di])
			{
				FloatType diff=f2[di]-f1[di];
				sm+=diff*diff;
			}
		}
		return std::sqrt(sm);
	}
};

//premature optimization is the root of all evil but it's possible to implement this recursively and compile-time.

template<class FloatType,size_t P>
inline FloatType Lnorm(const FloatType* f1,const FloatType* f2,size_t n,const std::vector<bool>& mask)
{
	return Lnorm_impl<FloatType,P>::call(f1,f2,n,mask);
}

#endif
