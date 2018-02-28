#ifndef ABSTRACT_RANGE_SEARCH_HPP
#define ABSTRACT_RANGE_SEARCH_HPP

#include<vector>
#include<algorithm>

template<class FloatType>
class abstract_range_search
{
protected:
	std::vector<FloatType> data;
public:
	const size_t num_features;
	const size_t feature_size;

protected:
	template<class FeatureFunc>
	static std::vector<FloatType> build_features_array(FeatureFunc f,size_t n,size_t m)
	{
		std::vector<FloatType> dout(n*m);
		for(size_t i=0;i<n;i++)
		{
			for(size_t j=0;j<m;j++)
			{
				dout[i*m+j]=f(i,j);
			}
		}
		return dout;
	}

	abstract_range_search(const std::vector<FloatType>& tdata,
		size_t tfeature_size):
		data(tdata),
		num_features(tdata.size()/tfeature_size),
		feature_size(tfeature_size)
	{
		//std::vector<size_t> sindices(num_features);
		//std::iota(sindices.begin(),sindices.end(),0);
	}
	
	template<class FeatureFunc>
	abstract_range_search(FeatureFunc f,size_t tnum_features,size_t tfeature_size):
		abstract_range_search(build_features_array(f,tnum_features,tfeature_size),tfeature_size)
	{}

public:
	virtual std::vector<size_t> range_query(
			const FloatType* lower,const FloatType* upper,
		std::vector<bool> mask=std::vector<bool>()) const=0;
	/*std::nth_element to sort by nearest neighbors.
	*/
	std::vector<size_t> reduce_nearest(
		size_t k,
		const FloatType* feature,
		std::vector<size_t> candidates,
		std::function<FloatType (const FloatType*,const FloatType*,size_t,const std::vector<bool>&)> metric,
		const std::vector<bool>& mask=std::vector<bool>()) const
	{
		std::vector<std::pair<size_t,FloatType> > distances2(candidates.size());
		for(size_t i=0;i<candidates.size();i++)
		{
			size_t c=candidates[i];
			distances2[i]=std::make_pair(c,metric(feature,&data[c*feature_size],feature_size,mask));
		}
		if(k < candidates.size())
		{
			std::nth_element(distances2.begin(),distances2.begin()+k,distances2.end(),
			[](const std::pair<size_t,FloatType>& d1,
			   const std::pair<size_t,FloatType>& d2) {return d1.second < d2.second;}
			);
			distances2.resize(k);
			candidates.resize(k);
		}
		std::transform(distances2.cbegin(),distances2.cend(),candidates.begin(),
			       [](const std::pair<size_t,FloatType>& d) { return d.first; }
		     );
		return candidates;
	}

	std::vector<std::size_t> nearest_query(size_t k,const FloatType* feature,
		std::function<FloatType (const FloatType*,const FloatType*,size_t,const std::vector<bool>&)> metric,
		std::vector<bool> mask=std::vector<bool>(),double rstart=0.01,double growthrate=10.0,size_t max_iters=~size_t(0)) const
	{
		if(mask.size()==0)
		{
			mask.resize(feature_size,true);
		}

		//technically the actual epsilon growthrate should be pow(growthrate,1.0/feature_size) for volume awareness.
		//these are doubles for precision.
		//they should also be relative to the size of the space (calculate extrema when building)

		//initial starting value should cover 1% of the data (if 2k/num_features==1%)...use rstart for that.
		//so basically compute this correctly.

		double epsilon=rstart;
		double epsilonrate=std::pow(growthrate,1.0/static_cast<double>(feature_size));

		std::vector<FloatType> upper(feature_size),lower(feature_size);
		std::vector<size_t> found;
		for(size_t current_iter=0;current_iter < max_iters;current_iter++)
		{
			for(size_t di=0;di<feature_size;di++)
			{
				upper[di]=feature[di]+epsilon;
				lower[di]=feature[di]-epsilon;
			}
			found=range_query(lower.data(),upper.data(),mask);
			if(found.size() >= k)
			{
				return reduce_nearest(k,feature,found,metric,mask);
			}
			epsilon*=epsilonrate;
		}
		return std::vector<size_t>();
	}
};

#endif
