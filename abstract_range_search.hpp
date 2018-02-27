#ifndef ABSTRACT_RANGE_SEARCH_HPP
#define ABSTRACT_RANGE_SEARCH_HPP

#include<vector>

template<class FloatType>
class abstract_static_range_search
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
		featuresize(tfeature_size)
	{
		//std::vector<size_t> sindices(num_features);
		//std::iota(sindices.begin(),sindices.end(),0);
	}
	
	template<class FeatureFunc>
	abstract_range_search(FeatureFunc f,size_t tnum_features,size_t tfeature_size):
		absract_range_search(build_features_array(f,tnum_features,tfeature_size),tfeature_size)
	{}

public:
	/*std::nth_element to sort by nearest neighbors.
	*/
	std::vector<size_t> reduce_knearest(
		size_t k,
		const FloatType* feature,
		std::vector<size_t> candidates,
		std::function<FloatType (const FloatType*,const FloatType*,size_t,const std::vector<bool>&)> metric,
		const std::vector<bool>& mask=std::vector<bool>())
	{
		std::vector<FloatType> distances(candidates.size());
		for(size_t i=0;i<candidates.size();i++)
		{
			distances[i]=metric(feature,&data[i*feature_size],feature_size,mask);
		}
		if(k < candidates.size())
		{
			std::nth_element(candidates.begin(),candidates.begin()+k,candidates.end(),
				[&distances](const size_t dex1,const size_t dex2) {return distances[dex1] < distances[dex2];}
			);
			candidates.resize(k);
		}
		return candidates;
	}
};

#endif
