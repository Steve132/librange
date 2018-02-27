#ifndef KDFOREST_HPP
#define KDFOREST_HPP

#include "abstract_range_search.hpp"


template<class FloatType>
class kdforest: public abstract_range_search<FloatType>
{
public:
	const size_t num_forests;
private:
	std::vector<std::vector<size_t>> forest_indices;
	std::vector<std::vector<size_t>> forest_dimension_indices;
	std::vector<size_t> forest_sizes;
	size_t max_size;
private:
	static size_t compute_num_forests(size_t n,size_t feature_size,size_t tnum_forests)
	{
		size_t recommended_max_size;
		for(recommended_max_size=0;(size_t(1) << recommended_max_size) <= n;recommended_max_size++);
		recommended_max_size-=1;

		if(tnum_forests==0)
		{
			return feature_size/recommended_max_size + ((feature_size % recommended_max_size) ? 1 : 0);
		}

		size_t requested_max_size=feature_size / tnum_forests;
		
		if(requested_max_size > recommended_max_size)
		{
			throw std::runtime_error("Requested kdtree size is more than the recommended size");
		}
		return feature_size/requested_max_size+((feature_size % requested_max_size) ? 1 : 0);
	}
public:
	kdforest(const std::vector<FloatType>& tdata,size_t tfeature_size,size_t tnum_forests=0):
		abstract_range_search(tdata,tfeature_size),
		num_forests(compute_num_forests(num_features,tfeature_size,tnum_forests)),
		forest_indices(num_forests),
		forest_dimension_indices(num_forests),
		forest_sizes(num_forests)
	{
		std::vector<size_t> sindices(num_features);
		std::iota(sindices.begin(),sindices.end(),0);

		size_t mainsize=feature_size / num_forests;
		size_t lastsize=feature_size % num_forests;
		std::fill(forest_sizes.begin(),forest_sizes.end(),mainsize);
		forest_sizes.back()=lastsize;

		max_size=1;
		
		for(size_t fi=0;fi<num_forests;fi++)
		{
			forest_indices[fi]=sindices;
			kdsort(data+fi*mainsize,forest_sizes[fi],
				forest_indices.begin(),forest_indices.end(),
				forest_dimension_indices.begin(),
				max_size,0,feature_size);
		}
	}
	
	template<class FeatureFunc>
	kdforest(FeatureFunc f,size_t tnum_features,size_t tfeature_size,size_t tnum_forests=0):
		kdforest(build_features_array(f,tnum_features,tfeature_size),tfeature_size,tnum_forests)
	{}

	std::vector<size_t> range_query(const FloatType* lower,const FloatType* upper,std::vector<bool> mask=std::vector<bool>())
	{
		if(mask.size() == 0)
		{
			mask=std::vector<bool>(feature_size,true);
		}
		std::vector<size_t> outrange;
		outrange.reserve(1024);

		bool inited=false;
		std::unordered_set<size_t> candidate(16);

		
		for(size_t fi=0;fi<num_forests;fi++)
		{
			outrange.clear();
			size_t countdown=~size_t(0);
			kdrangesearch(data+fi*forest_sizes[0],forest_sizes[fi],
			forest_indices[fi].cbegin(),forest_indices[fi].cend(),
			forest_dimension_indices[fi].cbegin(),forest_dimension_indices[fi].cend(),
			lower,upper,
			std::vector<bool>(mask.cbegin()+forest_sizes[0]*fi,mask.cbegin()+forest_sizes[0]*fi+forest_sizes.back()),
			std::back_inserter(outrange),countdown,
			max_size,0,feature_size);

			if(!inited)
			{
				candidate.insert(outrange.cbegin(),outrange.cend());
			}
			else
			{
				std::unordered_set<size_t> newc(16);
				for(size_t dex : outrange)
				{
					if(candidate.count(dex))
					{
						newc.insert(dex);
					}
				}
				candidate=newc;
			}
		}
		return std::vector<size_t>(candidate.cbegin(),candidate.cend());
	}
};

#endif