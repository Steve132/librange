//there should be a version of this which uses a linked list or append able vector of features.
//there are m multimaps that map floatvalue->linked list iterator.
//the linked list iterator contains a link to the current multimap iterator for each one if you remove it(not necessary if no removals)
//a search uses multimap.lower_bound and multimap.upper_bound to do the iteration WHILE constructing an unordered_set in place.

#include "abstract_range_search.hpp"
#include<algorithm>

template<class FloatType>
class projection_tree: public abstract_range_search<FloatType>
{
private:
	std::vector<std::vector<size_t>> element_ranges;
private:

	std::function<bool (size_t,size_t)> get_featurecomparator(size_t fi)
	{
		return [&data,feature_size,fi](size_t dex1,size_t dex2)
		{
			FloatType f1=data[dex1*feature_size+fi];
			FloatType f2=data[dex2*feature_size+fi];
			return f1 < f2;
		}
	}
	
	std::pair<size_t,size_t> get_bounds_indices(size_t fi,FloatType fbegin,FloatType fend)
	{
		auto fcfunc=[&data,feature_size,fi](FloatType f1,size_t dex2)
			{
				FloatType f2=data[dex2*feature_size+fi];
				return f1 < f2;
			};
		const std::vector<size_t>& er=element_ranges[fi];
		auto ebegin=std::lower_bound(er.cbegin(),er.cend(),fbegin,fcfunc);
		auto eend=std::upper_bound(er.cbegin(),er.cend(),fend,fcfunc);
		return std::make_pair(ebegin-er.cbegin(),eend-er.cbegin());
	}
	
	template<class Func>
	std::vector<size_t> merge_func_calls(size_t fi,Func f,const std::vector<bool>& mask=std::vector<bool>())
	{
		const std::vector<bool>* maskptr=&mask;
		std::vector<bool> rmo;
		if(mask.size() == 0)
		{
			rmo=std::vector<bool>(feature_size,true);
			maskptr=&rmo;
		}
		
		bool initialized=false;
		std::vector<size_t> fs;
		for(size_t fi=0;fi<feature_size;fi++)
		{
			if(maskptr->operator[fi])
			{
				std::vector<size_t> thfs=f(fi);
				if(initialized)
				{
					std::vector<size_t> copyfs(fs);
					auto last=std::set_intersection(
						copyfs.cbegin(),copyfs.cend(),
						thfs.cbegin(),thfs.cend(),
						fs.begin()
					);
					fs.resize(last-fs.begin());
				}
				else
				{
					fs=thfs;
				}
			}
		}
		return fs;
	}
	
	std::vector<size_t> find_one_in_range(size_t fi,FloatType fbegin,FloatType fend)
	{
		std::pair<size_t,size_t> bi=get_bounds_indices(fi,fbegin,fend);
		std::vector<size_t> set_output(bi.first+element_ranges[fi].cbegin(),eend+element_ranges[fi].cbegin());//should this just be a regular unordered set instead of a realloc and intersect?
		std::sort(set_output.begin(),set_output.end());
		return set_output;
	}
	std::vector<size_t> find_one_nearest(size_t fi,FloatType f)
	{
		std::pair<size_t,size_t> bi=get_bounds_indices(fi,f,f);
		if(bi->first > 0)
		{
			bi->first--;
		}
		if(bi->second < num_features)
		{
			bi->second++;
		}
		std::vector<size_t> set_output(bi.first+element_ranges[fi].cbegin(),eend+element_ranges[fi].cbegin());//should this just be a regular unordered set instead of a realloc and intersect?
		std::sort(set_output.begin(),set_output.end());
		return set_output;
	}
public:
	
	projection_set(const std::vector<FloatType>& tdata,size_t tfeature_size):
		abstract_range_search(tdata,tfeature_size),
		element_ranges(tfeature_size)
	{
		std::vector<size_t> sindices(num_features);
		std::iota(sindices.begin(),sindices.end(),0);
		
		for(size_t fi=0;fi<feature_size;fi++)
		{
			element_ranges[fi]=sindices;
			std::sort(element_ranges[fi].begin(),element_ranges[fi].end(),
				[&data,feature_size,fi](size_t dex1,size_t dex2)
				{
					FloatType f1=data[dex1*feature_size+fi];
					FloatType f2=data[dex2*feature_size+fi];
					return f1 < f2;
				}
			);
		}
	}
	
	template<class FeatureFunc>
	projection_set(FeatureFunc f,size_t tnum_features,size_t tfeature_size):
		projection_set(build_features_array(f,tnum_features,tfeature_size),tfeature_size)
	{}
	

	std::vector<size_t> find_in_range(
		const FloatType* lower,const FloatType* upper,
		const std::vector<bool>& mask=std::vector<bool>())
	{
		return merge_func_calls(fi,
			[lower,upper](const size_t fi2) { return find_one_in_range(fi2,lower[fi2],upper[fi2]);},
		mask
		);
	}
	
	std::vector<size_t> find_nearest(const FloatType* feature,const std::vector<bool>& mask=std::vector<bool>())
	{
		return merge_func_calls(fi,
			[lower,upper](const size_t fi2) { return find_one_nearest(fi2,feature[fi2]);},
		mask
		);
	}
	
	std::vector<size_t> find_near_using_subsearch_percentage(
		const FloatType* feature,
		FloatType subsearch_percentage=0.2,
		const std::vector<bool>& mask=std::vector<bool>()
	)
	{
		std::vector<FloatType> lower(feature_size),upper(feature_size);
		for(size_t fi=0;fi<feature_size;fi++)
		{
			FloatType l=data[element_ranges[fi].front()*feature_size+fi];
			FloatType u=data[element_ranges[fi].back()*feature_size+fi];
			FloatType radius=(l-u)*subsearch_percentage;
			lower=feature[fi]-radius/2.0;
			upper=feature[fi]+radius/2.0;
		}
		return find_in_range(lower,upper,mask);
	}
};
	
