#ifndef KDARRAY_HPP
#define KDARRAY_HPP

#include<algorithm>
#include<vector>
#include<memory>

template<class FloatType,class IndexIterator,class DimIndexIterator>
void kdsort(const FloatType* data,size_t feature_size,
	IndexIterator ibe,IndexIterator ied,
	DimIndexIterator dbe,size_t max_size=1,size_t level=0,size_t stride=0)
{
	size_t current_dim=level % feature_size;
	size_t n=ied-ibe;
	if(n==0)
	{
		return;
	}
	if(n<=max_size)
	{
		std::fill_n(dbe,n,current_dim);
		return;
	}
	if(stride==0)
	{
		stride=feature_size;
	}

	std::nth_element(ibe,ibe+n/2,ied,[data,feature_size,current_dim,stride](const size_t dex1,const size_t dex2)
	{
		FloatType f1=data[dex1*stride+current_dim];
		FloatType f2=data[dex2*stride+current_dim];
		return f1 < f2;
	});
	*(dbe+n/2)=current_dim;
	kdsort(data,feature_size,ibe,ibe+n/2,dbe,max_size,level+1,stride);
	kdsort(data,feature_size,ibe+n/2+1,ied,dbe+n/2+1,max_size,level+1,stride);
}

inline size_t max_feature_size(size_t n)
{
	size_t mf;
	for(mf=1;n > 0;n>>=1) mf++;
	return mf;
}

void printpoint(const float* l,size_t fs);

template<class FloatType>
inline bool checkinside(const FloatType* point,const FloatType* lower,const FloatType* upper,size_t feature_size,const std::vector<bool>& mask)
{
	/*std::cout << "now checking" << std::endl;
	std::cout << "point:"; printpoint(point,feature_size);
	std::cout << "lower:"; printpoint(lower,feature_size);
	std::cout << "upper:"; printpoint(upper,feature_size);*/

	for(size_t di=0;di<feature_size;di++)
	{
		if(mask[di])
		{
			FloatType fd=point[di];
			if((fd >= upper[di]) || (fd < lower[di]))
			{
				//std::cout << "False" << std::endl;
				return false;
			}
		}
	}
	//std::cout << "True" << std::endl;
	return true;
}

template<class FloatType,class IndexIterator,class OutputIterator,class DimIndexIterator>
void kdrangesearch(const FloatType* data,size_t feature_size,
	IndexIterator ibe,IndexIterator ied,
	DimIndexIterator dbe,DimIndexIterator ded,
	const FloatType* lower,const FloatType* upper,
	const std::vector<bool>& mask,
	OutputIterator& output,size_t& remaining_outputs,size_t max_size=1,size_t level=0,size_t stride=0)
{
	if(remaining_outputs==0) return;
	size_t n=ied-ibe;
	size_t n2=n/2;
	if(n==0) return;
	if(stride==0) stride=feature_size;
	if(n<=max_size)
	{
		std::copy_if(ibe,ibe+n,output,
			[=](size_t idx)
			{
				if(remaining_outputs==0)
					return false;
				return checkinside(&data[stride*idx],lower,upper,feature_size,mask);
			}
		);
		return;
	}
	size_t current_dim=*(dbe+n2);
	
	bool masked=!mask[current_dim];
	
	size_t current_dex=*(ibe+n2);

	const FloatType* tptr=&data[stride*current_dex];
	FloatType midline=tptr[current_dim];
	//std::cout << "midline is:" << midline << std::endl;
	
	if(checkinside(tptr,lower,upper,feature_size,mask))
	{
		*output++=current_dex;
		remaining_outputs--;
	}
	
	if(masked || midline >= lower[current_dim])
	{
		kdrangesearch(data,feature_size,ibe,ibe+n2,dbe,dbe+n2,lower,upper,mask,output,remaining_outputs,max_size,level+1,stride);
	}
	if(masked || midline < upper[current_dim])
	{
		kdrangesearch(data,feature_size,ibe+n2+1,ied,dbe+n2+1,ded,lower,upper,mask,output,remaining_outputs,max_size,level+1,stride);
	}
}
//this is a little esoteric, but basically it remaps the indices with a floating point reordering of the data array such that 
//the feature data is laid out linearly in memory (could be better for cache performance...maybe..
//for converting the data order to match what would be an in-order traversal (cache performance again)

//this should help a lot with cache.
//untested.

template<class FloatType,class IndexIterator,class DimIndexIterator>
void kdrotate(FloatType* data,size_t feature_size,
	IndexIterator ibe,IndexIterator ied,size_t max_size=1,size_t level=0,size_t stride=0,size_t* initial_index=nullptr,FloatType* olddata=nullptr)
{
	size_t n=ied-ibe;
	std::unique_ptr<size_t> initial_index_sp;
	std::unique_ptr<FloatType[]> olddata_sp;

	if(n==0)
	{
		return;
	}
	if(stride==0)
	{
		stride=feature_size;
	}
	if(level==0)
	{
		if(initial_index==nullptr)
		{
			initial_index=new size_t;
			initial_index_sp.reset(initial_index);
			*initial_index=0;
		}
		if(olddata==nullptr)
		{
			olddata=new FloatType[n*feature_size];
			olddata_sp.reset(olddata);
			if(stride==feature_size)
			{
				std::copy(data,data+n*feature_size,olddata);
			}
			else
			{
				for(size_t i=0;i<n;i++)
				{
					std::copy(data+i*stride,data,data+i*stride+feature_size,olddata+i*feature_size);
				}
			}
		}
	}
	auto writebackfunc=[=](size_t dexin)
	{
		size_t dexout=(*initial_index)++;
		std::copy(olddata+dexin*feature_size,olddata+(dexin)*feature_size+feature_size,data+dexout*stride);
		return dexout;
	};
	if(n<=max_size)
	{
		for(size_t i=0;i<n;i++)
		{
			*(ibe+i)=writebackfunc(*(ibe+i));
		}
		return;
	}
	*(ibe+n/2)=writebackfunc(*(ibe+n/2));

	kdsort(data,feature_size,ibe,ibe+n/2,max_size,level+1,stride,initial_index,olddata);
	kdsort(data,feature_size,ibe+n/2+1,max_size,level+1,stride,initial_index,olddata);
}
#endif
