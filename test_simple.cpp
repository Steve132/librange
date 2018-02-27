#include<iostream>
#include "kdarray.hpp"
#include<array>
#include<random>

void printpoint(const float* l,size_t fs)
{
	std::cout << "{" << l[0];
	for(size_t i=1;i<fs;i++)
	{
		std::cout << "," << l[i];
	}
	std::cout << "}";
}

void simple2dtest(size_t n=12)
{
	std::vector<float> testdata;
	
	for(size_t i=0;i<n;i++)
	{
		testdata.push_back(i);
		testdata.push_back(n-i);
	}

	std::vector<size_t> indices(n);
	std::iota(indices.begin(),indices.end(),0);
	std::vector<size_t> dimsout(n);

	kdsort(&testdata[0],2,
	indices.begin(),indices.end(),
	dimsout.begin());


	for(size_t i=0;i<n;i++)
	{
		std::cout << i << ":" << indices[i] << ":" << dimsout[i] << std::endl;
	}

	std::array<float,2> ll{0.0,5.0},ur{2.2,14.0};
	std::vector<size_t> cand;
	size_t rem=~size_t(0);
	auto outb=std::back_inserter(cand);
	
	kdrangesearch(&testdata[0],2,
	indices.cbegin(),indices.cend(),
	dimsout.cbegin(),dimsout.cend(),
	&ll[0],&ur[0],
	std::vector<bool>(2,true),
	outb,rem);


	for(size_t i=0;i<cand.size();i++)
	{
		std::cout << "#" << cand[i] << std::endl;
	}

}

void simplerandtest(size_t n=100000,size_t d=4)
{
	std::vector<float> testdata(n*d);

	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0,1.0);
	
	std::generate(testdata.begin(),testdata.end(),[&](){ return dis(gen); });	

	std::vector<size_t> indices(n);
	std::iota(indices.begin(),indices.end(),0);
	std::vector<size_t> dimsout(n);

	kdsort(&testdata[0],d,
	indices.begin(),indices.end(),
	dimsout.begin());

	size_t selected=16;
	float epsilon=0.0001f;
	float* ptr=&testdata[d*selected];
	std::cout << "searching for"; printpoint(ptr,d); std::cout << std::endl;
	std::vector<float> ll(ptr,ptr+d),ur(ptr,ptr+d);
	
	for(size_t fi=0;fi<d;fi++)
	{
		ll[fi]-=epsilon;
		ur[fi]+=epsilon;
	}
	
	std::vector<size_t> cand;
	size_t rem=~size_t(0);
	auto outb=std::back_inserter(cand);
	
	kdrangesearch(&testdata[0],d,
	indices.cbegin(),indices.cend(),
	dimsout.cbegin(),dimsout.cend(),
	&ll[0],&ur[0],
	std::vector<bool>(d,true),
	outb,rem);


	for(size_t i=0;i<cand.size();i++)
	{
		std::cout << "#" << cand[i] << std::endl;
	}

}


int main(int argc,char** argv)
{
	simplerandtest(100,2);
	return 0;
}
