#ifndef _MYUTILS_SORT_H_
#define _MYUTILS_SORT_H_

#include <vector>
#include <functional>

namespace myutils {

/*	WARNING: this class has very limited utility.
	
	Syntax:
	sort(sortme.begin(),sortme.end(),sort_by_vector<T>(sortby));

	where sortby is the vector of interest, if sortme is a vector
	that starts of as the indeces of sortby, i.e. 0,1,2,...,size()-1
	then following the sort, it will be reordered according to sortby.
	*/
template<typename T>
class sort_by_vector : public std::binary_function<int,int,bool>
{
	const vector<T> &sort_by;
public:
	sort_by_vector(const vector<T> &sort_by_in) : sort_by(sort_by_in) {}

	bool operator()(int a, int b) const
	{
		return (sort_by.at(a)<sort_by.at(b));
	}
};

};		// namespace myutils

#endif	// _MYUTILS_SORT_H_