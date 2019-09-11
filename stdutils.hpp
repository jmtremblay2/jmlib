#ifndef JMLIB_STDUTILS
#define JMLIB_STDUTILS

#include <map>
#include <set>
#include <iostream>
#include <vector>

/**
 * Basic operations on STL 
 */
namespace std{

/**
 * check if two vectors have the same size
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return true if v1 and v2 have the same size, false otherwise
 */
template<class T>
bool sameSize(const vector<T>& v1, const vector<T>& v2){
	return v1.size() == v2.size();
}
/**
 * takes a map and returns a set of its keys

 * @param m a map from keys type T to values type V
 *
 * @return a set of the keys in m
 */
template<class T, class V>
set<T> toSet(const map<T,V>& m){
	typename map<T,V>::const_iterator it;
	set<T> s;
	for(it = m.begin(); it != m.end(); ++it)
		s.insert(it->first);
	return s;
}

/**
 * takes a map and returns a set of its values
 *
 * @param m a map
 *
 * @return a set of the values of m
 */
template<class T, class V>
set<V> valuesToSet(const map<T,V>& m){
	typename map<T,V>::const_iterator it;
	set<V> s;
	for(it = m.begin(); it != m.end(); ++it)
		s.insert(it->second);
	return s;
}

/**
 * takes a map and returns a set of its keys
 *
 * @param m a map
 *
 * @return a set of the keys of m
 */
template<class T, class V>
set<T> valuesToSet(const map<T,V>& m){
	typename map<T,V>::const_iterator it;
	set<T> s;
	for(it = m.begin(); it != m.end(); ++it)
		s.insert(it->first);
	return s;
}

/**
 * takes a vector of T and returns a set of its values
 *
 * @param v a vector of Ts
 *
 * @return a set of values in v
 */
template<class T> 
set<T> toSet(const vector<T>& v){
	set<T> s;
	for(int i = 0; i < v.size(); ++i)
		s.insert(v.at(i));
	return s;
	}

/**
 * turns a vector of vectors into a vector of all elements
 *
 * @param x a vector of vector of elements
 *
 * @return a vector of all subvectors concatenated
 */
template<class T>
vector<T> toVector(const vector<vector<T> >& x){
	int size = 0;
	for(unsigned int i = 0; i < x.size(); ++i)
		size += x.at(i).size();
	std::vector<T> v(size);
	int index = 0;
	for(unsigned int i = 0; i < x.size(); ++i)
		for(unsigned int j = 0; j < x.at(i).size(); ++j)
			v.at(index++) = x.at(i).at(j);
	return v;
	
}

/**
 * turns a vector into a vector of vector
 *
 * (cuts a sausage)
 *
 * @param v the vector to be cut
 * @param sizes the vector of size for each piece
 *
 * @return v separated into pieces
 */
template <class T>
vector<vector<T> > toVectorVector(
		const vector<T>& v,
		const vector<int> sizes){
	int totSize = 0;
	for(unsigned int i = 0; i < sizes.size(); ++i)
		totSize += sizes.at(i);
	if((unsigned int) totSize != v.size()) throw string("incorect size or vector length");
	vector<vector<T> > vv(sizes.size());

	int index = 0;
	for(unsigned int i = 0; i < sizes.size(); ++i){
		vv.at(i) = vector<T>(sizes.at(i));
		for(int j = 0; j < sizes.at(i); ++j)
			vv.at(i).at(j) = v.at(index++);
	}
	
	return vv;
}

/**
 * apply sqrt element-wise to vectors
 *
 * @param v a vector
 *
 * @return a vector of the square roots of v
 */
template<class T>
vector<T> sqrt(const vector<T>& v){
	vector<T> ret(v.size());
	for(int i = 0; i < v.size(); ++i)
		ret.at(i) = sqrt(v.at(i));
	return ret;
}

/**
 * apply log element-wise to vectors
 *
 * @param v a vector 
 *
 * @return a vector of the logs of v
 */
template<class T>
vector<T> log(const vector<T>& v){
	vector<T> ret(v.size());
	for(int i = 0; i < v.size(); ++i)
		ret.at(i) = log(v.at(i));
	return ret;
}

/**
 * apply exp element-wise to vectors
 *
 * @param v a vector 
 *
 * @return a vector of the exponentials of v
 */
template<class T>
vector<T> exp(const vector<T>& v){
	vector<T> ret(v.size());
	for(int i = 0; i < v.size(); ++i)
		ret.at(i) = exp(v.at(i));
	return ret;
}

/**
 * calculates the mathematical union of two sets
 *
 * @param s1 a set
 * @param s2 a set
 *
 * @return the union of s1 and s2
 */
template<class T>
set<T> operator+(const set<T>& s1, const set<T>& s2){
	typename set<T>::const_iterator it;
	set<T> newSet(s1);
	for(it = s2.begin(); it != s2.end(); ++it)
		newSet.insert(*it);
	return newSet;
}

/**
 * setUnion of two mathematical set, same as operator +
 *
 * @param s1 a set
 * @param s2 a set
 *
 * @return the mathematical union of s1 and s2
 */
template<class T>
set<T> setUnion(const set<T>& s1, const set<T>& s2){
	return s1 + s2;	
}

/**
 * set inclusion
 *
 * checks if first argument is a subset of second argument
 *
 * @param s1 the set that should be included in the second one
 * @param s2 the set that should include the first one
 *
 * @return true is s1 is a subset of s2, false otherwise
 */
template<class T>
bool operator<=(const set<T>& s1, const set<T>& s2){
	typename set<T>::const_iterator it;
	for(it = s1.begin(); it != s1.end(); ++it){
		if(s2.find(*it) == s2.end())
			return false;
		}
	return true;
}

/**
 * prints the elements of a set
 *
 * @param out an output stream
 * @param s a set
 *
 * @return the output stream passed as argument
 */
template<class T>
ostream& operator<<(ostream& out, const set<T>& s){
	out << "{";
	typename set<T>::const_iterator it;
	for(it = s.begin(); it != s.end(); ++it)
		out << (*it) << ", ";
	out << "}";
	return out;
}


/**
 * prints the (key,value) pairs in a map
 *
 * the key type T and value type U must have an operator<<
 *
 * @param out an output stream
 * @param m a map from keys T to values U
 * @return the output stream passed as argument
 */
template<class T, class U>
ostream& operator<<(ostream& out, const map<T,U>& m){
	typename map<T,U>::const_iterator it;
	out << "{";
	for(it = m.begin(); it != m.end(); ++it)
		out << it->first << ":" << it->second << ", ";
	out << "}";
	return out;
	}

/**
 * prints the elements in a vector
 *
 * the type of vector elements must have a << operator
 *
 * @param out an output stream
 * @param v the vector to print
 *
 * @return the stream
 */
template<class T>
ostream& operator<<(ostream& out, const vector<T>& v){
	typename vector<T>::const_iterator it;
	out << "(";
	for(it = v.begin(); it != v.end(); ++it)
		out << (*it) << ", ";
	out << ")";
	return out;

}

/**
 * element-wise += operator for two vectors
 *
 * @param v1 a vector to modify
 * @param v2 a vector to add
 *
 * @return a reference to the updated v1
 */
template<class T>
vector<T>& operator+=(vector<T>& v1,	const vector<T>& v2){
	if(v1.size() != v1.size())
		throw string("wrong sizes");
	for(int i = 0; i < v1.size(); ++i)
		v1.at(i) += v2.at(i);
	return v1;
	}

/**
 * element-wise += operator for one vector and one element
 *
 * the element t is added to all the elements of v
 *
 * @param v a vector
 * @param t an element
 *
 * @return a reverence to the updated v
 */
template<class T>
vector<T>& operator+=(vector<T>& v, const T& t){
	for(int i = 0; i < v.size(); ++i)
		v.at(i) += t;
	return v;
}

/**
 * element-wise addition of two vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 * 
 * @return element wise addition of v1 and v2
 */
template<class T>
vector<T> operator+(const vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("wrong sizes in vector addition");
	vector<T> ret(v1.size());
	ret += v2;
	return ret;
}

/**
 * element wise addition of a vector and one element
 *
 * @param v a vector
 * @param t an element
 *
 * @return a vector of each element of v minus t
 */
template<class T>
vector<T> operator+(const vector<T>& v, const T& t){
	vector<T> ret(v);
	ret += t;
	return ret;
}

/**
 * element-wise /= of two vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return a rerefence to the updated v1
 */
template<class T>
vector<T>& operator/=(vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("sizes must match for vector /= operator");
	for(int i = 0; i < v1.size(); ++i)
		v1.at(i) /= v2.at(i);
	return v1;
}

/**
 * element-wise /= operator between vector and element
 *
 * @param v a vector
 * @param t an element
 *
 * @return a reference to the updated v
 */
template<class T>
vector<T>& operator/=(vector<T>& v, const T& t){
	for(int i = 0; i < v.size(); ++i)
		v.at(i) /= t;
	return v;
	}	

/**
 * element-wise division of vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 * 
 * @return the element wise division of v1 / v2
 */
template<class T>
vector<T> operator/(const vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("dimension mismatch for vector division");
	vector<T> ret(v1);
	ret /= v2;
	return ret;
}

/** 
 * element-wise division between a vector an one element
 *
 * @param v a vector
 * @param t an element
 *
 * @return each element of v divided by t
 */
template<class T>
vector<T> operator/(const vector<T>& v, const T& t){
	vector<T> ret(v);
	ret /= t;
	return ret;
}

/**
 * element-wise *= for two vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return a reference to the updated v1
 */
template<class T>
vector<T>& operator*=(vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("sizes must match for vector *= operator");
	for(int i = 0; i < v1.size(); ++ i)
		v1.at(i) *= v2.at(i);
	return v1;
}

/**
 * element-wise *= between a vector and one element
 *
 * @param v a vector
 * @param t an element
 *
 * @return each element of v multiplied by t
 */
template<class T>
vector<T>& operator*=(vector<T>& v, const T& t){
	for(int i = 0; i < v.size(); ++i)
		v.at(i) *= t;
	return t;
}

/**
 * element-wise multiplication of vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return the element wise product of x and y
 */
template<class T>
vector<T> operator*(const vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("dimensions of vector multiplication  must match");
	vector<T> ret(v1);
	v1 *= v2;
	return v1;
}
	
/**
 * element wise multiplication of a vector and one element
 *
 * @param v a vector
 * @param t an element
 *
 * @return each element of v time t
 */
template<class T>
vector<T> operator*(const vector<T>& v, const T& t){
	vector<T> ret(v);
	//ret *= t;
	for(int i = 0; i < ret.size(); ++i)
		ret.at(i) *= t;
	return ret;
}

/*
template<class T>
vector<T> operator*(const T& t, const vector<T>& v){
	return v*t;
}
*/
template<class T>
vector<T> operator*(const T& t, const vector<T>& v){
	return v*t;
}

/**
 * element-wise -= of vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return a reference to the updated v1;
 */
template<class T>
vector<T>& operator-=(vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("dimensions of vector -= must match");

	vector<T> ret(v1);
	for(int i = 0; i < ret.size(); ++i)
		ret.at(i) -= v2.at(i);
	return ret;	
}

/**
 * element wise -= of a vector and one element
 *
 * @param v a vector
 * @param t an element
 *
 * @return a reference to the updated v
 */
template<class T>
vector<T>& operator-=(vector<T>& v, const T& t){
	for(int i = 0; i < v.size(); ++i)
		v.at(i) -= t;
	return v;
}

/**
 * element-wise substraction of two vectors
 *
 * @param v1 a vector
 * @param v2 a vector
 *
 * @return the element wise substraction of v1 - v2
 */
template<class T>
vector<T> operator-(const vector<T>& v1, const vector<T>& v2){
	if( ! sameSize(v1,v2))
		throw string("dimensions mismatch for vector substraction");
	vector<T> ret(v1);
	v1 -= v2;
	return v1;
}

/**
 * element wise substraction of a vector and one element
 *
 * @param v a vector
 * @param t an element
 * 
 * @return each elemene of v plus t
 */
template<class T>
vector<T> operator-(const vector<T>& v, const T& t){
	vector<T> ret(v);
	ret -= t;
	return ret;
}

}//namespace
#endif
