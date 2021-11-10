/*
 * parameters.h
 *
 *  Created on: 2016年3月29日
 *      Author: hyliu
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_
#include "dataio/inputlines.h"
#include "dstl/stlutil.h"
#include "dstl/globalinstances.h"
#include <string>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace NSPdataio {

/*!This class is usually used to store controlling paratemers of a given type.
 * An instance stores a set of parameters as key:value pairs
 * The key is of type std::string
 * The value is of type PARA
 */
template<typename PARA>
class Parameters {
public:
	PARA operator[](const std::string & key) {
		return parameters_.at(key);
	}
	//   void addparameter(const std::string & key, PARA val) {
	//       parameters_.insert(std::make_pair(key,val));
	//   }

	/*! reads the key value pair from a string of the form "key=value"
	 * Does not do anything if the input string does not contain "=".
	 * Whitespace characters surrounding the "=" will be removed.
	 * If the key has already been stored, the new value replaces the previous value.
	 * Otherwise a new key:value entry is inserted in the stored set.
	 */
	void parsetxt(const std::string & line) {
		int key_len = line.find_first_of('=');
		if (key_len == std::string::npos)
			return;
		//       key_len +=1;
		std::string key = boost::trim_copy(line.substr(0, key_len));
		PARA val;
		getparametervalue(line.substr(key_len + 1), &val);
		if (defined(key)) {
			replace(key, val);
		} else {
			insert(key, val);
		}
	}

	/*!insert a new entry in the stored set
	 *
	 */
	void insert(const std::string &key, PARA val) {
#ifdef _OPENMP
		omp_set_nest_lock(&lock_);
#endif
		parameters_.insert(std::make_pair(key, val));
#ifdef _OPENMP
		omp_unset_nest_lock(&lock_);
#endif
	}

	/*!replace the value part of an existing entry
	 *
	 */
	void replace(const std::string &key, PARA val) {
#ifdef _OPENMP
		omp_set_nest_lock(&lock_);
#endif
		parameters_.erase(key);
		insert(key, val);
#ifdef _OPENMP
		omp_unset_nest_lock(&lock_);
#endif
	}

	/*!Check if the stored set contains an entry with the given key
	 *
	 */
	bool defined(const std::string & key) {
		return parameters_.find(key) != parameters_.end();
	}

	/*!Print stored key:value entries to std::out.
	 * Each entry occupies one line.
	 * The format is "key""\t=""value
	 */
	void print(std::ostream &ofs = std::cout) {
		for (auto &pair : parameters_) {
			ofs << pair.first << "\t=";
			print(pair.second);
			ofs << std::endl;
		}
	}

private:
	/*!the stored set of parameters
	 *
	 */
	std::map<std::string, PARA> parameters_;
#ifdef _OPENMP
	omp_nest_lock_t lock_;
#endif
	/*!transform a string into a value of type T using boost::lexical_cast
	 *
	 */
	template<typename T>
	void getparametervalue(const std::string & str, T * val) {
		//      std::cout <<"val: " <<str <<std::endl;
		(*val) = boost::lexical_cast<T>(boost::trim_copy(str));
	}

	/*!removing and leading and pending white spaces from a string
	 *
	 */
	void getparametervalue(const std::string & str, std::string * val) {
		//      std::cout <<"val: " <<str <<std::endl;
		(*val) = boost::trim_copy(str);
	}

	/*! transform a string with usual delimiters to a list of numerical values.
	 *
	 */
	template<typename T>
	void getparametervalue(const std::string &str, std::vector<T> *val) {
		std::vector<std::string> words(parseline(str, std::vector<int>()));
		(*val) = wordsToNumbers<T>(words);
	}

	/*! transform a string with usual delimiters to a list of delimiter-free strings.
	 *
	 */
	void getparametervalue(const std::string &str,
			std::vector<std::string> *val) {
		(*val) = parseline(str, std::vector<int>());
	}

	/*! Print a list of values to std::cout with "\t" as delimiter
	 *
	 */
	template<typename T>
	void print(const std::vector<T> & val,std::ostream &ofs=std::cout) {
		for (auto & v : val) {
			ofs << "\t" << v;
		}
	}

	/*!print a single value to std::cout with a leading "\t".
	 *
	 */
	template<typename T>
	void print(T val,std::ostream &ofs=std::cout) {
		ofs << "\t" << val;
	}
};

/*! This class combines several Parameters objects to store controlling parameters of different types.
 * As in Parameters, each parameter is a key:value pair.
 * The key is of type std::string
 * The value can be of any of the following types:
 *  int, double, long, string
 * Besides single values, the value identified by a key can also be a
 * std::vector that contains multiple values.
 */
struct ParameterSet {
	/*!Parameters objects that store the key:value pairs of respective types
	 *
	 */
	Parameters<int> IntPar;
	Parameters<double> DoublePar;
	Parameters<long> LongPar;
	Parameters<std::vector<int>> IntVecPar;
	Parameters<std::vector<long>> LongVecPar;
	Parameters<std::vector<double>> DoubleVecPar;
	Parameters<std::string> StringPar;
	Parameters<std::vector<std::string>> StringVecPar;

	/*!Known keys or valid keys
	 *These are usually used to check if a key identifier read from an input file is
	 *a valid one.
	 */
	std::set<std::string> IntKeys;
	std::set<std::string> LongKeys;
	std::set<std::string> LongVecKeys;
	std::set<std::string> DoubleKeys;
	std::set<std::string> IntVecKeys;
	std::set<std::string> DoubleVecKeys;
	std::set<std::string> StringKeys;
	std::set<std::string> StringVecKeys;

	/*!an initialization flag
	 */
	bool initialized { false };

	/*!Whether the read of unknown keys causes program abortion
	 *
	 */
	bool readunknownabort { true };

	/*!Define valid keys and the type of values associated with a key.
	 */
	void InitIntKeys(const std::vector<std::string> & keys) {
		for (auto &key : keys) {
			IntKeys.insert(key);
		}
	}
	void InitLongKeys(const std::vector<std::string> &keys) {
		for (auto &key : keys) {
			LongKeys.insert(key);
		}
	}
	void InitLongVecKeys(const std::vector<std::string> & keys) {
		for (auto &key : keys) {
			LongVecKeys.insert(key);
		}
	}
	void InitIntVecKeys(const std::vector<std::string> & keys) {
		for (auto &key : keys) {
			IntVecKeys.insert(key);
		}
	}
	void InitDoubleKeys(const std::vector<std::string> & keys) {
		for (auto &key : keys) {
			DoubleKeys.insert(key);
		}
	}
	void InitDoubleVecKeys(const std::vector<std::string> & keys) {
		for (auto &key : keys) {
			DoubleVecKeys.insert(key);
		}
	}
	void InitStringKeys(const std::vector<std::string> &keys) {
		for (auto &key : keys) {
			StringKeys.insert(key);
		}
	}
	void InitStringVecKeys(const std::vector<std::string> &keys) {
		for (auto &key : keys) {
			StringVecKeys.insert(key);
		}
	}

	/*!Store a key:value pair for value type int
	 *
	 */
	void insertparameter(const std::string &str, int val) {
		IntKeys.insert(str);
		IntPar.insert(str, val);
	}

	/*!store a key:value pair for value type double
	 *
	 */
	void insertparameter(const std::string &str, double val) {
		DoubleKeys.insert(str);
		DoublePar.insert(str, val);
	}
	void insertparameter(const std::string &str, long val) {
		LongKeys.insert(str);
		LongPar.insert(str, val);
	}
	void insertparameter(const std::string &str, std::string val) {
		StringKeys.insert(str);
		StringPar.insert(str, val);
	}

	/*!Check if a value of some type have already been stored for the enquired key
	 *
	 */
	bool keydefined(const std::string & key) {
		return (IntPar.defined(key) || IntVecPar.defined(key)
				|| DoublePar.defined(key) || DoubleVecPar.defined(key)
				|| LongPar.defined(key) || LongVecPar.defined(key)
				|| StringPar.defined(key) || StringVecPar.defined(key));
	}

	/*!Read and parse an input string of the form "key=value".
	 * The key contained in the string is compared to known valid keys.
	 * If it is valid for a certain value type,
	 * the corresponding key:value pair is stored, and the function returns true.
	 * Otherwise nothing is stored and the function returns false.
	 */
	bool readline(const std::string & line) {
		int key_len = line.find_first_of('=');
		if (key_len == std::string::npos)
			return false;
//        key_len +=1;
		std::string key = boost::trim_copy(line.substr(0, key_len));
		bool success = true;
		if (IntKeys.find(key) != IntKeys.end()) {
			IntPar.parsetxt(line);
		} else if (IntVecKeys.find(key) != IntVecKeys.end()) {
			IntVecPar.parsetxt(line);
		} else if (DoubleKeys.find(key) != DoubleKeys.end()) {
			DoublePar.parsetxt(line);
		} else if (DoubleVecKeys.find(key) != DoubleVecKeys.end()) {
			DoubleVecPar.parsetxt(line);
		} else if (StringKeys.find(key) != StringKeys.end()) {
			StringPar.parsetxt(line);
		} else if (StringVecKeys.find(key) != StringVecKeys.end()) {
			StringVecPar.parsetxt(line);
		} else if (LongKeys.find(key) != LongKeys.end()) {
			LongPar.parsetxt(line);
		} else if (LongVecKeys.find(key) != LongVecKeys.end()) {
			LongVecPar.parsetxt(line);
		} else {
			success = false;
		}
		return success;
	}

	/*!read a set of parameters from an input file.
	 * Each line of the input file is expected to be of the form "key=value".
	 * The key should be valid, i.e.,they should be one of those
	 * predefined in the program.
	 *
	 */
	void readparameters(const std::string &filename) {
		TextLines txts;
		std::cout << "Reading parameters from file " << filename << std::endl;
		txts.init(filename);
		for (auto &line : txts.lines()) {
			if (!readline(line)) {
				std::cout << "uninterpretable line encountered: " << std::endl;
				std::cout << line << std::endl;
			}
		}
	}

	void undefinedkeyerror(const std::string &key) {
		std::cout << "Key " << key << " undefined." << std::endl;
		abort();
	}

	/*!Access the value from the key
	 *
	 */
	int getval(const std::string &key, int *val) {
		if (IntKeys.find(key) != IntKeys.end()) {
			*val = IntPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	int getval(const std::string &key, long *val) {
		if (LongKeys.find(key) != LongKeys.end()) {
			*val = LongPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	double getval(const std::string &key, double *val) {
		if (DoubleKeys.find(key) != DoubleKeys.end()) {
			*val = DoublePar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	std::string getval(const std::string &key, std::string *val) {
		if (StringKeys.find(key) != StringKeys.end()) {
			*val = StringPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}

	std::vector<int> getval(const std::string &key, std::vector<int> *val) {
		if (IntVecKeys.find(key) != IntVecKeys.end()) {
			*val = IntVecPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	std::vector<long> getval(const std::string &key, std::vector<long> *val) {
		if (LongVecKeys.find(key) != LongVecKeys.end()) {
			*val = LongVecPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	std::vector<double> getval(const std::string &key,
			std::vector<double> *val) {
		if (DoubleVecKeys.find(key) != DoubleVecKeys.end()) {
			*val = DoubleVecPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	std::vector<std::string> getval(const std::string &key,
			std::vector<std::string> *val) {
		if (StringVecKeys.find(key) != StringVecKeys.end()) {
			*val = StringVecPar[key];
			return *val;
		} else
			undefinedkeyerror(key);
	}
	/*!Print the key:value entries to std::cout
	 *
	 */
	void printparameters(std::ostream &ofs=std::cout) {
		IntPar.print(ofs);
		IntVecPar.print(ofs);
		DoublePar.print(ofs);
		DoubleVecPar.print(ofs);
		StringPar.print(ofs);
		StringVecPar.print(ofs);
	}
};
/**This class wraps ParameterSet with the so-called "singleton" program pattern.
 * That is, a global, single instance of the ParameterSet class is defined and can be accessed
 * from any place of the code.
 * This can be used to create a master set of controlling parameters for a program.
 */
class Controls {
public:
	/*! get a reference to the "singleton" ParameterSet
	 */
	static ParameterSet & getparameterset() {
		static ParameterSet parameterset_;
		return parameterset_;
	}
	/*!get the value for a given key
	 *
	 */
	template<typename VAL>
	static VAL getpar(const std::string &key, VAL*val) {
		return getparameterset().getval(key, val);
	}

	/*!Cannot copy the singleton
	 *
	 */
	Controls(Controls const&) = delete;
	void operator=(Controls const&) = delete;
private:
	Controls() {
		;
	}
};

/*!This class wraps the ParameterSet as "typed" singletons.
 * Namely, different singleton objects may be defined with different T
 * In this way, different types of controls can be defined.
 * Each type of control has its own set of known keys and value types.
 * it can be used to implement different controlling parameter sets,
 * each applies to a distinct functionality.
 */
template<typename T>
class TypedControls {
public:
	static ParameterSet & getparameterset() {
		static ParameterSet parameterset_;
		return parameterset_;
	}
	template<typename VAL>
	static VAL getpar(const std::string &key, VAL*val) {
		return getparameterset().getval(key, val);
	}
	TypedControls(TypedControls<T> const&) = delete;
	void operator=(TypedControls<T> const&) = delete;
private:
	TypedControls() {
		;
	}
};

typedef NSPdstl::GlobalInstances<ParameterSet> MultiInstanceParameterSet;

/*! The GlobalInstances<> template is for
 *  creating multiple globally accessible instances of a class.
 *  Different global instance have different name identifiers  of type std::string.
 *  They are stored in a static map object.
 *  Different control sets of the same type (i.e., of identical set of known keys and value
 *  types) could be defined and used in different parts of the program.
 *
 */
template<typename T>
class TypedMultiInstanceControls: public MultiInstanceParameterSet {
public:
	/*!Get a reference of either the global singleton instance or one of the named instances.
	 *
	 */
	static ParameterSet &getparameterset(const std::string & name =
			std::string()) {
		return getinstance(name);
	}
	static ParameterSet &getnamedparameterset(const std::string &name) {
		return getnamedinstance(name);
	}
	template<typename VAL>
	static VAL getpar(const std::string &key, VAL *val) {
		return getparameterset().getval(key, val);
	}
	template<typename VAL>
	static VAL getnamedsetpar(const std::string &name, const std::string &key,
			VAL *val) {
		return getnamedparameterset(name).getval(key, val);
	}
	static void initdefaultkeyvals(std::string controlname,const std::map<std::string,double> &doublepars,
			const std::map<std::string,std::string> &stringpars,
			const std::map<std::string,int> &intpars,
			const std::map<std::string, std::vector<double>> &doublevecpars,
			const std::map<std::string, std::vector<std::string>> &stringvecpars,
			const std::map<std::string, std::vector<int>> &intvecpars){
		ParameterSet &pset=getnamedparameterset(controlname);
		pset.InitDoubleKeys(NSPdstl::getkeyvec(doublepars));
		pset.InitIntKeys(NSPdstl::getkeyvec(intpars));
		pset.InitStringKeys(NSPdstl::getkeyvec(stringpars));
		pset.InitDoubleVecKeys(NSPdstl::getkeyvec(doublevecpars));
		pset.InitStringVecKeys(NSPdstl::getkeyvec(stringvecpars));
		pset.InitIntVecKeys(NSPdstl::getkeyvec(intvecpars));
		pset.initialized=true;
		for (auto &kv : doublepars) {
				pset.DoublePar.insert(kv.first, kv.second);
		}
		for (auto &kv : intpars) {
				pset.IntPar.insert(kv.first, kv.second);
		}
		for (auto &kv : stringpars) {
				pset.StringPar.insert(kv.first, kv.second);
		}
		for (auto &kv : doublevecpars) {
				pset.DoubleVecPar.insert(kv.first, kv.second);
		}
		for (auto &kv : intvecpars) {
				pset.IntVecPar.insert(kv.first, kv.second);
		}
		for (auto &kv : stringvecpars) {
				pset.StringVecPar.insert(kv.first, kv.second);
		}
	}
	static int adjustvalues(std::string controlname,const std::vector<std::string> &controllines) {
		ParameterSet & pset = getparameterset(controlname);
		assert(pset.initialized);
		int nsuccess = 0;
		for (auto &line : controllines)
			if (pset.readline(line))
				nsuccess++;
			else
				std::cout<<"Unsupported key-value found: " <<line<<std::endl;
		return nsuccess;
	}

private:
	TypedMultiInstanceControls() {;}
};
}

#endif /* PARAMETERS_H_ */
