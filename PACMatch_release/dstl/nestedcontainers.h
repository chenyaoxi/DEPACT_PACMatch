/*
 * nestedcontainers.h
 *
 *  Created on: 2015年10月23日
 *      Author: hyliu
 */

#ifndef NESTEDCONTAINERS_H_
#define NESTEDCONTAINERS_H_
#include <vector>
#include <iostream>
#include <functional>  //c++11 functional programming

/*!class templates that combine some STL containers
 *
 */
namespace NSPdstl {

/*! Wraps a stl container T of data D as a member.
 *
 * Classes derived form this template may somehow act like STL container-derived classes
 *  Implemented to avoid deriving from STL classes directly.
 *  T should be the actual container type, D should be the contained data type.
 *  For example: SLtWrapper<std::vector<double>,double>
 */
template<typename T, typename D> class StlWrapper
{
public:
    typedef typename T::iterator iterator;
    typedef typename T::const_iterator const_iterator;
    iterator begin() {
        return data_.begin();
    }
    iterator end() {
        return data_.end();
    }
    const_iterator cbegin() const {
        return data_.cbegin();
    }
    const_iterator cend() const {
        return data_.cend();
    }
    D & back() {
        return data_.back();
    }
    void push_back(const D & d) {
        data_.push_back(d);
    }
    void pop_back() {
        data_.pop_back();
    }
    typename T::size_type size() const {
        return data_.size();
    }
    void resize(int i) {
        data_.resize(i);
    }
    bool empty() const {
        return data_.empty();
    }

    /*! accessor for the wrapped container
     */
    T & data() {
        return data_;
    }

    /*! const accessor to the wrapped container
     *
     */
    const T & data() const {
        return data_;
    }
private:
    T data_;
};

/*! class template for a vector of vectors
 *
 */
template<typename T> class NestedVector: public StlWrapper<
                std::vector<std::vector<T> >, std::vector<T> >
{
public:
    /*! carry out operations based on each element of each vector
     *
     * @param f 		The functor to call for each element.
     * @param delim	The functor to call after finishing on  each vector.
     */
    void processElements(std::function<void(const T&)> f,
                    std::function<void()> delim =
                                    []() {std::cout <<std::endl;}) const {
        for (auto itout = this->cbegin();
                        itout != this->cend(); itout++) {
            for (auto itinner = itout->begin();
                            itinner != itout->end();
                            itinner++) {
                f(*itinner);
            }        //inner vector
            delim();
        }        //outter vector
    }        // processElements

    /*! print  the elements to std::cout. Elements separated by tab, vectors separated by neew lines.
     *
     */
    void printElements(char spacer = '\t') {
        processElements(
                        [spacer](const T& t) {std::cout <<t<<spacer;});
    }

    /*! random access operator not implemented in StlWrapper.
     *
     */
    std::vector<T> & operator[](int idx) {
        return nestedvector()[idx];
    }

    /*!random access operator, const version
     *
     */
    const std::vector<T> & operator[](int idx) const {
        return nestedvector()[idx];
    }

    /*!access the wrapped nested container
     *
     */
    std::vector<std::vector<T> > & nestedvector() {
        return this->data();
    }

    /*! access the wrapped nested container, const version
     *
     */
    const std::vector<std::vector<T> > & nestedvector() const {
        return this->data();
    }
};
//class NestedVector

/*! compute the average of a given column in a table
 *
 * @param table 	The input data
 * @param col     The column id
 * @return average and standard variation of all data in the column
 */
std::pair<double, double> columnaverage(
                const NestedVector<double> & table,
                int col);

}        //namespace stlextension

#endif /* NESTEDCONTAINERS_H_ */
