/*
 * inputlines.h
 *
 *  Created on: 2015年9月28日
 *      Author: hyliu
 */

#ifndef INCLUDE_INPUTLINES_H_
#define INCLUDE_INPUTLINES_H_

#include <string>
#include <vector>
#include <memory>
#include <boost/lexical_cast.hpp>
#include "macros.h"  //for disable copy and assign
#include "dstl/nestedcontainers.h"

namespace NSPdataio {

/*! Echos the error message to std::out and exits the program.
 *
 * This error handler could be replaced by more sophisticated exception mechanisms in the future.
 */
void attendError(const std::string &errormessage);
std::pair<std::string,std::string> keyvalstrings(const std::string &line);
/*!parse one line into words
 *
 *Using white space separators  if cols empty. Otherwise cols gives the set of column widths.
 */
std::vector<std::string> parseline(
                const std::string & inputline,
                const std::vector<int> & cols);

/*! Convert a integer-selecting string to a list of integers. For example, the string
 * "5, 6,-8, 11-13" would to converted to the integers {5,6,7,8,9,11,12,13}
 *
 */
std::vector<int> stringtoselection(const std::string &str);

/*! Convert the words in the input string vector into a vector of numerical values.
 *
 *  Words that cannot be converted are ignored. Needs boost library.
 * @param	line	Input vector containing words to be converted.
 * @return  		A vector contain wanted numerical values.
 */
template<typename T>
std::vector<T> wordsToNumbers(
                const std::vector<std::string> & line) {
    std::vector<T> numbers;
    for (auto it = line.begin(); it != line.end(); it++) {
        T d;
        try {
            d = boost::lexical_cast<T>(*it);
        } catch (const boost::bad_lexical_cast &e) {
            continue;
        }
        numbers.push_back(d);
    }
    return numbers;
}

/*! Functors to transform a NestedVector<T> into NestedVector<D>.
 *
 */
template<typename T, typename D>
class convertelements
{
public:

    /*! initialize with an input pointer to destination NestedVector<D>
     *
     */
    convertelements(NSPdstl::NestedVector<D> *pt) :
                    table_(pt), currentline_() {
    }

    /*! convert one element t to d and store  it to destination
     *
     * If conversion failed, the input element is simply ignored.
     * @param t  The element to be converted
     */
    void convertOneEle(const T & t) {
        D d;
        try {
            d = boost::lexical_cast<D>(t);
            currentline_.push_back(d);
        } catch (const boost::bad_lexical_cast &e) {
            ;
        }
    }

    /*!Define the operation at the end of an vector.
     *
     * Push the non-empty converted vector into destination. Clear the buffer.
     *
     */
    void lineEnd() {
        if (currentline_.size() != 0)
            table_->push_back(std::move(currentline_));
        currentline_.clear();
    }
private:
    /*!Pointer to destination NestedVector
     */
    NSPdstl::NestedVector<D> *table_;

    /*!buffer to store one converted vector
     */
    std::vector<D> currentline_;
};

/*!All lines read from an input file. Each line as a string
 *
 */
class TextLines: public NSPdstl::StlWrapper<
                std::vector<std::string>, std::string>
{
public:
    TextLines() {
        ;
    }
    /*!read in content from the input
     *
     * @param filename The inputfilename
     */
    void init(const std::string &filename);

    /*! access the stored lines as a vector of strings
     *
     */
    std::vector<std::string> & lines() {
        return NSPdstl::StlWrapper<std::vector<std::string>,
                        std::string>::data();
    }

    /*! access stored lines, constant version
     *
     */
    const std::vector<std::string> & lines() const {
        return NSPdstl::StlWrapper<std::vector<std::string>,
                        std::string>::data();
    }
};
/*! Parsed text input, each line into a vector of input words.
 *
 * Provide converter into numeric values
 */
class InputLines: public NSPdstl::NestedVector<
                std::string>
{
public:
    /*! Provide default constructor
     *
     */
    InputLines() {
        ;
    }

    /*! read and parse lines from file
     *
     *  The input file are scanned line by line. Each line is parsed into
     *  a vector of words, using white spaces as delimiters, or
     *  	use based on a given set of column widths;
     *  . Each entry in the
     *  InputLines vector contains one parsed line.
     *   Lines start with the special comment char are ignored.
     *   If any error happens handling the input, an error handler is called.
     * @param	filename	The input file name.
     * @param  commentstarter			The special character to start a comment line.
     * @param  cols  	If given, the line will be segmented with fixed column widths.
     */
    void init(const std::string &filename,
                    char commentstarter,
                    const std::vector<int> &cols =
                                    std::vector<int>());

    /*! extract the numerical values in the parsed lines.
     *
     * @return The extracted values from one line are stored in one vector of doubles.
     * 			All vectors are store in the returned vector of vectors.
     */
    NSPdstl::NestedVector<double> extractTable() const;

    /*! extracts lines that contain a label word.
     *
     * @param label 	The label word.
     * @return 		The extracted lines.
     */
    std::unique_ptr<InputLines> extractLines(
                    const std::string & label) const;

    /*! removes lines that contain a label word.
     *
     * @param label	 lines containing the label word will be removed.
     */
    void removeLines(const std::string & label);
private:
    DISALLOW_COPY_AND_ASSIGN(InputLines);
};
//class InputLines

}//namespace dataio

#endif /* INCLUDE_INPUTLINES_H_ */
