/* 
 * File:   StringUtil.h
 * Author: Ha Minh Lam
 *
 * Created on 16 June 2011, 11:33
 */

#ifndef STRING_H
#define	STRING_H

#include <string>
#include <sstream>

using namespace std;

class String : string {
public:
    /**
     * Indicates if the given string is an integer or not.
     * @param str The string
     * @return <b>true</b> if the given string is an integer, <b>false</b> otherwise.
     */
    static const bool isInteger(string str);
    
    /**
     * Parse the given integer to string. If the number of digits in the given 
     * integer less than <b>minDigitNum</b>, some extra 0's will added to the 
     * left side of this integer.
     * @param integer
     * @param minDigitNum
     * @return The string that equivalent to the given integer.
     */
    static string formatInt(const int& integer, int minDigitNum);

    /**
     * Delete white spaces at the beginning and the end of the given string
     * @param inStr
     * @return A new string with all white spaces at the beginning and end removed.
     */
    static string trim(string inStr);
    
    /**
     * Delete all white spaces occur in the given string.
     * @param inStr
     * @return A new string that has no white space.
     */
    static string deleteAllSpace(string inStr);
    
private:
    String();
    String(const String& orig);
    virtual ~String();

};


#endif	/* STRING_H */

