/*
 *    This program is a derivative work of CasADi.
 *    The original program has been altered starting from February 15, 2019.
 *    The license of this file was changed from LGPL to GPL on February 15, 2019.
 *
 *    Copyright (C) 2019 Jonas Koenemann
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
 *                            K.U. Leuven. All rights reserved.
 *    Copyright (C) 2011-2014 Greg Horn
 *
 *    This program is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License and GNU Lesser General Public License along with this program;
 *    if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

 // history:
 //   JK 2019/02/16: remove implementations to convert to interface


#ifndef CASADI_GENERIC_EXPRESSION_HPP
#define CASADI_GENERIC_EXPRESSION_HPP

#include "calculus.hpp"

namespace casadi {


  /** \brief Expression interface
  *
  This is a common base class for SX, MX and Matrix<>, introducing a uniform syntax and implementing
  common functionality using the curiously recurring template pattern (CRTP) idiom.\n

  \author Joel Andersson
  \date 2012
*/
template<typename ExType>
class GenericExpression {
  public:

  /** \brief Addition: (x,y) -> x + y */
  ExType plus(const ExType &x, const ExType &y);
  ExType operator+(const ExType &x, const ExType &y);

  /** \brief Subtraction: (x,y) -> x - y */
  ExType minus(const ExType &x, const ExType &y);
  ExType& operator-=(const ExType &y);

  /** \brief Elementwise multiplication: (x,y) -> x .* y */
  ExType times(const ExType &x, const ExType &y);
  ExType& operator*=(const ExType &y);

  /** \brief Elementwise division: (x,y) -> x ./ y */
  ExType rdivide(const ExType &x, const ExType &y);
  ExType& operator/=(const ExType &y);

  /** \brief Logical less than: (x,y) -> x < y */
  ExType lt(const ExType &x, const ExType &y);
  ExType operator<(const ExType &x, const ExType &y);

  /** \brief Logical less or equal to: (x,y) -> x <= y */
  ExType le(const ExType &x, const ExType &y);
  ExType operator<=(const ExType &x, const ExType &y);

  /** \brief Logical greater than: (x,y) -> x > y */
  ExType gt(const ExType &x, const ExType &y);
  ExType operator>(const ExType &x, const ExType &y);

  /** \brief Logical greater or equal to: (x,y) -> x <= y */
  ExType ge(const ExType &x, const ExType &y);
  ExType operator>=(const ExType &x, const ExType &y);

  /** \brief Logical equal to: (x,y) -> x == y */
  ExType eq(const ExType &x, const ExType &y);
  ExType operator==(const ExType &x, const ExType &y);

  /** \brief Logical not equal to: (x,y) -> x != y */
  ExType ne(const ExType &x, const ExType &y);
  inline ExType operator!=(const ExType &x, const ExType &y);

  /** \brief Logical `and`
   * Returns (an expression evaluating to) 1 if both
   * expressions are nonzero and 0 otherwise
   */
  ExType logic_and(const ExType &x, const ExType &y);
  ExType operator&&(const ExType &x, const ExType &y);

  /** \brief  Logical `or`
   * returns (an expression evaluating to) 1 if at
   * least one expression is nonzero and 0 otherwise
   */
  ExType logic_or(const ExType &x, const ExType &y);
  ExType operator||(const ExType &x, const ExType &y);

   /** \brief  Logical `not` x -> !x
    * Returns (an expression evaluating to) 1 if
    * expression is zero and 0 otherwise
    */
  ExType logic_not(const ExType& x);
  ExType operator!() const;

  /** \brief Absolute value: x -> abs(x) */
  ExType abs(const ExType& x);
  ExType fabs(const ExType& x);

  /** \brief Square root: x -> sqrt(x) */
  ExType sqrt(const ExType& x);

  /** \brief Square: x -> x^2 */
  ExType sq(const ExType& x);

  /** \brief Sine: x -> sin(x) */
  ExType sin(const ExType& x);

  /** \brief Cosine: x -> cos(x) */
  ExType cos(const ExType& x);

  /** \brief Tangent: x -> tan(x) */
  ExType tan(const ExType& x);

  /** \brief Arc tangent: x -> atan(x) */
  ExType atan(const ExType& x);

  /** \brief Arc sine: x -> asin(x) */
  ExType asin(const ExType& x);

  /** \brief Arc cosine: x -> acos(x) */
  ExType acos(const ExType& x);

  /** \brief Hyperbolic tangent: x -> tanh(x) */
  ExType tanh(const ExType& x);

  /** \brief Hyperbolic sin: x -> sinh(x) */
  ExType sinh(const ExType& x);

  /** \brief Hyperbolic cosine: x -> cosh(x) */
  ExType cosh(const ExType& x);

  /** \brief Inverse hyperbolic tangent: x -> atanh(x) */
  ExType atanh(const ExType& x);

  /** \brief Inverse hyperbolic sin: x -> asinh(x) */
  ExType asinh(const ExType& x);

  /** \brief Inverse hyperbolic cosine: x -> acosh(x) */
  ExType acosh(const ExType& x);

  /** \brief Elementwise exponential: x -> exp(x) */
  ExType exp(const ExType& x);;

  /** \brief Natural logarithm: x -> log(x) */
  ExType log(const ExType& x);;

  /** \brief Base-10 logarithm: x -> log10(x) */
  ExType log10(const ExType& x);

  /** \brief Round down to nearest integer: x -> floor(x) */
  ExType floor(const ExType& x);

  /** \brief Round up to nearest integer: x -> ceil(x) */
  ExType ceil(const ExType& x);

  /** \brief Error function: x -> erf(x) */
  ExType erf(const ExType& x);

  /** \brief Inverse error function: x -> erfinv(x) */
  ExType erfinv(const ExType& x);

  /** \brief Sign function:
      sign(x)   := -1 for x<0
      sign(x)   :=  1 for x>0,
      sign(0)   :=  0
      sign(NaN) :=  NaN
  */
  ExType sign(const ExType& x);

  /** \brief Elementwise power: (x,y) -> x.^y */
  ExType pow(const ExType& x, const ExType& y);

  /** \brief Remainder after division: (x,y) -> mod(x,y) */
  ExType mod(const ExType& x, const ExType& y);

  friend inline ExType fmod(const ExType& x, const ExType& y);

  /** \brief Two argument arc tangent: (x,y) -> atan2(x,y) */
  ExType atan2(const ExType& x, const ExType& y);

  /** \brief Conditional assignment: (x,y) -> x ? y : 0 */
  ExType if_else_zero(const ExType& x, const ExType& y);

  /** \brief Smallest of two values: (x,y) -> min(x,y) */
  ExType fmin(const ExType& x, const ExType& y);

  /** \brief Largest of two values: (x,y) -> max(x,y) */
  ExType fmax(const ExType& x, const ExType& y);

  /** \brief Check if two nodes are equivalent up to a given depth.
   * Depth=0 checks if the expressions are identical, i.e. points to the same node.
   *
   * a = x*x
   * b = x*x
   *
   *  is_equal(a,b,0)  will return false, but a.is_equal(a,b,1) will return true
   */
   friend inline bool is_equal(const ExType& x, const ExType& y, casadi_int depth=0);

   /// Copy sign
   ExType copysign(const ExType& x, const ExType& y);

   /// Elementwise power with const power
   ExType constpow(const ExType& x, const ExType& y);

   /// Debug printing
   ExType printme(const ExType& x, const ExType& y);

/** @} */
};

} // namespace casadi

#endif // CASADI_GENERIC_EXPRESSION_HPP
