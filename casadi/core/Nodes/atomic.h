/*
 *    This program is a derivative work of CasADi.
 *    The original program has been altered starting from February 15, 2019.
 *    The license of this file was changed from LGPL to GPL on February 16, 2019.
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


#ifndef ADLIB_ATOMIC_H_
#define ADLIB_ATOMIC_H_


namespace adlib
{


typedef std::vector<Expression> ExpressionList;


class Expression {
  const &T eval_fcn(const Assignment &a);
}

class Matrix : Expression {
  ExpressionList fcn();
  ExpressionList der();
}

class DoubleExpression : double, Expression
{
  DoubleExpression(double v) {
    this = v;
  }
  const Expression& fcn() {
    return this;
  }

  const Expression& der() {
    //return DoubleExpression(1); return DoubleExpression(0); ?
  }
}

class SymPrimitive : Expression
{
  const Expression& eval_fcn(const Assignment& a) {
    return a.get(this);
  }
}

// Operation with multiple inputs and multiple outputs
class MimoOperation : Expression
{
  AtomicOperation(const Expression& inputs);

  virtual const Expression& fcn()=0;
  virtual const Expression& der()=0;

  const Expression& eval_fcn(const Assignment& a) {
    evaluatedList = {};
    for in : inputs {
      evaluatedList.append(in.eval_fcn(a));
    }
    return this->fcn(evaluatedList);
  }

  private const Expression& eval_der();
  private ExpressionList inputs;
}

// Operation with single input and a single output
class UnaryOperation : Expression
{
  const Expression& eval_fcn(const Assignment& a) {
    in = input.eval_fcn(a);
    return this->fcn(in);
  }
  private Expression input;
}

// Operation with two inputs and a asingle output
class BinaryOperation : Expression
{
  const Expression& eval_fcn(const Assignment& a) {
    in1 = input.eval_fcn(a);
    in2 = input.eval_fcn(a);
    return this->fcn(in1, in2);
  }
  private Expression input1;
  private Expression input2;
}

typedef Expression T;

// operator functions
const &T assign(const T &x) { return Assign(x); }
const &T add(const T &x, const T &y) { return Add(x,y); }
const &T sub(const T &x, const T &y) { return Sub(x,y); }


// unary operations
class Assign : UnaryOperation
{
  Assign(const T& x) { input = x; }
  static const T& fcn(const T& x) { return x; }
  static const T& der(const T& x, const T& f) { return {1}; }
}

// binary operations
class Add : BinaryOperation
{
  Add(const T& x, const T& y) { input1 = x; input2 = y; }
  static const T& fcn(const T& x, const T& y) { return  x+y; }
  static const T& der(const T& x, const T& y, const T& f) { {1,1}; }
};







// TODO: complete operations
class Sub : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x-y;}
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=1; d[1]=-1; }
};

/// Multiplication
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x*y;}
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=y; d[1]=x; }
};

/// Division
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x/y;}
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=1/y; d[1]=-f/y; }
};

/// Negation
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = -x; }
  void der(const T& x, const T& f, T* d) { d[0]=-1; }
};

/// Natural exponent
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = exp(x); }
  void der(const T& x, const T& f, T* d) { d[0]=f; }
};

/// Natural logarithm
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = log(x); }
  void der(const T& x, const T& f, T* d) { d[0]=1/x; }
};

/// Power, defined only for x>=0
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = pow(x, y); }
  // See issue #104 why d[0] is no longer y*f/x
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=y*pow(x, y-1); d[1]=log(x)*f; }
};

/// Power, defined only for y constant
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = pow(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=y*pow(x, y-1); d[1]=0; }
};

/// Square root
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = sqrt(x); }
  void der(const T& x, const T& f, T* d) { d[0]=1/(twice(f)); }
};

/// Square
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = sq(x); }
  void der(const T& x, const T& f, T* d) { d[0]=twice(x); }
};

/// Times two
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = 2.*x; }
  void der(const T& x, const T& f, T* d) { d[0] = 2; }
};

/// Sine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = sin(x); }
  void der(const T& x, const T& f, T* d) { d[0]=cos(x); }
};

/// Cosine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = cos(x);}
  void der(const T& x, const T& f, T* d) { d[0]=-sin(x); }
};

/// Tangent
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = tan(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 1/sq(cos(x)); }
};

/// Arcus sine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = asin(x); }
  void der(const T& x, const T& f, T* d) { d[0]=1/sqrt(1-x*x); }
};

/// Arcus cosine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = acos(x);}
  void der(const T& x, const T& f, T* d) { d[0]=-1/sqrt(1-x*x); }
};

/// Arcus tangent
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = atan(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 1/(1+x*x); }
};

/// Less than
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x < y;}
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Less or equal to
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x <= y; }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Floor function
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = floor(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 0; }
};

/// Ceil function
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = ceil(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 0; }
};

/// Remainder of division
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = fmod(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=1; d[1]=(f-x)/y; }
};

/// Equal to
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x==y; }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Not equal to
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x!=y;}
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Logical not
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = !x; }
  void der(const T& x, const T& f, T* d) { d[0] = 0; }
};

/// Logical and
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x && y; }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Logical or
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x || y; }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=d[1]=0; }
};

/// Error function
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = erf(x); }
  void der(const T& x, const T& f, T* d) { d[0] = (2/sqrt(pi))*exp(-x*x); }
};

/// Absolute value
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = fabs(x); }
  void der(const T& x, const T& f, T* d) { d[0]=sign(x); }
};

/// Sign
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = sign(x); }
  void der(const T& x, const T& f, T* d) { d[0]=0; }
};

/// Copysign
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = copysign(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { T e = 1; d[0]=copysign(e, y); d[1]=0; }
};

/// Minimum
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = fmin(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=x<=y; d[1]=!d[0]; }
};

/// Maximum
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = fmax(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=x>=y; d[1]=!d[0]; }
};

/// Elementwise inverse
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = 1./x; }
  void der(const T& x, const T& f, T* d) { d[0] = -f*f; }
};

/// Hyperbolic sine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = sinh(x); }
  void der(const T& x, const T& f, T* d) { d[0] = cosh(x); }
};

/// Hyperbolic cosine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = cosh(x); }
  void der(const T& x, const T& f, T* d) { d[0] = -sinh(x); }
};

/// Hyperbolic tangent
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = tanh(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 1-f*f; }
};

/// Inverse hyperbolic sine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = asinh(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 1/sqrt(1+x*x); }
};

/// Inverse hyperbolic cosine
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = acosh(x);}
  void der(const T& x, const T& f, T* d) { d[0] = 1/sqrt(x-1)/sqrt(x+1); }
};

/// Inverse hyperbolic tangent
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = atanh(x); }
  void der(const T& x, const T& f, T* d) { d[0] = 1/(1-x*x); }
};

/// Inverse of error function
class Assign : Atomic
{
  void fcn(const T& x, T& f) { f = erfinv(x); }
  void der(const T& x, const T& f, T* d) { d[0] = (sqrt(pi)/2)*exp(f*f); }
};

/// Identity operator with the side effect of printing
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) {f = printme(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=1; d[1]=0; }
};

/// Arctan2
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = atan2(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { T t = x*x+y*y; d[0]=y/t; d[1]=-x/t; }
};

/// Conditional assignment
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = if_else_zero(x, y); }
  void der(const T& x, const T& y, const T& f, T* d) { d[0]=0; d[1]=if_else_zero(x, T(1)); }
};

/// Inverse of error function
class Assign : Atomic
{
  void fcn(const T& x, const T& y, T& f) { f = x; }
  void der(const T& x, const T& y, const T& f, T* d) { d[0] = 1; d[1] = 0; }
};

} // namespace adlib
#endif // ADLIB_ATOMIC_H_
