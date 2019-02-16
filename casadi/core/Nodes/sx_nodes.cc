
namespace casadi
{
  
PlusNode(const Symbolic &a, const Symbolic &b)
{
}
int PlusNode::eval(const double** arg, double** res, casadi_int* iw, double* w) const
{
  if (x.is_zero())
    return y;
  else if (y->is_zero()) // term2 is zero
    return x;
  else if (y.is_op(OP_NEG))  // x + (-y) -> x - y
    return x - (-y);
  else if (x.is_op(OP_NEG)) // (-x) + y -> y - x
    return y - x.dep();
  else if (x.is_op(OP_MUL) && y.is_op(OP_MUL) &&
           x.dep(0).is_constant() && static_cast<double>(x.dep(0))==0.5 &&
           y.dep(0).is_constant() && static_cast<double>(y.dep(0))==0.5 &&
           is_equal(y.dep(1), x.dep(1), SXNode::eq_depth_)) // 0.5x+0.5x = x
    return x.dep(1);
  else if (x.is_op(OP_DIV) && y.is_op(OP_DIV) &&
           x.dep(1).is_constant() && static_cast<double>(x.dep(1))==2 &&
           y.dep(1).is_constant() && static_cast<double>(y.dep(1))==2 &&
           is_equal(y.dep(0), x.dep(0), SXNode::eq_depth_)) // x/2+x/2 = x
    return x.dep(0);
  else if (x.is_op(OP_SUB) && is_equal(x.dep(1), y, SXNode::eq_depth_))
    return x.dep(0);
  else if (y.is_op(OP_SUB) && is_equal(x, y.dep(1), SXNode::eq_depth_))
    return y.dep(0);
  else if (x.is_op(OP_SQ) && y.is_op(OP_SQ) &&
           ((x.dep().is_op(OP_SIN) && y.dep().is_op(OP_COS))
            || (x.dep().is_op(OP_COS) && y.dep().is_op(OP_SIN)))
           && is_equal(x.dep().dep(), y.dep().dep(), SXNode::eq_depth_))
    return 1; // sin^2 + cos^2 -> 1
  break;
}

} // namespace casadi
