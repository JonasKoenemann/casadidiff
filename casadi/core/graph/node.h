/*
 *
 *    Copyright (C) 2019 Jonas Koenemann
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
 *    License along with this program;
 *    if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef CASADI_NODE_H_
#define CASADI_NODE_H_

namespace casadi {

class Node {

  // evaluate node numerically
  int MXNode::eval(const double** arg, double** res, casadi_int* iw, double* w) const;

  // evaluate node symbolically
  int MXNode::eval_sym(const Symbolic** arg, Symbolic** res, casadi_int* iw, Symbolic* w) const;

  // forwards differentiation
  void MXNode::ad_forward(const vector<vector<MX> >& fseed,
    vector<vector<MX> >& fsens) const;

  // backwards differentiation
  void MXNode::ad_reverse(const vector<vector<MX> >& aseed,
    vector<vector<MX> >& asens) const;

  // forwards sparsity pattern propagation
  int MXNode::sp_forward(const bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

  // backwards sparsity pattern
  int MXNode::sp_reverse(bvec_t** arg, bvec_t** res, casadi_int* iw, bvec_t* w) const;

  // code generation
  virtual void generate(CodeGenerator& g,
    const std::vector<casadi_int>& arg,
    const std::vector<casadi_int>& res) const;

  void display(std::ostream& stream);

}
} // namespace casadi


#endif // CASADI_NODE_H_
