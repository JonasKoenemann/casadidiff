/*
 *    This file is part of CasADi.
 *
 *    CasADi -- A symbolic framework for dynamic optimization.
 *    Copyright (C) 2010 by Joel Andersson, Moritz Diehl, K.U.Leuven. All rights reserved.
 *
 *    CasADi is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    CasADi is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with CasADi; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "integrator_internal.hpp"
#include <cassert>
#include "../stl_vector_tools.hpp"
#include "jacobian.hpp"
#include "../matrix/matrix_tools.hpp"
#include "../mx/mx_tools.hpp"
#include "../sx/sx_tools.hpp"
#include "mx_function.hpp"
#include "sx_function.hpp"

INPUTSCHEME(IntegratorInput)
OUTPUTSCHEME(IntegratorOutput)

using namespace std;
namespace CasADi{

IntegratorInternal::IntegratorInternal(const FX& fd, const FX& fq) : fd_(fd), fq_(fq){
  new_design_ = false;
  ctorInit();
}

IntegratorInternal::IntegratorInternal(const FX& f, const FX& g, const FX& h) : f_(f), g_(g), h_(h){
  new_design_ = true;
  ctorInit();
}

void IntegratorInternal::ctorInit(){
  // set default options
  setOption("name","unnamed_integrator"); // name of the function 
  
  // Additional options
  addOption("print_stats",                 OT_BOOLEAN,  false, "Print out statistics after integration");
  addOption("nrhs",                        OT_INTEGER, 1); // number of right hand sides
  addOption("t0",                          OT_REAL, 0.0); // start of the integration
  addOption("tf",                          OT_REAL, 1.0); // end of the integration
  
  // Negative number of parameters for consistancy checking
  np_ = -1;
}

IntegratorInternal::~IntegratorInternal(){ 
}

void IntegratorInternal::setDimensions(int nx, int np){
  nx_ = nx;
  np_ = np;
  
  // Allocate space for inputs
  input_.resize(INTEGRATOR_NUM_IN);
  input(INTEGRATOR_X0)  = DMatrix(nx_,1,0); // initial state value
  input(INTEGRATOR_XP0) = DMatrix(nx_,1,0); // initial state derivative value
  input(INTEGRATOR_P)   = DMatrix(np_,1,0); // parameter
  
  // Allocate space for outputs
  output_.resize(INTEGRATOR_NUM_OUT);
  output(INTEGRATOR_XF) = DMatrix(nx_,1,0);
  output(INTEGRATOR_XPF)= DMatrix(nx_,1,0);
}

void IntegratorInternal::evaluate(int nfdir, int nadir){
  
  // Reset solver
  reset(nfdir, nadir);

  // Integrate forward to the end of the time horizon
  integrate(tf_);

  // If backwards integration is needed
  if(nadir>0){
    
    // Re-initialize backward problem
    resetAdj();

    // Integrate backwards to the beginning
    integrateAdj(t0_);
  }
  
  // Print statistics
  if(getOption("print_stats")) printStats(std::cout);
}

void IntegratorInternal::init(){
  
  // Initialize the functions and get dimensions
  if(new_design_){
    
    // Initialize the functions
    casadi_assert(!f_.isNull());

    // Initialize, get and assert dimensions of the forward integration
    if(!f_.isInit()) f_.init();
    nxd_ = f_.input(DAE_F_XD).numel();
    nxa_ = f_.input(DAE_F_XA).numel();
    np_  = f_.input(DAE_F_P).numel();
    nxq_ = f_.output(DAE_F_QUAD).numel();
    casadi_assert_message(f_.output(DAE_F_ODE).numel()==nxd_,"Inconsistent dimensions");
    casadi_assert_message(f_.output(DAE_F_ALG).numel()==nxa_,"Inconsistent dimensions");

    // Make sure that both h and g are given, or neither
    casadi_assert_message(h_.isNull()==g_.isNull(),"Either both h and g should be given, or neither of them");
    if(h_.isNull()){
      nyd_ = 0;
      nyq_ = 0;
      nya_ = 0;
    } else {
      // Initialize, get and assert dimensions of the terminal constraint function
      if(!h_.isInit()) h_.init();
      casadi_assert_message(h_.input(DAE_H_XD).numel()==nxd_,"Inconsistent dimensions");
      casadi_assert_message(h_.input(DAE_H_XA).numel()==nxa_,"Inconsistent dimensions");
      casadi_assert_message(h_.input(DAE_H_P).numel()==np_,"Inconsistent dimensions");
      nyd_ = h_.output(DAE_H_YD).numel();
      nyq_ = h_.output(DAE_H_YQ).numel();
      nya_ = h_.output(DAE_H_YA).numel();
      
      // Initialize and assert the dimensions of the backward integration
      if(!g_.isInit()) g_.init();
      casadi_assert_message(g_.input(DAE_G_XD).numel()==nxd_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(DAE_G_XA).numel()==nxa_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(DAE_G_YD).numel()==nyd_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(DAE_G_YA).numel()==nya_,"Inconsistent dimensions");
      casadi_assert_message(g_.input(DAE_G_P).numel()==np_,"Inconsistent dimensions");
      casadi_assert_message(g_.output(DAE_G_ODE).numel()==nyd_,"Inconsistent dimensions");
      casadi_assert_message(g_.output(DAE_G_QUAD).numel()==nyq_,"Inconsistent dimensions");
      casadi_assert_message(g_.output(DAE_G_ALG).numel()==nya_,"Inconsistent dimensions");
    }
    
    // Allocate space for inputs
    input_.resize(NEW_INTEGRATOR_NUM_IN);
    input(NEW_INTEGRATOR_XD0) = f_.output(DAE_F_ODE);
    input(NEW_INTEGRATOR_XQ0) = f_.output(DAE_F_QUAD);
    input(NEW_INTEGRATOR_XA0) = f_.output(DAE_F_ALG);
    input(NEW_INTEGRATOR_P) = f_.input(DAE_F_P);
  
    // Allocate space for outputs
    output_.resize(NEW_INTEGRATOR_NUM_OUT);
    output(NEW_INTEGRATOR_XDF) = input(NEW_INTEGRATOR_XD0);
    output(NEW_INTEGRATOR_XQF) = input(NEW_INTEGRATOR_XQ0);
    output(NEW_INTEGRATOR_XAF) = input(NEW_INTEGRATOR_XA0);
    if(!g_.isNull()){
      output(NEW_INTEGRATOR_YD0) = g_.output(DAE_G_ODE);
      output(NEW_INTEGRATOR_YQ0) = g_.output(DAE_G_QUAD);
      output(NEW_INTEGRATOR_YA0) = g_.output(DAE_G_ALG);
    }
  }
  
  // Make sure that the dimensions have been set
  casadi_assert_message(np_>=0, "\"setDimensions\" has not been called.");
  
  // Call the base class method
  FXInternal::init();

  // read options
  nrhs_ = getOption("nrhs");
  
  // Give an intial value for the time horizon
  t0_ = getOption("t0");
  tf_ = getOption("tf");
}

void IntegratorInternal::deepCopyMembers(std::map<SharedObjectNode*,SharedObject>& already_copied){
  FXInternal::deepCopyMembers(already_copied);
  if(new_design_){
    f_ = deepcopy(f_,already_copied);
    g_ = deepcopy(g_,already_copied);
    h_ = deepcopy(h_,already_copied);
  } else {
    fd_ = deepcopy(fd_,already_copied);
    fq_ = deepcopy(fq_,already_copied);
  }
}

} // namespace CasADi


