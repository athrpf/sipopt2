// Copyright 2010, 2011, 2012 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-30-04

#include <assert.h>
#include <stdio.h>

#include "parametricTNLP.hpp"
#include "IpDenseVector.hpp"
#include "IpIpoptData.hpp"

using namespace Ipopt;


/* Constructor */
ParametricTNLP::ParametricTNLP() :
  nominal_eta1_(4.5),
  nominal_eta2_(1.0),
  eta_1_perturbed_value_(4.0),
  eta_2_perturbed_value_(1.0)
{
}

ParametricTNLP::~ParametricTNLP()
{
}

bool ParametricTNLP::get_nlp_info(Index& n, Index& np, Index& m, Index& nnz_jac_g,
				  Index& nnz_h_lag, Index& nnz_jac_g_p,
				  Index& nnz_h_lag_p, TNLP::IndexStyleEnum& index_style)
{
  // x1, x2, x3
  n = 3;

  // eta1, eta2
  np = 2;

  // 2 constraints
  m = 2;

  nnz_jac_g = 6;

  nnz_h_lag = 3;

  nnz_jac_g_p = 2;
  nnz_h_lag_p = 1;

  index_style = TNLP::FORTRAN_STYLE;

  return true;
}

bool ParametricTNLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
				     Index m, Number* g_l, Number* g_u)
{
  for (Index k=0; k<3; ++k) {
    x_l[k] = 0.0;
    x_u[k] = 1.0e19;
  }

  g_l[0] = 0.0;
  g_u[0] = 0.0;
  g_l[1] = 0.0;
  g_u[1] = 0.0;

  return true;
}

bool ParametricTNLP::get_starting_point(Index n, bool init_x, Number* x,
					Index np, bool init_p, Number* p,
					bool init_z, Number* z_L, Number* z_U,
					Index m, bool init_lambda,
					Number* lambda)
{
  if (init_x) {
    for (Index k=0; k<n; ++k) {
      x[k] = 5.0;
    }
  }
  if (init_p) {
    p[0] = 5.;
    p[1] = 1.0;
  }

  return true;
}
/*
bool ParametricTNLP::get_parameters(Index np, Number* p)
{
  assert(np ==2);
  p[0] = 5.0;
  p[1] = 1.0;
  return true;
  }*/

bool ParametricTNLP::eval_f(Index n, const Number* x, bool new_x,
			    Index np, const Number* p, bool new_p,
			    Number& obj_value)
{
  obj_value = 0;
  for (Index k=0; k<n; ++k) {
    obj_value += x[k]*x[k];
  }
  return true;
}

bool ParametricTNLP::eval_grad_f(Index n, const Number* x, bool new_x,
				 Index np, const Number* p, bool new_p,
				 Number* grad_f)
{
  grad_f[0] = 2*x[0];
  grad_f[1] = 2*x[1];
  grad_f[2] = 2*x[2];
  return true;
}

bool ParametricTNLP::eval_g(Index n, const Number* x, bool new_x,
			    Index np, const Number* p, bool new_p,
			    Index m, Number* g)
{
  Number x1, x2, x3, eta1, eta2;
  x1 = x[0];
  x2 = x[1];
  x3 = x[2];
  eta1 = p[0];
  eta2 = p[1];
  g[0] = 6*x1+3*x2+2*x3-eta1;
  g[1] = eta2*x1+x2-x3-1;
  return true;
}

bool ParametricTNLP::eval_jac_g(Index n, const Number* x, bool new_x,
				Index np, const Number* p, bool new_p,
				Index m, Index nele_jac, Index* iRow, Index *jCol,
				Number* values)
{
  assert(n==3);
  assert(np==2);
  assert(m==2);
  if (values == NULL) {
    iRow[0] = 1; // dg1/dx1
    jCol[0] = 1;
    iRow[1] = 1; // dg1/dx2
    jCol[1] = 2;
    iRow[2] = 1; // dg1/dx3
    jCol[2] = 3;
    iRow[3] = 2; // dg2/dx1
    jCol[3] = 1;
    iRow[4] = 2; // dg2/dx2
    jCol[4] = 2;
    iRow[5] = 2; // dg2/dx3
    jCol[5] = 3;
  }
  else {
    values[0] = 6.0;
    values[1] = 3.0;
    values[2] = 2.0;
    values[3] = p[1];
    values[4] = 1.0;
    values[5] = -1.0;
  }
  return true;
}

bool ParametricTNLP::eval_jac_g_p(Index n, const Number* x, bool new_x,
				  Index np, const Number* p, bool new_p,
				  Index m, Index nele_jac_p, Index* iRow, Index *jCol,
				  Number* values)
{
  if (values == NULL) {
    iRow[0] = 1;
    jCol[0] = 1;
    iRow[1] = 2;
    jCol[1] = 2;
  } else {
    values[0] = -1.0;
    values[1] = x[0];
  }
  return true;
}

bool ParametricTNLP::eval_L_xp(Index n, const Number* x, bool new_x,
			       Index np, const Number* p, bool new_p,
			       Number obj_factor, Index m,
			       const Number* lambda, bool new_lambda,
			       Index nele_hess_p, Index* iRow, Index* jCol,
			       Number* values)
{
  std::cout << "nele_hess_p: " << nele_hess_p << std::endl;
  assert(nele_hess_p == 1);
  if (values==NULL) { // return structure
    iRow[0] = 1;
    jCol[0] = 2;
  } else {
    values[0] = lambda[1];
  }
  return true;
}

bool ParametricTNLP::eval_h(Index n, const Number* x, bool new_x,
			    Index np, const Number* p, bool new_p,
			    Number obj_factor, Index m, const Number* lambda,
			    bool new_lambda, Index nele_hess, Index* iRow,
			    Index* jCol, Number* values)
{
  if (values == NULL) {
    iRow[0] = 1;
    jCol[0] = 1;

    iRow[1] = 2;
    jCol[1] = 2;

    iRow[2] = 3;
    jCol[2] = 3;
  }
  else {
    values[0] = 2.0*obj_factor;
    values[1] = 2.0*obj_factor;
    values[2] = 2.0*obj_factor;
  }
  return true;
}

void ParametricTNLP::finalize_solution(SolverReturn status,
				       Index n, const Number* x, const Number* z_L, const Number* z_U,
				       Index m, const Number* g, const Number* lambda,
				       Number obj_value,
				       const IpoptData* ip_data,
				       IpoptCalculatedQuantities* ip_cq)
{
  // Check whether sIPOPT Algorithm aborted internally
  //  bool sens_internal_abort;
  //options_->GetBoolValue("sens_internal_abort", sens_internal_abort, "");

  // Get access to the metadata, where the solutions are stored. The metadata is part of the DenseVectorSpace.
  //SmartPtr<const DenseVectorSpace> x_owner_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(ip_data->curr()->x()->OwnerSpace()));

  //if (!IsValid(x_owner_space)) {
  //  printf("Error IsValid(x_owner_space) failed\n");
  //  return;
  //}
  //std::string state;
  //std::vector<Number> sens_sol_vec;
  //state = "sens_sol_state_1";
  //sens_sol_vec = x_owner_space->GetNumericMetaData(state.c_str());

  // Print the solution vector
  printf("\n"
	 "                Nominal                    Perturbed\n");
  for (Index k=0; k<n; ++k) {
    //printf("x[%3d]   % .23f   % .23f\n", k, x[k], sens_sol_vec[k]);
    printf("x[%3d]   % .23f \n", k, x[k]);
  }
}
