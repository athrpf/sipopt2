// Copyright (C) 2012 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Hans Pirnay    2012-03-14

#ifndef __IPParaTNLPWrapper_HPP__
#define __IPParaTNLP_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
#include "IpException.hpp"
#include "IpAlgTypes.hpp"
#include "IpReturnCodes.hpp"
#include "IpParaTNLP.hpp"

#include <map>

namespace Ipopt
{
  // forward declarations
  class IpoptData;
  class IpoptCalculatedQuantities;
  class IteratesVector;

  /** Base class for all NLP's that use standard triplet matrix form
   *  and dense vectors.  This is the standard base class for all
   *  NLP's that use the standard triplet matrix form (as for Harwell
   *  routines) and dense vectors. The class TNLPAdapter then converts
   *  this interface to an interface that can be used directly by
   *  ipopt.
   *
   *  This interface presents the problem form:
   *
   *     min f(x)
   *
   *     s.t. gL <= g(x) <= gU
   *
   *          xL <=  x   <= xU
   *
   *  In order to specify an equality constraint, set gL_i = gU_i =
   *  rhs.  The value that indicates "infinity" for the bounds
   *  (i.e. the variable or constraint has no lower bound (-infinity)
   *  or upper bound (+infinity)) is set through the option
   *  nlp_lower_bound_inf and nlp_upper_bound_inf.  To indicate that a
   *  variable has no upper or lower bound, set the bound to
   *  -ipopt_inf or +ipopt_inf respectively
   */
  class ParaTNLPWrapper : public ParaTNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    ParaTNLPWrapper(SmartPtr<TNLP> tnlp) : tnlp_(tnlp)
    {}

    /** Default destructor */
    virtual ~ParaTNLPWrapper()
    {}
    //@}

    virtual bool get_nlp_info(Index& n, Index& np, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, Index& nnz_jac_g_p,
			      Index& nnz_h_lag_p, TNLP::IndexStyleEnum& index_style)
    {
      np = 0;
      nnz_jac_g_p = 0;
      nnz_h_lag_p = 0;
      return tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, index_style);
    }

    virtual bool get_var_con_metadata(Index n,
                                      StringMetaDataMapType& var_string_md,
                                      IntegerMetaDataMapType& var_integer_md,
                                      NumericMetaDataMapType& var_numeric_md,
                                      Index m,
                                      StringMetaDataMapType& con_string_md,
                                      IntegerMetaDataMapType& con_integer_md,
                                      NumericMetaDataMapType& con_numeric_md)

    {
      return tnlp_->get_var_con_metadata(n, var_string_md,
					 var_integer_md,
					 var_numeric_md,
					 m,
					 con_string_md,
					 con_integer_md,
					 con_numeric_md);
    }

    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u)
    {
      return get_bounds_info(n, x_l, x_u, m, g_l, g_u);
    }

    virtual bool get_parameters(Index np, Number* p)
    {
      return true;
    }

    virtual bool get_scaling_parameters(Number& obj_scaling,
                                        bool& use_x_scaling, Index n,
                                        Number* x_scaling,
                                        bool& use_g_scaling, Index m,
                                        Number* g_scaling)
    {
      return tnlp_->get_scaling_parameters(obj_scaling, use_x_scaling, n,
					   x_scaling, use_g_scaling, m, g_scaling);
    }

    virtual bool get_variables_linearity(Index n, TNLP::LinearityType* var_types)
    {
      return tnlp_->get_variables_linearity(n, var_types);
    }

    virtual bool get_constraints_linearity(Index m, TNLP::LinearityType* const_types)
    {
      return tnlp_->get_constraints_linearity(m, const_types);
    }

    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda,
                                    Number* lambda)
    {
      return tnlp_->get_starting_point(n, init_x, x, init_z, z_L, z_U,
				       m, init_lambda, lambda);
    }

    virtual bool get_warm_start_iterate(IteratesVector& warm_start_iterate)
    {
      return tnlp_->get_warm_start_iterate(warm_start_iterate);
    }

    virtual bool eval_f(Index n, const Number* x, bool new_x,
			Index np, const Number* p, bool new_p,
                        Number& obj_value)
    {
      return tnlp_->eval_f(n, x, new_x, obj_value);
    }

    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
			     Index np, const Number* p, bool new_p,
			     Number* grad_f)
    {
      return tnlp_->eval_grad_f(n, x, new_x, grad_f);
    }

    virtual bool eval_g(Index n, const Number* x, bool new_x,
			Index np, const Number* p, bool new_p,
                        Index m, Number* g)
    {
      return tnlp_->eval_g(n, x, new_x, m, g);
    }

    virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
			    Index np, const Number* p, bool new_p,
                            Index m, Index nele_jac, Index* iRow,
                            Index *jCol, Number* values)
    {
      return tnlp_->eval_jac_g(n, x, new_x, m, nele_jac, iRow, jCol, values);
    }

    virtual bool eval_jac_g_p(Index n, const Number* x, bool new_x,
			      Index np, const Number* p, bool new_p,
			      Index m, Index nele_jac, Index* iRow,
			      Index *jCol, Number* values)
    {
      return true;
    }

    virtual bool eval_h(Index n, const Number* x, bool new_x,
			Index np, const Number* p, bool new_p,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess,
                        Index* iRow, Index* jCol, Number* values)
    {
      return tnlp_->eval_h(n, x, new_x, obj_factor, m, lambda, new_lambda, nele_hess,
			   iRow, jCol, values);
    }

    virtual bool eval_L_xp(Index n, const Number* x, bool new_x,
			   Index np, const Number* p, bool new_p,
			   Number obj_factor, Index m,
			   const Number* lambda, bool new_lambda,
			   Index nele_hess_p, Index* iRow, Index* jCol,
			   Number* values)
    {
      return true;
    }
    //@}

    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value,
                                   const IpoptData* ip_data,
                                   IpoptCalculatedQuantities* ip_cq)
    {
      return tnlp_->finalize_solution(status, n, x, z_L, z_U, m, g, lambda, obj_value, ip_data, ip_cq);
    }

    virtual void finalize_metadata(Index n,
                                   const StringMetaDataMapType& var_string_md,
                                   const IntegerMetaDataMapType& var_integer_md,
                                   const NumericMetaDataMapType& var_numeric_md,
                                   Index m,
                                   const StringMetaDataMapType& con_string_md,
                                   const IntegerMetaDataMapType& con_integer_md,
                                   const NumericMetaDataMapType& con_numeric_md)
    {
      return tnlp_->finalize_metadata(n, var_string_md,
				      var_integer_md,
				      var_numeric_md,
				      m,
				      con_string_md,
				      con_integer_md,
				      con_numeric_md);
    }


    /** Intermediate Callback method for the user.  Providing dummy
     *  default implementation.  For details see IntermediateCallBack
     *  in IpNLP.hpp. */
    virtual bool intermediate_callback(AlgorithmMode mode,
                                       Index iter, Number obj_value,
                                       Number inf_pr, Number inf_du,
                                       Number mu, Number d_norm,
                                       Number regularization_size,
                                       Number alpha_du, Number alpha_pr,
                                       Index ls_trials,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
    {
      return tnlp_->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du,
					  mu, d_norm, regularization_size, alpha_du,
					  alpha_pr, ls_trials, ip_data, ip_cq);
    }
    //@}

    /** @name Methods for quasi-Newton approximation.  If the second
     *  derivatives are approximated by Ipopt, it is better to do this
     *  only in the space of nonlinear variables.  The following
     *  methods are call by Ipopt if the quasi-Newton approximation is
     *  selected.  If -1 is returned as number of nonlinear variables,
     *  Ipopt assumes that all variables are nonlinear.  Otherwise, it
     *  calls get_list_of_nonlinear_variables with an array into which
     *  the indices of the nonlinear variables should be written - the
     *  array has the lengths num_nonlin_vars, which is identical with
     *  the return value of get_number_of_nonlinear_variables().  It
     *  is assumed that the indices are counted starting with 1 in the
     *  FORTRAN_STYLE, and 0 for the C_STYLE. */
    //@{
    virtual Index get_number_of_nonlinear_variables()
    {
      return tnlp_->get_number_of_nonlinear_variables();
    }

    virtual bool get_list_of_nonlinear_variables(Index num_nonlin_vars,
        Index* pos_nonlin_vars)
    {
      return tnlp_->get_list_of_nonlinear_variables(num_nonlin_vars, pos_nonlin_vars);
    }
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    //ParaTNLPWrapper();

    /** Copy Constructor */
    ParaTNLPWrapper(const ParaTNLPWrapper&);

    /** Overloaded Equals Operator */
    void operator=(const ParaTNLPWrapper&);
    //@}

    SmartPtr<TNLP> tnlp_;
  };

} // namespace Ipopt

#endif
