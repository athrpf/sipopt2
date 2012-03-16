// Copyright 2010, 2011, 2012 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-01-05

#include "parametricTNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "SensApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "SensRegOp.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpDenseVector.hpp"


int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app_ipopt = new IpoptApplication();

  //SmartPtr<SensApplication> app_sens = new SensApplication(app_ipopt->Jnlst(),
  //app_ipopt->Options(),
  //app_ipopt->RegOptions());

  // Register sIPOPT options
  //RegisterOptions_sIPOPT(app_ipopt->RegOptions());
  //app_ipopt->Options()->SetRegisteredOptions(app_ipopt->RegOptions());

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app_ipopt->Initialize("");
  if (retval != Solve_Succeeded) {
    //printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }
  app_ipopt->Initialize();

  // create AmplSensTNLP from argc. This is an nlp because we are using our own TNLP Adapter
  SmartPtr<ParaTNLP> sens_tnlp = new ParametricTNLP();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  SmartPtr<const IteratesVector> curr = app_ipopt->IpoptDataObject()->curr();
  SmartPtr<const Vector> x = curr->x();
  SmartPtr<const Vector> y_c = curr->y_c();
  SmartPtr<const Vector> y_d = curr->y_d();

  // Get the parameter sensitivity matrices
  SmartPtr<IpoptNLP> ipopt_nlp = app_ipopt->IpoptNLPObject();
  SmartPtr<const Matrix> opt_jac_c_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_c_p(*x);
  SmartPtr<const Matrix> opt_jac_d_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_d_p(*x);
  SmartPtr<const Matrix> opt_h_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->h_p(*x, 1.0, *y_c, *y_d);
  opt_jac_c_p->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
  opt_jac_d_p->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");
  opt_h_p->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_h_p");

  // Set up the perturbed parameters
  SmartPtr<const DenseVector> p0 = dynamic_cast<const DenseVector*>(GetRawPtr(ipopt_nlp->p()));
  SmartPtr<DenseVector> dp = dynamic_cast<DenseVector*>(p0->MakeNewCopy());
  Number* dp_ptr = dp->Values();
  dp_ptr[0] = -0.5;
  dp_ptr[1] = 0.0;

  // Get the (factorized) KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app_ipopt->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver_ = pd_search->PDSolver();

  // Set up RHS and LHS for solve
  SmartPtr<IteratesVector> rhs = app_ipopt->IpoptDataObject()->curr()->MakeNewIteratesVector();
  rhs->Set(0.0);
  SmartPtr<Vector> rhs_x = x->MakeNew();
  opt_h_p->MultVector(1.0, *dp, 0.0, *rhs_x);
  rhs->Set_x_NonConst(*rhs_x);
  SmartPtr<Vector> rhs_c = y_c->MakeNew();
  opt_jac_c_p->MultVector(1.0, *dp, 0.0, *rhs_c);
  rhs->Set_y_c_NonConst(*rhs_c);
  SmartPtr<Vector> rhs_d = y_d->MakeNew();
  opt_jac_d_p->MultVector(1.0, *dp, 0.0, *rhs_d);
  rhs->Set_y_d_NonConst(*rhs_d);

  SmartPtr<IteratesVector> lhs = rhs->MakeNewIteratesVector();

  pd_solver_->Solve(-1.0, 0.0, *rhs, *lhs, false, false);

  lhs->Axpy(1.0, *curr);
  lhs->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "perturbed_x");
}
