// Copyright 2010, 2011 Hans Pirnay
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
  SmartPtr<TNLP> sens_tnlp = new ParametricTNLP();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  SmartPtr<IpoptNLP> ipopt_nlp = app_ipopt->IpoptNLPObject();
  SmartPtr<const Matrix> opt_jac_c_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_c_p(*app_ipopt->IpoptDataObject()->curr()->x());
  SmartPtr<const Matrix> opt_jac_d_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_d_p(*app_ipopt->IpoptDataObject()->curr()->x());
  SmartPrt<const Matrix> opt_h_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->h_p(*app_ipopt->IpoptDataObject()->curr()->x());
  opt_jac_c_p->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
  opt_jac_d_p->Print(*app_ipopt->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");

}
