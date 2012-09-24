// Copyright (C) 2004, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "AmplTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "IpoptConfig.h"
#ifdef HAVE_CSTRING
# include <cstring>
#else
# ifdef HAVE_STRING_H
#  include <string.h>
# else
#  error "don't have header file for string"
# endif
#endif

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

// for parametric stuff at the end
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"

using namespace Ipopt;
SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app);
SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app,
					                      SmartPtr<Matrix> sens_matrix);
// bewa01 starting to do stuff
bool doIntervallization(SmartPtr<IpoptApplication> app, SmartPtr<AmplSuffixHandler> suffix_handler,
			            SmartPtr<ParaTNLP> ampl_tnlp);

int main(int argc, char**args)
{

  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Check if executable is run only to print out options documentation
  if (argc == 2) {
    bool print_options = false;
    bool print_latex_options = false;
    if (!strcmp(args[1],"--print-options")) {
      print_options = true;
    }
    else if (!strcmp(args[1],"--print-latex-options")) {
      print_options = true;
      print_latex_options = true;
    }
    if (print_options) {
      SmartPtr<OptionsList> options = app->Options();
      options->SetStringValue("print_options_documentation", "yes");
      if (print_latex_options) {
        options->SetStringValue("print_options_latex_mode", "yes");
      }
      app->Initialize("");
      return 0;
    }
  }

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app->Initialize("");
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }

  // Add the suffix handler for scaling
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
  // Modified for warm-start from AMPL
  suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
  // Add the suffix for parameter-marking
  suffix_handler->AddAvailableSuffix("parameter", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("perturbed", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

  SmartPtr<ParaTNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()),
                                          app->Options(),
                                          args, suffix_handler);

  // Call Initialize again to process output related options
  retval = app->Initialize();
  if (retval != Solve_Succeeded) {
    printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
    exit(-101);
  }

  const int n_loops = 1; // make larger for profiling
  for (Index i=0; i<n_loops; i++) {
    retval = app->OptimizeTNLP(ampl_tnlp);
  }

  SmartPtr<Matrix> sens_matrix = getSensitivityMatrix(app);
  SmartPtr<Vector> delta_s = getDirectionalDerivative(app, sens_matrix);
  if (IsValid(delta_s))
    delta_s->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "delta_s");
  return 0;
}

SmartPtr<Matrix> getSensitivityMatrix(SmartPtr<IpoptApplication> app)
{
  // finalize_solution method in AmplTNLP writes the solution file

  SmartPtr<const Vector> x = app->IpoptDataObject()->curr()->x();
  SmartPtr<const Vector> y_c = app->IpoptDataObject()->curr()->y_c();
  SmartPtr<const Vector> y_d = app->IpoptDataObject()->curr()->y_d();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));
  SmartPtr<const Matrix> opt_jac_c_p = orig_nlp->jac_c_p(*x);
  SmartPtr<const Matrix> opt_jac_d_p = orig_nlp->jac_d_p(*x);
  SmartPtr<const Matrix> opt_h_p = orig_nlp->h_p(*x, 1.0, *y_c, *y_d);
  opt_jac_c_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
  opt_jac_d_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");
  opt_h_p->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_h_p");

  // Get the (factorized) KKT matrix
  SmartPtr<IpoptAlgorithm> alg = app->AlgorithmObject();
  SmartPtr<PDSearchDirCalculator> pd_search;
  pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
  SmartPtr<PDSystemSolver> pd_solver = pd_search->PDSolver();

  // Compute ds/dp from (del K)/(del s) ds/dp = (del K)/(del p)
  Index np = orig_nlp->p()->Dim();
  SmartPtr<MultiVectorMatrixSpace> mv_space = new MultiVectorMatrixSpace(np, *x->OwnerSpace());
  SmartPtr<MultiVectorMatrix> mv = dynamic_cast<MultiVectorMatrix*>(mv_space->MakeNew());
  Number* dp_values = new Number[np];
  for (int k=0; k<np; ++k) {
    // set up current dp vector as unit vector with entry at index k
    SmartPtr<DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()))->MakeNewDenseVector();
    for (int j=0; j<np; ++j)
      dp_values[j] = 0.0;
    dp_values[k] = 1.0;
    dp->SetValues(dp_values);
    // set up iterates vector and initialize - will be rhs for linear system
    SmartPtr<IteratesVector> it_vec = app->IpoptDataObject()->curr()->MakeNewIteratesVector();
    it_vec->Set(0.0);
    // set up x part of rhs iterates vector
    SmartPtr<Vector> x_it = x->MakeNew();
    opt_h_p->MultVector(1.0, *dp, 0.0, *x_it);
    it_vec->Set_x_NonConst(*x_it);
    // set up c part of rhs iterates vector
    SmartPtr<Vector> c_it = y_c->MakeNew();
    opt_jac_c_p->MultVector(1.0, *dp, 0.0, *c_it);
    it_vec->Set_y_c_NonConst(*c_it);
    // set up d part of rhs iterates vector
    SmartPtr<Vector> d_it = y_d->MakeNew();
    opt_jac_d_p->MultVector(1.0, *dp, 0.0, *d_it);
    it_vec->Set_y_d_NonConst(*d_it);
    // do actual backsolve
    SmartPtr<IteratesVector> lhs = it_vec->MakeNewIteratesVector();
    pd_solver->Solve(1.0, 0.0, *it_vec, *lhs);
    mv->SetVector(k, *lhs->x());
  }
  delete[] dp_values;
  mv->Print(*app->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "dxdp");
  SmartPtr<Matrix> retval(GetRawPtr(mv));
  return retval;
}

SmartPtr<Vector> getDirectionalDerivative(SmartPtr<IpoptApplication> app, SmartPtr<Matrix> sens_matrix) {
  SmartPtr<const IteratesVector> curr = app->IpoptDataObject()->curr();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // if perturbed values are given, compute the step and print it
  SmartPtr<const DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> dp_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(dp->OwnerSpace()));
  SmartPtr<DenseVector> delta_p;
  if (dp_space->HasNumericMetaData("perturbed")) {
    const std::vector<Number> perturbed = dp_space->GetNumericMetaData("perturbed");
    delta_p = dp->MakeNewDenseVector();
    const Number* dp_values = dp->Values();
    Number* new_values = new Number[dp->Dim()];
    for (int k=0; k<dp->Dim(); ++k) {
      new_values[k] = perturbed[k] - dp_values[k];
    }
    delta_p->SetValues(new_values);
    delete[] new_values;
    SmartPtr<Vector> delta_x = curr->x()->MakeNewCopy();

    sens_matrix->MultVector(1.0, *delta_p, 0.0, *delta_x);
    return delta_x;
  }
  else {
    return NULL;
  }
}


bool doIntervallization(SmartPtr<IpoptApplication> app, SmartPtr<AmplSuffixHandler> suffix_handler,  SmartPtr<ParaTNLP> ampl_tnlp)
{

  /////////////////////////////////////////////////////////////////////
  //bewa01 : trying to get a hold of the intervallisation suffix data//
  /////////////////////////////////////////////////////////////////////


  SmartPtr<const IteratesVector> curr = app->IpoptDataObject()->curr();

  SmartPtr<IpoptNLP> ipopt_nlp = app->IpoptNLPObject();
  SmartPtr<OrigIpoptNLP> orig_nlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp));

  // trying to determine parameter positions and intervalIDs in the NLP object
  SmartPtr<const DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> dp_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(dp->OwnerSpace()));

  // set up vars to store ampl tnlp information
  Index n = 0;
  Index np = 0;
  Index m = 0;
  Index nnz_jac_g = 0;
  Index nnz_h_lag = 0;
  Index nnz_jac_g_p = 0;
  Index nnz_h_lag_ = 0;
  AmplTNLP::IndexStyleEnum index_style;
  bool info_status = false;

  // fill the vars with available ampl tnlp information
  info_status = ampl_tnlp->get_nlp_info(n, np, m, nnz_jac_g,
                              nnz_h_lag, nnz_jac_g_p,
  				       nnz_h_lag_,index_style);
  Index nn = n+np;

  // output of variable indexes with or without parameter and intervalID tags
  const Index* parameter = suffix_handler->GetIntegerSuffixValues("parameter",
								 AmplSuffixHandler::Variable_Source);
  std::vector<Index> parameter_vec(nn);
  std::vector<Index> par_index;
  if (!parameter) {
    return 0;  // NO PARAMETERs???? HOW AM I SUPPOSED TO DO MY JOB WITHOUT PARAMETERs???
  }
  std::copy(parameter, parameter+nn, &parameter_vec[0]);
  //var_integer_md["parameter"] = parameter_vec;
  for (Index k_it=0;k_it<nn; k_it++){
    if (parameter_vec[k_it]) {
      // printf("Ipopt variable no %d ist ein Parameter.\n",k_it,parameter_vec[k_it]);
      par_index.push_back(k_it);
    }
  }

  const Index* intervalID = suffix_handler->GetIntegerSuffixValues("intervalID",
								 AmplSuffixHandler::Variable_Source);
  Index int_obj_idx = 0;
  std::vector<Index> intervalID_vec(nn);
  if (!intervalID) {
    return 0;  // NO INTERVAL IDs???? HOW AM I SUPPOSED TO DO MY JOB WITHOUT INTERVAL IDs???
  }
  Index nint_tmp = 0;
  std::copy(intervalID, intervalID+nn, &intervalID_vec[0]);
  //var_integer_md["intervalID"] = intervalID_vec;
  for (Index k_it=0;k_it<nn; k_it++){
    if (intervalID_vec[k_it]){
      //      printf("IntervalID no %d is %d\n",k_it,intervalID_vec[k_it]);
      if (parameter_vec[k_it]!=1 && !int_obj_idx)
	int_obj_idx = k_it;
      if (intervalID_vec[k_it]>nint_tmp)
	nint_tmp = intervalID_vec[k_it];
    }
  }
  const Index nint = nint_tmp;

  const Index* includeID = suffix_handler->GetIntegerSuffixValues("includeID",
								 AmplSuffixHandler::Variable_Source);
  std::vector<Index> includeID_vec(nn);
  if (!includeID) {
    return 0;  // NO INTERVAL IDs???? HOW AM I SUPPOSED TO DO MY JOB WITHOUT INTERVAL IDs???
  }
  Index ninc_tmp = 0;
  std::copy(includeID, includeID+nn, &includeID_vec[0]);
  //var_integer_md["includeID"] = includeID_vec;
  for (Index k_it=0;k_it<nn; k_it++){
    if (includeID_vec[k_it]){
      //      printf("IncludeID no %d is %d\n",k_it,includeID_vec[k_it]);
	ninc_tmp = includeID_vec[k_it];
    }
  }
  const Index ninc = ninc_tmp;

  printf("\n\n\n\n ninc: %d \n\n\n\n",ninc);
 // output of info variable content
  //  printf("\n\n The values of the infovariables are now:\n n = %d \n np = %d \n m = %d \n nnz_jac_g = %d \n nnz_h_lag = %d \n nnz_jac_g_p = %d \n nnz_h_lag_ = %d\n\n",n,np,m,nnz_jac_g,nnz_h_lag,nnz_jac_g_p,nnz_h_lag_);

  AmplTNLP::StringMetaDataMapType var_string_md;
  AmplTNLP::IntegerMetaDataMapType var_integer_md;
  AmplTNLP::NumericMetaDataMapType var_numeric_md;
  AmplTNLP::StringMetaDataMapType para_string_md;
  AmplTNLP::IntegerMetaDataMapType para_integer_md;
  AmplTNLP::NumericMetaDataMapType para_numeric_md;
  AmplTNLP::StringMetaDataMapType con_string_md;
  AmplTNLP::IntegerMetaDataMapType con_integer_md;
  AmplTNLP::NumericMetaDataMapType con_numeric_md;

  bool var_con_metadata_status = 0;
  var_con_metadata_status = ampl_tnlp->get_var_con_metadata(n,var_string_md,var_integer_md,
				var_numeric_md,np,para_string_md,para_integer_md,para_numeric_md,
				m,con_string_md,con_integer_md,con_numeric_md);


  // get parameter names
  const std::vector<std::string> parnames = dynamic_cast<const DenseVectorSpace*>
                       (GetRawPtr(orig_nlp->p()->OwnerSpace()))->GetStringMetaData("idx_names");

  //bewa01: this should always be interchangable with np and hence deleted or edited someday
  Index i_p = parnames.size();
  std::vector<std::string> par_names_tmp;
  for (int i=0;i<i_p;i++)
    par_names_tmp.push_back(parnames[i].c_str());
  const std::vector<std::string> par_names = par_names_tmp;
  // get parameter values
  const Number* p_val = dp->Values();
  std::vector<Number> par_values(i_p);
  std::copy(p_val, p_val+i_p,&par_values[0]);

  // output of gathered parameter information

  /*  for (Index i=0; i<i_p;i++)
    printf("Der %d. Parametereintrag schimpft sich %s. Sein Wert ist %f. Er hat im Ipopt Problem den Variablenplatz %d und die Intervalnummer %d.\n",i+1,parnames[i].c_str(), par_values[i], par_index[i], intervalID_vec[par_index[i]]);
  SmartPtr<const Vector> x = curr->x();
  std::vector<Number> var_values(nn);
  */ //  printf("\n x hat %d Einträge. \n", x->Dim());

  IntervallInfo IntInfo = IntervallInfo(nint, ninc,  par_names, par_values);


  std::vector<std::string> * ipnames = new std::vector<std::string>;
  std::vector<Number> * ipvalues = new std::vector<Number>;
  /*
  IntInfo.GetParameters(&*ipnames,&*ipvalues);
  // int i=0;
    if (ipnames)
     for (int i=0;i<ipnames->size();i++)
      printf("\n Die gespeicherten Namen lauten: %s", ipnames->at(i).c_str());

  if (ipvalues)
         for (int i=0;i<ipvalues->size();i++)
      printf("\n Die gespeicherten Werte lauten: %f", ipvalues->at(i));
  */
  const Index nop =3;
  IntInfo.AddRandomInts(nop);

  IntInfo.GetParameters(&*ipnames,&*ipvalues);
  IntInfo.WriteIntFile();


  /*   // int i=0;
  if (ipnames)
     for (int i=0;i<ipnames->size();i++)
      printf("\n Die gespeicherten Namen lauten: %s", ipnames->at(i).c_str());

  if (ipvalues)
         for (int i=0;i<ipvalues->size();i++)
      printf("\n Die gespeicherten Werte lauten: %f", ipvalues->at(i));
*/
  //get the value of arbitrary objective variable at the solution point
  // Number * p_vv = NULL;
  // bool get_x_status = false;
  // *p_vv = &var_values[0];
  //  get_x_status = ampl_tnlp->get_var_and_para_x(const nn, p_vv);

  // SmartPtr<const IteratesVector> curr = app->IpoptDataobject()->curr();


  //  SmartPtr<const DenseVector> dx = x->MakeNewDenseVector();
  // SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const VectorSpace*>(GetRawPtr(x->Ownerspace()));
  // get var values

  // const Number* v_values = x->Values();

  // const std::vector<std::string> varnames = x_space->GetStringMetaData("idx_names");
    // dynamic_cast<const DenseVectorSpace*>
    //                 (GetRawPtr(orig_nlp->p()->OwnerSpace()))->GetStringMetaData("idx_names");

  //bewa01: this should be interchangable with n and hence deleted or edited someday
  //Index i_v = varnames.size();

  // std::copy(x,x+n,&var_values[0]);
  // std::copy(v_values,v_values+nn,&var_values[0]);
  // get parameter values
  // const Number* p_val = dp->Values();
  // std::vector<Number> par_values(i_p);
  // std::copy(p_val, p_val+i_p,&par_values[0]);


  /*

  SmartPtr<const DenseVector> dp = dynamic_cast<const DenseVector*>(GetRawPtr(orig_nlp->p()));
  SmartPtr<const DenseVectorSpace> dp_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(dp->OwnerSpace()));




    const std::vector<int> pariIDs = dynamic_cast<const DenseVectorSpace*>
                       (GetRawPtr(orig_nlp->p()->OwnerSpace()))->GetIntegerMetaData("intervalID");

 const std::vector<std::string> varnames = dynamic_cast<const DenseVectorSpace*>
                       (GetRawPtr(orig_nlp->x()->OwnerSpace()))->GetStringMetaData("idx_names");
  Index i_v = varnames.size();
  for (Index i=0; i<i_v;i++)
    printf("\n Der %d. Variableneintrag schimpft sich %s.",i+1,varnames[i].c_str());

  */





  //////////////////////////end of bewa striking//////////////
  return 1;
}
