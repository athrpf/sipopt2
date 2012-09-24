Ipopt trunk: 

******************************************************************************
This program contains Ipopt, a library for large-scale nonlinear optimization.
 Ipopt is released as open source code under the Eclipse Public License (EPL).
         For more information visit http://projects.coin-or.org/Ipopt
******************************************************************************

NOTE: You are using Ipopt by default with the MUMPS linear solver.
      Other linear solvers might be more efficient (see Ipopt documentation).


This is Ipopt version trunk, running with linear solver mumps.

Number of nonzeros in equality constraint Jacobian...:     3591
Number of nonzeros in inequality constraint Jacobian.:        2
Number of nonzeros in Lagrangian Hessian.............:     2160

Number of nonzeros in equality parameter Jacobian..:      720
Number of nonzeros in inequality parameter Jacobian:        0
Number of nonzeros in parameter Lagrangian Hessian.:      964

Total number of variables............................:      483
                     variables with only lower bounds:        1
                variables with lower and upper bounds:      480
                     variables with only upper bounds:        1
Total number of equality constraints.................:      481
Total number of inequality constraints...............:        1
        inequality constraints with only lower bounds:        0
   inequality constraints with lower and upper bounds:        0
        inequality constraints with only upper bounds:        1

Total number of parameters...........................:        4

iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
   0  0.0000000e+00 9.90e-01 9.77e-01  -1.0 0.00e+00    -  0.00e+00 0.00e+00   0
   1  1.7359560e-03 7.85e-01 1.70e+02  -1.7 3.10e+01    -  1.01e-02 2.03e-01f  1
   2  1.9419467e-03 7.83e-01 1.70e+02  -1.7 2.87e+01    -  1.15e-01 1.73e-03h  1
   3 -1.0567039e-01 6.86e-01 5.15e+02  -1.7 2.83e+01    -  1.16e-02 1.26e-01f  1
   4 -1.1785856e-01 6.63e-01 4.89e+02  -1.7 8.52e+00    -  3.37e-04 3.29e-02f  1
   5 -2.8010797e-01 3.78e-01 5.93e+02  -1.7 7.32e+00    -  4.25e-02 4.57e-01f  1
   6 -6.8870719e-02 1.18e-02 1.92e+02  -1.7 5.89e-01    -  5.69e-02 9.94e-01h  1
   7 -1.9786141e-01 9.08e-04 6.70e+01  -1.7 1.54e-01    -  5.56e-01 1.00e+00f  1
   8 -2.1309517e-01 9.32e-06 6.36e-01  -1.7 1.86e-02    -  1.00e+00 1.00e+00h  1
   9 -2.1289171e-01 7.50e-09 2.01e-04  -1.7 4.41e-04    -  1.00e+00 1.00e+00h  1
iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls
  10 -2.3274369e-01 9.38e-12 1.96e-04  -3.8 1.99e-02    -  1.00e+00 1.00e+00f  1
  11 -2.3302794e-01 2.67e-08 7.86e-05  -5.7 9.52e-04    -  1.00e+00 1.00e+00f  1
  12 -2.3311456e-01 1.12e-08 1.69e-07  -5.7 6.16e-04    -  1.00e+00 1.00e+00h  1
  13 -2.3311824e-01 6.19e-12 1.49e-08  -8.6 1.46e-05    -  1.00e+00 1.00e+00h  1
  14 -2.3311826e-01 5.27e-16 1.33e-13  -9.0 1.32e-07    -  1.00e+00 1.00e+00h  1

Number of Iterations....: 14

                                   (scaled)                 (unscaled)
Objective...............:  -2.3311825758574073e-01   -2.3311825758574073e-01
Dual infeasibility......:   1.3342052078713678e-13    1.3342052078713678e-13
Constraint violation....:   5.2735593669694936e-16    5.2735593669694936e-16
Complementarity.........:   9.0923476497266106e-10    9.0923476497266106e-10
Overall NLP error.......:   9.0923476497266106e-10    9.0923476497266106e-10


Number of objective function evaluations             = 15
Number of objective gradient evaluations             = 15
Number of equality constraint evaluations            = 15
Number of inequality constraint evaluations          = 15
Number of equality constraint Jacobian evaluations   = 15
Number of inequality constraint Jacobian evaluations = 15
Number of Lagrangian Hessian evaluations             = 14
Total CPU secs in IPOPT (w/o function evaluations)   =      0.276
Total CPU secs in NLP function evaluations           =      0.136

EXIT: Optimal Solution Found.

 AmplTNLP::doIntervallization started.




 ninc: 1 




 intervalls.inc has been opened successfully. 

 AmplTNLP::doIntervallization closed.
 
Ipopt trunk: Optimal Solution Found
time = 0.481428

xCU [*] :=
1  0.2
;

