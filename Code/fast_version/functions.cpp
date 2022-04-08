#include <Rcpp.h>

using namespace Rcpp;

NumericVector E_s_Rho_X_cpp(NumericVector beta_rho, NumericMatrix yx_all, int rho_pwr) {
  int nrow = yx_all.nrow(), ncol = yx_all.ncol();
  NumericVector out(nrow);
  
  for (int i = 0; i < nrow; i++) {
    double sumTemp = 0;
    for (int j = 0; j < ncol; j++) {
      sumTemp += exp(rho_pwr*(yx_all(i,j)*beta_rho(0)+yx_all(i, j)*yx_all(i, j)*beta_rho(1)));
    }
    out(i) = sumTemp/ncol;
  }
  return out;
}

double E_s_Rho_cpp(NumericVector beta_rho, NumericMatrix sDat)
{
  int nrow = sDat.nrow();
  NumericVector yVec = sDat( _ , 0);
  double tmpSum = 0, yVal;
  for (int i = 0; i < nrow; i++) {
    yVal = yVec(i);
    tmpSum += exp(yVal*beta_rho(0)+yVal*yVal*beta_rho(1));
  }
  tmpSum /= nrow;
  
  return tmpSum;
}

NumericVector ComputeTau_cpp(double c_ps, NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double p1) {
  int nlen = e_s_rho_x.length();
  NumericVector out(nlen);
  double e_t_rho_x;
  for (int i=0; i<nlen; i++) {
    e_t_rho_x = e_s_rho2_x(i)/e_s_rho_x(i);
    out(i) = (e_t_rho_x/p1/c_ps)/(e_t_rho_x/c_ps/p1+1/(1-p1));
  }
  return out;
}

double E_t_Tau_cpp(int n, int m, double c_ps, NumericVector e_s_rho_x, NumericVector e_s_rho2_x, double p1) {
  NumericVector tauVec = ComputeTau_cpp(c_ps, e_s_rho_x, e_s_rho2_x, p1);
  double out = 0;
  for (int i=0; i<m; i++) {
    out += tauVec(i+n);
  }
  out /= m;
  return out;
}

NumericVector E_t_d_log_Rho_d_Beta_cpp(NumericVector beta_rho, NumericMatrix sDat, double c_ps) {
  int nrow = sDat.nrow();
  NumericVector yVec = sDat( _ , 0);
  
  double rhoVal;
  double out1 = 0, out2 = 0;
  for (int i=0; i<nrow; i++) {
    rhoVal = exp(beta_rho(0)*yVec(i)+beta_rho(1)*yVec(i)*yVec(i));
    out1 += yVec(i)*rhoVal;
    out2 += yVec(i)*yVec(i)*rhoVal;
  }
  out1 /= (nrow*c_ps);
  out2 /= (nrow*c_ps);
  
  NumericVector out(2);
  out(0) = out1;
  out(1) = out2;
  
  return out;
}

NumericMatrix E_t_d_log_Rho_d_Beta_X_cpp(NumericVector beta_rho, NumericMatrix yx_all, NumericVector e_s_rho_x) {
  int i, j;
  int nrow = yx_all.nrow(), ncol = yx_all.ncol();
  double tmp1, tmp2, tmp;
  NumericVector ySample(ncol);
  NumericMatrix out(nrow, 2);
  
  for (i=0; i<nrow; i++) {
    ySample = yx_all(i, _ );
    tmp1 = 0;
    tmp2 = 0;
    for (j=0; j<ncol; j++) {
      tmp = exp(beta_rho(0)*ySample(j)+beta_rho(1)*ySample(j)*ySample(j));
      tmp1 += ySample(j)*tmp;
      tmp2 += ySample(j)*ySample(j)*tmp;
    }
    tmp1 /= (ncol*e_s_rho_x(i));
    tmp2 /= (ncol*e_s_rho_x(i));
    
    out(i, 0) = tmp1;
    out(i, 1) = tmp2;
  }
  
  return out;
}

NumericMatrix ComputeS_cpp(NumericVector beta_rho, NumericMatrix yx_all, NumericVector e_s_rho_x, NumericMatrix sDat, double c_ps) {
  NumericMatrix e_t_d_log_rho_d_beta_x = E_t_d_log_Rho_d_Beta_X_cpp(beta_rho, yx_all, e_s_rho_x);
  NumericVector e_t_d_log_rho_d_beta = E_t_d_log_Rho_d_Beta_cpp(beta_rho, sDat, c_ps);
  
  int nrow = e_t_d_log_rho_d_beta_x.nrow();
  for (int i=0; i<nrow; i++) {
    e_t_d_log_rho_d_beta_x(i,0) = e_t_d_log_rho_d_beta_x(i,0)-e_t_d_log_rho_d_beta(0);
    e_t_d_log_rho_d_beta_x(i,1) = e_t_d_log_rho_d_beta_x(i,1)-e_t_d_log_rho_d_beta(1);
  }
  
  return e_t_d_log_rho_d_beta_x;
}

NumericVector S_Eff_Multiplier_cpp(NumericVector beta_rho, int n, int m, NumericVector yn, double c_ps) {
  double p1 = n / (double) (n+m);
  NumericVector out(n+m);
  for (int i=0; i<n; i++) {
    out(i) = 1/p1*exp(beta_rho(0)*yn(i)+beta_rho(1)*yn(i)*yn(i))/c_ps;
  }
  for (int i=0; i<m; i++) {
    out(i+n) = -1/(1-p1);
  }
  
  return out;
}

NumericMatrix ComputeB1_cpp(NumericVector beta_rho, NumericMatrix yx_all,
                            NumericVector e_s_rho_x, NumericVector e_s_rho2_x,
                            double c_ps, int n, int m,
                            NumericMatrix sDat) {
  double p1 = n/(double) (n+m);
  NumericVector tauVec = ComputeTau_cpp(c_ps, e_s_rho_x, e_s_rho2_x, p1);
  NumericMatrix sVec = ComputeS_cpp(beta_rho, yx_all, e_s_rho_x, sDat, c_ps);
  
  NumericVector e_t_tau_s(2);
  double e_t_tau_s_1=0, e_t_tau_s_2=0;
  double e_t_tau = E_t_Tau_cpp(n, m, c_ps, e_s_rho_x, e_s_rho2_x, p1);
  
  NumericMatrix out(n+m, 2);
  
  for (int i=0; i<m; i++) {
    e_t_tau_s_1 += tauVec(i+n)*sVec(i+n,0);
    e_t_tau_s_2 += tauVec(i+n)*sVec(i+n,1);
  }
  
  e_t_tau_s_1 /= (m*(e_t_tau-1));
  e_t_tau_s_2 /= (m*(e_t_tau-1));
  
  for (int i=0; i<n+m; i++) {
    out(i,0) = -(1-p1)*(1-tauVec(i))*(sVec(i,0)-e_t_tau_s_1);
    out(i,1) = -(1-p1)*(1-tauVec(i))*(sVec(i,1)-e_t_tau_s_2);
  }
  
  return out;
}

NumericMatrix ComputeSEff_cpp(NumericVector beta_rho, NumericMatrix yx_all,
                              NumericVector e_s_rho_x, NumericVector e_s_rho2_x,
                              double c_ps, int n, int m, NumericVector yn,
                              NumericMatrix sDat) {
  NumericVector smulti = S_Eff_Multiplier_cpp(beta_rho, n, m, yn, c_ps);
  NumericMatrix b1Mat = ComputeB1_cpp(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, n, m, sDat);
  
  NumericMatrix out(n+m, 2);
  for (int i=0; i<n+m; i++) {
    out(i,0) = smulti(i)*b1Mat(i,0);
    out(i,1) = smulti(i)*b1Mat(i,1);
  }
  
  return out;
}

// [[Rcpp::export]]
double ComputeEquation_cpp(NumericVector beta_rho, NumericMatrix yx_all, NumericMatrix sDat, int n, int m, double p1) {
  NumericVector e_s_rho_x = E_s_Rho_X_cpp(beta_rho, yx_all, 1), e_s_rho2_x = E_s_Rho_X_cpp(beta_rho, yx_all, 2);
  double c_ps = E_s_Rho_cpp(beta_rho, sDat);
  NumericVector yVec = sDat( _ , 0);
  NumericMatrix sEffMat = ComputeSEff_cpp(beta_rho, yx_all, e_s_rho_x, e_s_rho2_x, c_ps, n, m, yVec, sDat);
  
  double out1=0, out2=0, out;
  for (int i=0; i<n+m; i++) {
    out1 += sEffMat(i,0);
    out2 += sEffMat(i,1);
  }
  out1 /= n+m;
  out2 /= n+m;
  out = out1*out1+out2*out2;
  
  return out;
}

// double my_f(const gsl_vector *v, void *params) {
  
//   double beta1, beta2;
//   beta1 = gsl_vector_get(v, 0);
//   beta2 = gsl_vector_get(v, 1);
  
//   NumericVector beta_vec(2);
//   beta_vec(0) = beta1;
//   beta_vec(1) = beta2;
  
//   struct optim_params *dat = (struct optim_params *) params;
//   NumericMatrix yx_all = dat->yx_all;
//   NumericMatrix sDat = dat->sDat;
//   int n = dat->n, m = dat->m;
//   double p1 = dat->p1;
  
//   double out;
//   out = ComputeEquation_cpp(beta_vec, yx_all, sDat, n, m, p1);
  
//   return out;
// }

// // [[Rcpp::export]]

// NumericVector findBetaLabelShift(List dat) {
//   NumericMatrix yx_all = dat[0];
//   NumericMatrix sDat = dat[1];
//   int n = dat[2], m = dat[3];
//   double p1 = dat[4];
//   NumericVector beta_init = dat[5];
  
//   struct optim_params par = {yx_all, sDat, n, m, p1};
  
//   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
//   gsl_multimin_fminimizer *s = NULL;
//   gsl_vector *ss, *x;
//   gsl_multimin_function minex_func;
  
//   int iter = 0;
//   int status;
//   double size;
  
//   x = gsl_vector_alloc (2);
  
//   gsl_vector_set (x, 0, beta_init(0));
//   gsl_vector_set (x, 1, beta_init(1));
  
//   ss = gsl_vector_alloc (2);
//   gsl_vector_set_all (ss, 1.0);
  
//   /* Initialize method and iterate */
//   minex_func.n = 2;
//   minex_func.f = my_f;
//   minex_func.params = &par;
  
//   s = gsl_multimin_fminimizer_alloc (T, 2);
//   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
//   do
//   {
//     iter++;
//     status = gsl_multimin_fminimizer_iterate(s);
    
//     if (status)
//       break;
    
//     size = gsl_multimin_fminimizer_size (s);
//     status = gsl_multimin_test_size (size, 1e-5);
//     // printf("%5d %.5f %.5f %10.5f\n",
//     //        iter,
//     //        gsl_vector_get (s->x, 0),
//     //        gsl_vector_get (s->x, 1),
//     //        gsl_multimin_fminimizer_minimum(s));
//   }
//   while (status == GSL_CONTINUE && iter < 500);
  
//   NumericVector MLEs(2);
//   if (status == GSL_SUCCESS)
//   {
//     MLEs[0] = gsl_vector_get (s->x, 0);
//     MLEs[1] = gsl_vector_get (s->x, 1);
//   }
//   else
//   {
//     MLEs[0] = -99;
//     MLEs[1] = -99;
//   }
  
//   gsl_vector_free(x);
//   gsl_vector_free(ss);
//   gsl_multimin_fminimizer_free (s);
  
//   return MLEs;
// }
