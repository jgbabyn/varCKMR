#include <TMB.hpp>

template<class Type>
Type smooth_choose(Type n, Type k,bool give_log){
  Type leftover = n-k;
  Type ret = lfactorial(n) - (lfactorial(leftover)+lfactorial(k));

  if(give_log == false){
    ret = exp(ret);
  }

  return ret;
}

template<class Type>
Type smooth_dhyper(Type x, Type m,Type n, Type k, bool give_log){
  Type tot = m+n;
  Type ret = (smooth_choose(m,x,true)+smooth_choose(n,k-x,true)) - smooth_choose(tot,k,true);
  if(give_log == false){
    ret = exp(ret);
  }
  return ret;
}

template<class Type>
Type exp_pairs_given_nb(Type n,Type N,Type lambda,Type theta){
  vector<Type> densnb(2001);
  vector<Type> fecprop(2001);

  Type nb_var = lambda+pow(lambda,2)/theta;

  for(int i = 0; i < densnb.size(); ++i){
    densnb(i) = dnbinom2(Type(i),lambda,nb_var);
    fecprop(i) = N*densnb(i);
  }

  Type totballs = N*lambda;

  vector<Type> exppairs(2001);
  Type exp_res = 0.0;
  for(int x = 0; x < densnb.size(); ++x){
    exp_res += fecprop(x)*((n*x)/(totballs)-smooth_dhyper(Type(1.0),Type(x),totballs-x,n,false));
  }
  return exp_res;
}



template<class Type>
Type objective_function<Type>::operator() ()
{
  PARAMETER(log_theta);
  PARAMETER(log_N);
  PARAMETER(log_lambda);

  Type theta = exp(log_theta);
  Type N = exp(log_N);
  Type lambda = exp(log_lambda);

  DATA_IMATRIX(datan);
  Type nll = 0.0;

  vector<Type> exppairs(datan.rows());
  for(int i = 0; i < datan.rows(); ++i){
    Type n = Type(datan(i,0));
    Type obs = Type(datan(i,4));
    Type expp = exp_pairs_given_nb(n,N,lambda,theta);
    exppairs(i) = expp;
    nll -= dpois(obs,expp,true);
  }

  REPORT(exppairs);
  ADREPORT(exppairs);
  return nll;


}
