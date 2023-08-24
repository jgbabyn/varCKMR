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
Type objective_function<Type>::operator() ()
{

    PARAMETER(log_V_A);
    PARAMETER(log_N);
    PARAMETER(log_E_A);
  
    Type V_A = exp(log_V_A);
    Type N = exp(log_N);
    Type E_A = exp(log_E_A);
  
    DATA_IMATRIX(datan);
    Type nll = 0.0;

    vector<Type> exppairs(datan.rows());
    for(int i = 0; i < datan.rows(); ++i){
      Type n = Type(datan(i,0));
      Type obs = Type(datan(i,4));
      Type pairs = smooth_choose(n,Type(2.0),false);
      Type expp = 1/N*(1+1/V_A)*pairs;
      exppairs(i) = expp;
      nll -= dpois(obs,expp,true);
    }

    REPORT(exppairs);
    ADREPORT(exppairs);
    REPORT(V_A);
    REPORT(N);
    REPORT(E_A);
    return nll;

}
