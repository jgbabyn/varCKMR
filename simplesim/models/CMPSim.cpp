#include <TMB.hpp>

//used to easily get pmf from CMP dist with given mean and overdispersion.
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    PARAMETER(mean);
    PARAMETER(nu);
    DATA_INTEGER(n);

    vector<Type> sims(n);
    SIMULATE{
      sims = rcompois2(n,mean,nu);
      REPORT(sims);
    }

    return 0;
  }
