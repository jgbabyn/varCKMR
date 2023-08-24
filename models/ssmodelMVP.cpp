#include <TMB.hpp>
#include <Eigen/Eigenvalues>

template<class Type>
vector<Type>  get_ls(vector<Type> &Z){
  vector<Type> surv = exp(-Z);
  vector<Type> ls(Z.size());
  ls(0) = 1;
  for(int a = 1; a < surv.size();++a){
    ls(a) = ls(a-1)*surv(a-1);
  }
  return ls;
}

template<class Type>
vector<Type> death_probs(vector<Type> &Z){
  vector<Type> surv = exp(-Z);
  vector<Type> dprobs(surv.size());
  vector<Type> ells = get_ls(Z);
  for(int a = 0; a < surv.size()-1; ++a){
    dprobs(a) = ells(a)*(1-surv(a));
  }
  dprobs(surv.size()-1) = ells(surv.size()-1);
  return dprobs;
}

template<class Type>
Type lifetime_reproductive_mean(vector<Type> &Z, vector<Type> &fec){
  vector<Type> dprob = death_probs(Z);
  vector<Type> fecS(fec.size());
  fecS(0) = fec(0);
  for(int a = 1; a < fec.size(); ++a){
    fecS(a) = fec(a)+fecS(a-1);
  }  
  Type mean_fec = 0;
  for(int a = 0; a < fec.size(); ++a){
    mean_fec += dprob(a)*fecS(a);
  }
  return mean_fec;
}

//The good one
template<class Type>
Type lifetime_reproductive_var2(vector<Type> &Z, vector<Type> &fec, Type theta){
  vector<Type> surv = exp(-Z);
  vector<Type> dprob = death_probs(Z);
  vector<Type> nb_var(fec.size());
  vector<Type> nb_varS(fec.size());
  vector<Type> fecS(fec.size());
  for(int a = 0; a < fec.size(); ++a){
    nb_var(a) = fec(a)+pow(fec(a),2)/theta;
  }

  nb_varS(0) = nb_var(0);
  fecS(0) = fec(0);
  for(int a = 1; a < fec.size(); ++a){
    nb_varS(a) = nb_var(a)+nb_varS(a-1);
    fecS(a) = fec(a)+fecS(a-1);
  }

  Type firstsum = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    firstsum += dprob(a)*nb_varS(a);
  }

  Type secsum = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    secsum += pow(fecS(a),2)*(1-dprob(a))*dprob(a);
  }

  Type lastsum = 0.0;
  for(int a = 1; a < fec.size(); ++a){
    for(int j = 0; j < a; ++j){
      lastsum += fecS(a)*dprob(a)*fecS(j)*dprob(j);
    }
  }

  Type tot = firstsum+secsum-2*lastsum;
  return tot;

}


template<class Type>
Type TRO(matrix<Type> &log_N,vector<Type> &fec,int juv_by){
  vector<Type> log_Nyj = log_N.col(juv_by);
  Type TRO = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    TRO += fec(a)*exp(log_Nyj(a));
  }
  return TRO;
}

template<class Type>
vector<Type> RO(matrix<Type> &log_N,vector<Type> &fec,int juv_by){
  vector<Type> log_Nyj = log_N.col(juv_by);
  vector<Type> RO(log_Nyj.size());
  for(int a = 0; a < fec.size(); ++a){
    RO(a) = fec(a)*exp(log_Nyj(a));
  }
  return RO;
}


template<class Type>
Type yearly_variance(vector<Type> &fec,vector<Type> nb_var,matrix<Type> &log_N, int by){
    int A = log_N.rows();
    vector<Type> ROy = RO(log_N,fec,by+1);
    Type TROy = TRO(log_N,fec,by+1);

    vector<Type> Nya(A);
    Type Ny = 0.0;
    for(int a = 0; a < A; ++a){
      Nya(a) = exp(log_N(a,by));
      Ny += exp(log_N(a,by));
    }

    Type firstsum = 0.0;
    for(int a = 0; a < A; ++a){
      firstsum += nb_var(a)*(Nya(a)/Ny);
    }

    Type secsum = 0.0;
    for(int a = 0; a < A; ++a){
      secsum += pow(fec(a),2)*(1-Nya(a)/Ny)*(Nya(a)/Ny);
    }

    Type lastsum = 0.0;
    for(int a = 1; a < A; ++a){
      for(int j = 0; j < a; ++j){
	lastsum += fec(a)*(Nya(a)/Ny)*fec(j)*(Nya(j)/Ny);
      }
    }

    Type tot = firstsum + secsum-2*lastsum;
    return tot;

}

template<class Type>
vector<Type> reproductive_value(vector<Type> &fec,vector<Type> &Z,Type growth_rate){
  vector<Type> surv = exp(-Z);
  vector<Type> v_is(fec.size()+1);
  for(int i = 0; i < fec.size()+1; ++i){
    v_is(i) = 0.0;
  }
  v_is(0) = 1.0;
  Type l_i = 0.0;
  for(int a = 1; a < fec.size(); ++a){
    for(int j = a; j < fec.size(); ++j){
      Type prod = 1.0;
      for(int h = 0; h < j; ++h){
	prod *= surv(h);
	if(h == j-1){
	  l_i = prod;
	}
      }
      v_is(a) += prod*fec(j)*pow(growth_rate,-j-1);

    }
    v_is(a) *= pow(growth_rate,a)/l_i;
  }
  return v_is;
}


template<class Type>
matrix<Type> make_leslie(vector<Type> &fec, vector<Type> &Z){
  vector<Type> surv = exp(-Z);
  matrix<Type> leslie(fec.size(),fec.size());
  leslie.setZero();
  leslie.row(0) = fec;
  for(int i = 0; i < fec.size()-1; ++i){
    leslie(i+1,i) = surv(i);
  }
  return leslie;
}

template<class Type>
Type gen_length(vector<Type> &fec, vector<Type> &Z,Type growth_rate){
  vector<Type> surv = exp(-Z);
  Type top = 0.0;
  Type bot = 0.0;
  for(int i = 0; i < fec.size();++i){
    Rcout << "i is: " << i << std::endl;
    Type prod = 1.0;
    for(int j = 0; j < i; ++j){
      prod *= surv(j);
    }
    //Rcout << "prod is: " << prod << std::endl;
    top += (i+1)*fec(i)*prod*pow(growth_rate,-(i+1));
    //bot += fec(i)*prod;
  }
  Type ret = top;
  return ret;
}


template<class Type>
Type weighted_lambda(vector<Type> &log_Nyj, vector<Type> &fec){
  Type w_lambda = 0.0;
  Type tot_Nyj = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    tot_Nyj += exp(log_Nyj(a));
  }
  for(int a = 0; a < fec.size(); ++a){
    w_lambda += fec(a)*exp(log_Nyj(a));
  }
  w_lambda = w_lambda/tot_Nyj;
  return w_lambda;
}

template<class Type>
vector<Type> init_N(vector<Type> &Z,Type total_init_N){
  vector<Type> ret_N(Z.size());
  ret_N(0) = 1.0;
  for(int a = 1; a < Z.size();++a){
    ret_N(a) = ret_N(a-1)*exp(-Z(a-1));
  }
  ret_N /= ret_N.sum();
  ret_N *= total_init_N;
  return ret_N;
}


template<class Type>
Type Pr_POP(matrix<Type> &log_N,vector<Type> &fec, int j_by,int p_by,int p_sy){

  //get total number of individals of sex j_sex born from parents of sex p_sex
  vector<Type> log_Nyj = log_N.col(j_by);
  Type TROyj = TRO(log_N,fec,j_by);

  int p_age_yj = j_by-p_by-1;
  Type eprob = (2*fec(p_age_yj))/(TROyj);

  //If the parent is sampled before the juvenile is born then they have to survive
  if(p_sy < j_by){
    int p_age_sy = p_sy-p_by-1;
    vector<Type> log_Nsy = log_N.col(p_sy);
    //Just the number at age in juv birth year over number at age in parent sample year
    Type Pr_surv = exp(log_Nyj(p_age_yj))/exp(log_Nsy(p_age_sy));
    eprob *= Pr_surv;
  }

  return eprob;

}

template<class Type>
Type Pr_HSP(matrix<Type> &log_N,vector<Type> &fec, int i_by, int o_by){
  int A = log_N.rows();
  int a_diff = i_by-o_by;
  Type eprob = 0.0;

  for(int d = a_diff; d < A; ++d){
    vector<Type> log_Nyi = log_N.col(i_by);
    vector<Type> log_Nyo = log_N.col(o_by);

    Type TROi = TRO(log_N,fec,i_by);
    Type TROo = TRO(log_N,fec,o_by);

    //parents age at o's birth
    //I think this only works because maturity is age 1?
    int p_age_o = d - a_diff;

    Type prob_i = (2*fec(d))/(TROi);
    Type prob_o = (2*fec(p_age_o))/(TROo);
    eprob += exp(log_Nyo(p_age_o))*prob_o*prob_i;

    //if the age difference is greater than 0 the parent must survive to i_by
    if(a_diff > 0){
      Type Pr_surv = exp(log_Nyi(d))/exp(log_Nyo(p_age_o));
      eprob *= Pr_surv;
    }

  }

  return eprob;
}  

template<class Type>
Type Pr_HSPSC(matrix<Type> &log_N,vector<Type> &fec,vector<Type> wlamb,vector<Type> &nb_var,int by){
  int A = log_N.rows();
  Type var_y = yearly_variance(fec,nb_var,log_N,by);

  //find yearly mean 
  vector<Type> ROy = RO(log_N,fec,by+1);
  Type Nytot = 0.0;
  for(int a = 0; a < A; ++a){
    Nytot += exp(log_N(a,by));
  }

  Type mean_y = wlamb[by];



  Type eprob = 0.0;

    eprob = (1/(Nytot/2.0))*(1+(var_y-mean_y)/pow(mean_y,2));


  return eprob;

}


  ///Because I sillyly hard coded this to use NB variance I need a new version
template<class Type>
Type lifetime_reproductive_var3(vector<Type> &Z, vector<Type> &fec, Type theta,Type var_power){
  vector<Type> surv = exp(-Z);
  vector<Type> dprob = death_probs(Z);
  vector<Type> nb_var(fec.size());
  vector<Type> nb_varS(fec.size());
  vector<Type> fecS(fec.size());
  for(int a = 0; a < fec.size(); ++a){
    nb_var(a) = theta*pow(fec(a),var_power);
  }

  nb_varS(0) = nb_var(0);
  fecS(0) = fec(0);
  for(int a = 1; a < fec.size(); ++a){
    nb_varS(a) = nb_var(a)+nb_varS(a-1);
    fecS(a) = fec(a)+fecS(a-1);
  }

  Type firstsum = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    firstsum += dprob(a)*nb_varS(a);
  }

  Type secsum = 0.0;
  for(int a = 0; a < fec.size(); ++a){
    secsum += pow(fecS(a),2)*(1-dprob(a))*dprob(a);
  }

  Type lastsum = 0.0;
  for(int a = 1; a < fec.size(); ++a){
    for(int j = 0; j < a; ++j){
      lastsum += fecS(a)*dprob(a)*fecS(j)*dprob(j);
    }
  }

  Type tot = firstsum+secsum-2*lastsum;
  return tot;

}


  template<class Type>
  Type objective_function<Type>::operator() ()
  {

    //Parameters
    PARAMETER(log_init_tot_N);
    PARAMETER_VECTOR(log_M);
    PARAMETER_VECTOR(log_fecundity);
    //Theta here is a constant c
    PARAMETER(log_theta);
    PARAMETER(var_power);

    //Transformations
    Type init_tot_N = exp(log_init_tot_N);
    vector<Type> Z = exp(log_M);
    vector<Type> fecundity = exp(log_fecundity);
    Type theta = exp(log_theta);

    //Data
    DATA_INTEGER(Y);
    DATA_INTEGER(A);
    DATA_INTEGER(skipsame);

    DATA_IMATRIX(POPs);
    DATA_IMATRIX(SIBsXC);
    DATA_IMATRIX(SIBsSC);
    DATA_IMATRIX(SIBsSCP);

    //Pop Dynamics
    matrix<Type> log_N(A,Y);
    matrix<Type> N(A,Y);
    N = log_N.array().exp();

    //Death matrix
    matrix<Type> D(A,Y);

    vector<Type> init_N_v = init_N(Z,init_tot_N);
    log_N.col(0) = log(init_N_v);

    for(int y = 1; y < Y; ++y){
      Type TROy = TRO(log_N,fecundity,y-1);
      log_N(0,y) = log(TROy);
      for(int a = 1; a < A; ++a){
	log_N(a,y) = log_N(a-1,y-1)-Z(a-1);
      }
    }

    for(int y = 0; y < Y; ++y){
      for(int a = 0; a < A; ++a){
	D(a,y) = exp(log_N(a,y))-exp(log_N(a,y)-Z(a));
	if(a == A-1){
	  D(a,y) = exp(log_N(a,y));
	}
      }
    }

    //Average lambdas each year, weighted by pop
    vector<Type> w_lambdas(Y);
    //The variance of reproductive success at age
    //here it's not really neg bin...
    vector<Type> nb_var(A);

    for(int a = 0; a < A; ++a){
      nb_var(a) = theta*pow(fecundity(a),var_power);
    }

    for(int y = 0; y < Y; ++y){
      vector<Type> log_Ny = log_N.col(y);
      w_lambdas(y) = weighted_lambda(log_Ny,fecundity);
    }

    Type nll = 0.0;

    //Observation Model

    //POPs
    for(int i = 0; i < POPs.rows(); ++i){
      int juv_by = POPs(i,0);
      int juv_sy = POPs(i,1);
      int par_by = POPs(i,2);
      int par_sy = POPs(i,3);
      Type n_comps = Type(POPs(i,4));
      Type n_POPs = Type(POPs(i,5));

      Type eprob = Pr_POP(log_N,fecundity,juv_by,par_by,par_sy);
      Type eobs = eprob*n_comps;
      nll -= dpois(n_POPs,eobs,true);

    }

    //Sibs Cross Cohort
    for(int i = 0; i < SIBsXC.rows(); ++i){
      int i_by = SIBsXC(i,0);
      int i_sy = SIBsXC(i,1);
      int o_by = SIBsXC(i,2);
      int o_sy = SIBsXC(i,3);
      Type n_comps = Type(SIBsXC(i,4));
      Type n_SIBsXC = Type(SIBsXC(i,5));

      Type eprob = Pr_HSP(log_N,fecundity,i_by,o_by);
      Type eobs = eprob*n_comps;
      nll -= dpois(n_SIBsXC,eobs,true);

    }

    //Sibs Same Cohort
    vector<Type> SSCexp(SIBsSCP.rows());
    vector<Type> SSCprob(SIBsSCP.rows());
    vector<Type> varys(SIBsSCP.rows());
      for(int i = 0; i < SIBsSCP.rows(); ++i){
      int i_by = SIBsSCP(i,0);
      varys(i) = yearly_variance(fecundity,nb_var,log_N,i_by);
      Type n_comps = Type(SIBsSCP(i,2));
      Type n_SIBsSCP = Type(SIBsSCP(i,1));

      Type eprob = Pr_HSPSC(log_N,fecundity,w_lambdas,nb_var,i_by);
      Type eobs = eprob*n_comps;
      eobs = 2.0*eobs;
      SSCexp(i) = eobs;
      SSCprob(i) = eprob;
      if(skipsame == false){
	nll -= dpois(n_SIBsSCP,eobs,true);
      }

    }

      REPORT(varys);
      REPORT(SSCexp);
      REPORT(SSCprob);


    matrix<Type> les = make_leslie(fecundity,Z);
    REPORT(les);

    using namespace Eigen;
    ComplexEigenSolver<Matrix<Type, Dynamic, Dynamic> > es(les);
    vector<Type> EvR = es.eigenvalues().real();
    vector<Type> EvRA = EvR.abs();
    Type growth_rate = max(EvRA);
    matrix<Type> VR = es.eigenvectors().real();

    ComplexEigenSolver<Matrix<Type, Dynamic, Dynamic> > es2(les.transpose());
    vector<Type> EvRL = es2.eigenvalues().real();
    vector<Type> EvRAL = EvRL.abs();
    matrix<Type> VL = es2.eigenvectors().real();
    vector<Type> rp = VL.col(A-1);
    vector<Type> rp2 = rp/rp(0);
    vector<Type> rp3 = rp2/sum(rp2);


    vector<Type> repro_value = reproductive_value(fecundity,Z,growth_rate);
    REPORT(repro_value);
    vector<Type> ls = get_ls(Z);
    REPORT(ls);
    REPORT(growth_rate);
    vector<Type> var_a(A);
    for(int a = 0; a < A; ++a){
      var_a(a) = 2*fecundity(a)+pow(2*fecundity(a),2)/theta;
    }
    vector<Type> var_ab(A);
    var_ab(0) = var_a(0);
    for(int a = 1; a < A; ++a){
      var_ab(a) = var_a(a)+var_ab(a-1);
    }
    REPORT(var_a);
    REPORT(var_ab);
    vector<Type> surv = exp(-Z);
    REPORT(surv);


    //Adjusted fecundity because you gotta account for the number of kids everybody has
    vector<Type> doub_fec(fecundity.size());
    for(int a = 0; a < fecundity.size(); ++a){
      doub_fec(a) = 2.0*fecundity(a);
    }


    Type Vk = lifetime_reproductive_var3(Z,doub_fec,theta,var_power);
    Type kbar = lifetime_reproductive_mean(Z,doub_fec);
    REPORT(kbar);
    ADREPORT(kbar);
    REPORT(Vk);
    ADREPORT(Vk);

    //N_e guesstimate
    Type glen = gen_length(fecundity,Z,growth_rate);
    REPORT(glen);
    vector<Type> N_e(Y);
    for(int y = 0; y < Y; ++y){
      Type Nby = exp(log_N(0,y));
      Type topS = (2.0*Nby-1)*kbar*glen;
      Type botS = 2.0*(1.0+Vk/kbar);
      N_e(y) = topS/botS;
    }



    REPORT(EvR);
    REPORT(D);
    REPORT(VR);
    REPORT(VL);
    REPORT(rp2);
    REPORT(rp3);
    REPORT(log_N);
    REPORT(N);
    REPORT(fecundity);
    REPORT(theta);
    REPORT(Z);
    REPORT(w_lambdas);
    REPORT(nb_var);
    REPORT(N_e);
    ADREPORT(N_e);
    ADREPORT(growth_rate);
    ADREPORT(log_N);
    //REPORT(same_cohort_exp);
    return nll;


  }
