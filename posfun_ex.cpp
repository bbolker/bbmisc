// https://github.com/kaskr/adcomp/issues/7
// example from James Thorson
#include <TMB.hpp>
template<class Type>
Type posfun(Type x, Type eps, Type &pen) {
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  Type xp = -(x/eps-1);
  return CppAD::CondExpGe(x,eps,x,
			  eps*(1/(1+xp+pow(xp,2)+pow(xp,3)+pow(xp,4)+pow(xp,5))));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_SCALAR(eps);
    DATA_VECTOR(x);

    PARAMETER(p);
    PARAMETER(Dummy);
    Type pen;
    pen=0;

    Type nll;
    nll=0;

    Type var;
    var = p*(1-p);
    var = posfun(var, eps, pen);
    nll += pen;

    nll += -sum( dnorm(x, p, sqrt(var),true) );
    nll += -dnorm(Dummy, Type(0.0), Type(1.0), true);
    return nll;
}
