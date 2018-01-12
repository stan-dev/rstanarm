#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4bernoulli_mod) {


    class_<rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> >("model_bernoulli")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_bernoulli_namespace::model_bernoulli, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4binomial_mod) {


    class_<rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> >("model_binomial")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_binomial_namespace::model_binomial, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4continuous_mod) {


    class_<rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> >("model_continuous")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_continuous_namespace::model_continuous, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4count_mod) {


    class_<rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> >("model_count")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_count_namespace::model_count, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4lm_mod) {


    class_<rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> >("model_lm")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_lm_namespace::model_lm, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4polr_mod) {


    class_<rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> >("model_polr")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_polr_namespace::model_polr, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4spatial_mod) {


    class_<rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> >("model_spatial")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_spatial_namespace::model_spatial, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
