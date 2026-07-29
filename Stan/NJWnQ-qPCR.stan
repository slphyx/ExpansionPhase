functions {

  // PDF of normal distribution
  real PDF_normal(real mu, real sig, real x) {
    return exp(-0.5 * square((x - mu) / sig)) / (sqrt(2 * pi()) * sig);
  }

  // for generating the initial parasite distribution
  vector InitDist(real N0, real mu, real sigma, real PMR) {
    real sumP;
    vector[48] p;
    vector[48] ls;

    for (age in 1:48)
      p[age] = (PDF_normal(mu, sigma, age - 48) / PMR)
               + PDF_normal(mu, sigma, age)
               + (PDF_normal(mu, sigma, age + 48) * PMR);

    sumP = sum(p);

    for (age in 1:48)
      ls[age] = pow(10, N0) * p[age] / sumP;

    return ls;
  }

  // parasite aging with immunity
  vector Aging(vector ls, real PMR, real Q, real lag, real k, int t) {
    vector[48] tmpvec;

    if (t <= lag)
      tmpvec[1] = ls[48] * PMR;
    else
      tmpvec[1] = ls[48] * PMR * exp(-k * (t - lag));

    for (age in 2:48)
      tmpvec[age] = (1.0 - Q) * ls[age - 1];

    return tmpvec;
  }

  // qPCR observation model:
  // rho = fraction of sequestered parasites that contribute to qPCR
  // real CountPCR(vector ls, real rho) {
  //   return sum(ls[1:26]) + rho * sum(ls[27:48]);
  // }
  real CountPCR(vector ls) {
    return sum(ls[1:26]);
  }

  // qPCR version of NJWnQ
  // h0[n]  = lower hourly index (1-based)
  // frac[n] = fractional part of the observation time in hours
  vector NJWnQ_qPCR(
      real N0,
      real mu,
      real sigma,
      real PMR,
      real Q,
      real lag,
      real k,
      // real rho,
      array[] int h0,
      vector frac,
      int nsteps
  ) {
    int n_obs = size(h0);
    vector[n_obs] out;
    vector[nsteps + 1] val;
    vector[48] ls;

    ls = InitDist(N0, mu, sigma, PMR);

    // time 0
    // val[1] = log10(CountPCR(ls, rho));
    val[1] = log10(CountPCR(ls));

    // hourly latent simulation
    for (t in 1:nsteps) {
      ls = Aging(ls, PMR, Q, lag, k, t);
      // val[t + 1] = log10(CountPCR(ls, rho));
      val[t + 1] = log10(CountPCR(ls));
    }

    // interpolate to exact qPCR observation times
    for (n in 1:n_obs) {
      if (frac[n] == 0)
        out[n] = val[h0[n]];
      else
        out[n] = (1 - frac[n]) * val[h0[n]] + frac[n] * val[h0[n] + 1];
    }

    return out;
  }

} // end functions

data {
  int<lower=1> n_obs;
  vector[n_obs] obs_time_days;   // e.g. 1, 1.5, 2, 4.3, 6.3
  vector[n_obs] log10_qpcr;      // qPCR values already on log10 scale
}

transformed data {
  vector[n_obs] obs_time_hours;
  array[n_obs] int h0;
  vector[n_obs] frac;
  int nsteps_qpcr;

  for (n in 1:n_obs) {
    obs_time_hours[n] = 24.0 * obs_time_days[n];

    // lower hourly index, 1-based
    // example: 31.2 hours -> floor = 31 -> h0 = 32
    h0[n] = to_int(floor(obs_time_hours[n])) + 1;

    // fractional part within the hour
    frac[n] = obs_time_hours[n] - floor(obs_time_hours[n]);
  }

  // simulate through the last needed hourly step
  nsteps_qpcr = to_int(ceil(max(obs_time_hours)));
}

parameters {
  real N0;
  real mu;
  real<lower=0> sigma;
  real<lower=0> PMR;
  real<lower=0> Q;
  real<lower=0> lag;
  real<lower=0> k;
  // real<lower=0,upper=1> rho;
  real<lower=0> sigma_obs;
}

model {
  vector[n_obs] njwmodel;

  // priors
  N0 ~ uniform(1, 4);
  mu ~ uniform(1, 48);
  sigma ~ uniform(1, 8);
  PMR ~ uniform(1, 30);
  Q ~ uniform(0, 1);
  lag ~ uniform(0, 200);
  k ~ uniform(0, 2);
  // rho ~ beta(2, 2);
  sigma_obs ~ uniform(0.05, 1.5);

  njwmodel = NJWnQ_qPCR(
    N0, 
    mu, 
    sigma, 
    PMR, 
    Q, 
    lag, 
    k, 
    // rho, 
    h0,
    frac, 
    nsteps_qpcr);
  
  // likelihood
  log10_qpcr ~ normal(njwmodel, sigma_obs);
}

generated quantities {
  vector[n_obs] y_pred;
  vector[n_obs] log_lik;
  vector[n_obs] njwmodel;

  njwmodel = NJWnQ_qPCR(
    N0, 
    mu, 
    sigma, 
    PMR, 
    Q, 
    lag, 
    k, 
    // rho, 
    h0, 
    frac, 
    nsteps_qpcr);
  
  for (n in 1:n_obs) {
    y_pred[n] = normal_rng(njwmodel[n], sigma_obs);
    log_lik[n] = normal_lpdf(log10_qpcr[n] | njwmodel[n], sigma_obs);
  }
}