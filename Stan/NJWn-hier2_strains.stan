functions {
  // PDF of normal distribution
  real PDF_normal(real mu, real sig, real x) {
    return (1 / (exp(((x - mu) * (x - mu)) / (2 * sig * sig)) * sqrt(2 * pi()) * sig));
  }

  // for generating the initial parasite distribution
  vector InitDist(real N0, real mu, real sigma, real PMR) {
    real sumP;
    vector[48] p;
    vector[48] ls;

    for (age in 1:48) {
      p[age] = (PDF_normal(mu, sigma, age - 48) / PMR) + PDF_normal(mu, sigma, age) + (PDF_normal(mu, sigma, age + 48) * PMR);
    }

    sumP = sum(p);

    for (age in 1:48) {
      ls[age] = pow(10, N0) * p[age] / sumP;
    }

    return ls;
  }

  // parasite aging
  vector RotateRight(vector ls, real PMR) {
    vector[48] tmpvec;
    tmpvec[1] = ls[48] * PMR;
    for (age in 2:48) {
      tmpvec[age] = ls[age - 1];
    }
    return tmpvec;
  }

  // the main function of NJW model 
  vector NJWn(real N0, real mu, real sigma, real PMR, int npoints) {
    vector[48] ls;
    vector[npoints] vecout;
    int runs;
    int pointindx;

    runs = (npoints + 1) * 24;
    pointindx = 1;

    ls = InitDist(N0, mu, sigma, PMR);
    for (t in 1:runs) {
      // aging 
      ls = RotateRight(ls, PMR);

      // record the circulating parasites every 24 hrs
      if (t % 24 == 0 && pointindx <= npoints) {
        // circulating parasites = sum from age 1-26 hrs
        vecout[pointindx] = log10(sum(ls[1:26]));
        pointindx = pointindx + 1;
      }
    }
    return vecout;
  }

  // parasite aging function for NJWnQ
  vector Aging(vector ls, real PMR, real Q, real lag, real k, int t) {
    vector[48] tmpvec;

    if (t <= lag) {
      tmpvec[1] = ls[48] * PMR;
    } else {
      tmpvec[1] = ls[48] * PMR * exp(-k * (t - lag));
    }

    for (age in 2:48) {
      tmpvec[age] = (1.0 - Q) * ls[age - 1];
    }

    return tmpvec;
  }

  // the main model function (NJW + Immunity)
  vector NJWnQ(real N0, real mu, real sigma, real PMR, real Q, real lag, real k, int npoints) {
    vector[48] ls;
    vector[npoints] vecout;
    int runs;
    int pointindx;

    runs = (npoints + 1) * 24;
    pointindx = 1;

    ls = InitDist(N0, mu, sigma, PMR);
    for (t in 1:runs) {
      // aging
      ls = Aging(ls, PMR, Q, lag, k, t);

      // record the circulating parasites every 24 hrs
      if (t % 24 == 0 && pointindx <= npoints) {
        // circulating parasites = sum from age 1-26 hrs
        vecout[pointindx] = log10(sum(ls[1:26]));
        pointindx = pointindx + 1;
      }
    }
    return vecout;
  }

  // Partial sum function for reduce_sum for NJWn
  real partial_sum_reduce(array[] int slice, int start, int end, array[] real N0, array[] real mu, array[] real sigma, 
    array[] real PMR, array[] int N, matrix y, array[] real sig) {

    real lp = 0.0;
    matrix[end-start+1, 14] temp_njwmodel;

    for (m in start:end) {
      temp_njwmodel[m - start + 1, 1:N[m]] = to_row_vector(NJWn(N0[m], mu[m], sigma[m], PMR[m], N[m]));

      for (t in 1:N[m]) {
        lp += normal_lpdf(y[m, t] | temp_njwmodel[m - start + 1, t], sig[m]);
      }
    }
    return lp;
  }
}

data {
  int<lower=1> N_pid; // Number of datasets
  array[N_pid] int N_obs;
  array[N_pid] int N_mis;
  array[N_pid, 14] int ii_obs;
  array[N_pid, 4] int ii_mis;
  matrix[N_pid, 14] y_obs;  
  array[N_pid] int strain_indx; // strain id {"Santee_C" -> 1, "McLendon" -> 2, "El_Limon" -> 3}
}

transformed data {
  array[N_pid] int N;
  for (m in 1:N_pid) {
    N[m] = N_obs[m] + N_mis[m];
  }
}

parameters {
  // Parameters for each dataset
  array[N_pid] real<lower=1, upper=4> N0;
  array[N_pid] real<lower=1, upper=48> mu;
  array[N_pid] real<lower=1, upper=48> sigma;
  array[N_pid] real<lower=1, upper=50> PMR;
  array[N_pid] real<lower=0, upper=4> sig;

  matrix[N_pid, 14] y_mis;

  // Hyperparameters for population-level distributions
  array[3] real<lower=0, upper=4> N0_pop_mean;
  array[3] real<lower=0.001> N0_pop_sd;

  array[3] real<lower=1> mu_pop_mean;
  array[3] real<lower=0> mu_pop_sd;

  array[3] real<lower=1> sigma_pop_mean;
  array[3] real<lower=0.01> sigma_pop_sd;

  array[3] real<lower=0> PMR_pop_mean;
  array[3] real<lower=0> PMR_pop_sd;

}

transformed parameters {
  matrix[N_pid, 14] y;

  for (m in 1:N_pid) {
    y[m, ii_obs[m, 1:N_obs[m]]] = y_obs[m, 1:N_obs[m]];

    if (N_mis[m] > 0) {
      y[m, ii_mis[m, 1:N_mis[m]]] = y_mis[m, 1:N_mis[m]];
    }
  }
}

model {
  matrix[N_pid, 14] njwqmodel;

  // Priors for hyperparameters
  target += gamma_lpdf(N0_pop_mean | 10, 0.15);
  target += gamma_lpdf(N0_pop_sd | 2.5, 0.1);

  target += gamma_lpdf(mu_pop_mean | 5, 5);
  target += gamma_lpdf(mu_pop_sd | 5, 1);

  target += uniform_lpdf(sigma_pop_mean | 1, 48);
  target += gamma_lpdf(sigma_pop_sd | 5, 1.5);

  target += gamma_lpdf(PMR_pop_mean | 2.5, 1);
  target += gamma_lpdf(PMR_pop_sd | 5, 0.1);

  for (m in 1:N_pid) {
    // Priors for dataset-level parameters
    target += lognormal_lpdf(N0[m] | N0_pop_mean[strain_indx[m]], N0_pop_sd[strain_indx[m]]);
    target += normal_lpdf(mu[m] | mu_pop_mean[strain_indx[m]], mu_pop_sd[strain_indx[m]]);
    target += normal_lpdf(sigma[m] | sigma_pop_mean[strain_indx[m]], sigma_pop_sd[strain_indx[m]]);
    target += lognormal_lpdf(PMR[m] | PMR_pop_mean[strain_indx[m]], PMR_pop_sd[strain_indx[m]]);
    target += gamma_lpdf(sig[m] | 2, 1);

  }

  // Likelihood
  array[N_pid] int indx_slice = linspaced_int_array(N_pid, 1, N_pid);
  target += reduce_sum(partial_sum_reduce, indx_slice, 1, N0, mu, sigma, PMR, N, y, sig);
}

generated quantities {
  matrix[N_pid, 14] log_lik;
  matrix[N_pid, 14] njwmodel;
  matrix[N_pid, 14] y_pred;

  for (m in 1:N_pid) {
    njwmodel[m, 1:N[m]] = to_row_vector(NJWn(N0[m], mu[m], sigma[m], PMR[m], N[m]));

    for (n in 1:N[m]) {
      y_pred[m, n] = normal_rng(njwmodel[m, n], sig[m]);
      log_lik[m, n] = normal_lpdf(y[m, n] | njwmodel[m, n], sig[m]);
    }
  }
}
