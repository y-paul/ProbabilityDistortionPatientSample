functions {
  
  real pwf(real p, real gam) {
    
    real wp;
    
    if (p > 0) { wp = exp(-pow(-log(p), gam)); } else { wp = 0; }
    
    return wp;
  }
}

data {
  // ---------- Patient group ----------
  int<lower=0> n_pat;                       // trials
  int<lower=1> N_pat;                       // subjects
  int<lower=1, upper=N_pat> sub_pat[n_pat]; // subject index

  matrix[n_pat, 2] xa_pat;
  matrix[n_pat, 2] pa_pat;
  matrix[n_pat, 2] xb_pat;
  matrix[n_pat, 2] pb_pat;

  matrix<lower=0, upper=1>[n_pat, 2] fa_pat;
  matrix<lower=0, upper=1>[n_pat, 2] fb_pat;

  // choices: col1 = description, col2 = experience
  int<lower=0, upper=1> choices_pat[n_pat, 2];

  // ---------- Control group ----------
  int<lower=0> n_con;
  int<lower=1> N_con;
  int<lower=1, upper=N_con> sub_con[n_con];

  matrix[n_con, 2] xa_con;
  matrix[n_con, 2] pa_con;
  matrix[n_con, 2] xb_con;
  matrix[n_con, 2] pb_con;

  matrix<lower=0, upper=1>[n_con, 2] fa_con;
  matrix<lower=0, upper=1>[n_con, 2] fb_con;

  int<lower=0, upper=1> choices_con[n_con, 2];
}

parameters {
  // ---------- patient: individual z-scores ----------
  matrix[2, N_pat] lam_z_pat;
  matrix[2, N_pat] gam_z_pat;
  matrix[2, N_pat] theta_z_pat;

  // pop means (probit scale) per condition: [1]=experience, [2]=description
  row_vector[2] lam_mu_phi_pat;
  row_vector[2] gam_mu_phi_pat;
  row_vector[2] theta_mu_phi_pat;

  // scales and correlations (over conditions)
  vector<lower=0>[2] lam_sigma_pat;
  vector<lower=0>[2] gam_sigma_pat;
  vector<lower=0>[2] theta_sigma_pat;

  cholesky_factor_corr[2] L_lam_omega_pat;
  cholesky_factor_corr[2] L_gam_omega_pat;
  cholesky_factor_corr[2] L_theta_omega_pat;

  // ---------- control: individual z-scores ----------
  matrix[2, N_con] lam_z_con;
  matrix[2, N_con] gam_z_con;
  matrix[2, N_con] theta_z_con;

  // pop means (probit scale) per condition
  row_vector[2] lam_mu_phi_con;
  row_vector[2] gam_mu_phi_con;
  row_vector[2] theta_mu_phi_con;

  // scales and correlations (over conditions)
  vector<lower=0>[2] lam_sigma_con;
  vector<lower=0>[2] gam_sigma_con;
  vector<lower=0>[2] theta_sigma_con;

  cholesky_factor_corr[2] L_lam_omega_con;
  cholesky_factor_corr[2] L_gam_omega_con;
  cholesky_factor_corr[2] L_theta_omega_con;
}

transformed parameters {
  // ---------- patient: individual-level probit parameters ----------
  matrix[N_pat, 2] lam_phi_pat =
    (diag_pre_multiply(lam_sigma_pat, L_lam_omega_pat) * lam_z_pat)';
  matrix[N_pat, 2] gam_phi_pat =
    (diag_pre_multiply(gam_sigma_pat, L_gam_omega_pat) * gam_z_pat)';
  matrix[N_pat, 2] theta_phi_pat =
    (diag_pre_multiply(theta_sigma_pat, L_theta_omega_pat) * theta_z_pat)';

  vector[N_pat] lambda_ex_pat = Phi(lam_mu_phi_pat[1] + lam_phi_pat[,1]) * 5;
  vector[N_pat] lambda_de_pat = Phi(lam_mu_phi_pat[2] + lam_phi_pat[,2]) * 5;
  vector[N_pat] gamma_ex_pat  = Phi(gam_mu_phi_pat[1] + gam_phi_pat[,1]) * 5;
  vector[N_pat] gamma_de_pat  = Phi(gam_mu_phi_pat[2] + gam_phi_pat[,2]) * 5;
  vector[N_pat] theta_ex_pat  = Phi(theta_mu_phi_pat[1] + theta_phi_pat[,1]) * 5;
  vector[N_pat] theta_de_pat  = Phi(theta_mu_phi_pat[2] + theta_phi_pat[,2]) * 5;

  // trial-level SV differences
  vector[n_pat] de_sv_diff_pat;
  vector[n_pat] ex_sv_diff_pat;

  for (i in 1:n_pat) {
    int s = sub_pat[i];

    // DESCRIPTION (choices_pat[,1])
    {
      real vA1 = xa_pat[i,1];
      real vA2 = - lambda_de_pat[s] * xa_pat[i,2];
      real wA1 = pwf(pa_pat[i,1], gamma_de_pat[s]);
      real wA2 = pwf(pa_pat[i,2], gamma_de_pat[s]);
      real svA = vA1 * wA1 + vA2 * wA2;

      real vB1 = xb_pat[i,1];
      real vB2 = - lambda_de_pat[s] * xb_pat[i,2];
      real wB1 = pwf(pb_pat[i,1], gamma_de_pat[s]);
      real wB2 = pwf(pb_pat[i,2], gamma_de_pat[s]);
      real svB = vB1 * wB1 + vB2 * wB2;

      de_sv_diff_pat[i] = theta_de_pat[s] * (svA - svB);
    }

    // EXPERIENCE (choices_pat[,2])
    {
      real vA1 = xa_pat[i,1];
      real vA2 = - lambda_ex_pat[s] * xa_pat[i,2];
      real wA1 = pwf(fa_pat[i,1], gamma_ex_pat[s]);
      real wA2 = pwf(fa_pat[i,2], gamma_ex_pat[s]);
      real svA = vA1 * wA1 + vA2 * wA2;

      real vB1 = xb_pat[i,1];
      real vB2 = - lambda_ex_pat[s] * xb_pat[i,2];
      real wB1 = pwf(fb_pat[i,1], gamma_ex_pat[s]);
      real wB2 = pwf(fb_pat[i,2], gamma_ex_pat[s]);
      real svB = vB1 * wB1 + vB2 * wB2;

      ex_sv_diff_pat[i] = theta_ex_pat[s] * (svA - svB);
    }
  }

  // ---------- control: individual-level probit parameters ----------
  matrix[N_con, 2] lam_phi_con =
    (diag_pre_multiply(lam_sigma_con, L_lam_omega_con) * lam_z_con)';
  matrix[N_con, 2] gam_phi_con =
    (diag_pre_multiply(gam_sigma_con, L_gam_omega_con) * gam_z_con)';
  matrix[N_con, 2] theta_phi_con =
    (diag_pre_multiply(theta_sigma_con, L_theta_omega_con) * theta_z_con)';

  vector[N_con] lambda_ex_con = Phi(lam_mu_phi_con[1] + lam_phi_con[,1]) * 5;
  vector[N_con] lambda_de_con = Phi(lam_mu_phi_con[2] + lam_phi_con[,2]) * 5;
  vector[N_con] gamma_ex_con  = Phi(gam_mu_phi_con[1] + gam_phi_con[,1]) * 5;
  vector[N_con] gamma_de_con  = Phi(gam_mu_phi_con[2] + gam_phi_con[,2]) * 5;
  vector[N_con] theta_ex_con  = Phi(theta_mu_phi_con[1] + theta_phi_con[,1]) * 5;
  vector[N_con] theta_de_con  = Phi(theta_mu_phi_con[2] + theta_phi_con[,2]) * 5;

  vector[n_con] de_sv_diff_con;
  vector[n_con] ex_sv_diff_con;

  for (i in 1:n_con) {
    int s = sub_con[i];

    // DESCRIPTION (choices_con[,1])
    {
      real vA1 = xa_con[i,1];
      real vA2 = - lambda_de_con[s] * xa_con[i,2];
      real wA1 = pwf(pa_con[i,1], gamma_de_con[s]);
      real wA2 = pwf(pa_con[i,2], gamma_de_con[s]);
      real svA = vA1 * wA1 + vA2 * wA2;

      real vB1 = xb_con[i,1];
      real vB2 = - lambda_de_con[s] * xb_con[i,2];
      real wB1 = pwf(pb_con[i,1], gamma_de_con[s]);
      real wB2 = pwf(pb_con[i,2], gamma_de_con[s]);
      real svB = vB1 * wB1 + vB2 * wB2;

      de_sv_diff_con[i] = theta_de_con[s] * (svA - svB);
    }

    // EXPERIENCE (choices_con[,2])
    {
      real vA1 = xa_con[i,1];
      real vA2 = - lambda_ex_con[s] * xa_con[i,2];
      real wA1 = pwf(fa_con[i,1], gamma_ex_con[s]);
      real wA2 = pwf(fa_con[i,2], gamma_ex_con[s]);
      real svA = vA1 * wA1 + vA2 * wA2;

      real vB1 = xb_con[i,1];
      real vB2 = - lambda_ex_con[s] * xb_con[i,2];
      real wB1 = pwf(fb_con[i,1], gamma_ex_con[s]);
      real wB2 = pwf(fb_con[i,2], gamma_ex_con[s]);
      real svB = vB1 * wB1 + vB2 * wB2;

      ex_sv_diff_con[i] = theta_ex_con[s] * (svA - svB);
    }
  }
}

model {
  // ---------- patient priors ----------
  gam_mu_phi_pat   ~ std_normal();
  to_vector(gam_z_pat) ~ std_normal();
  gam_sigma_pat    ~ normal(0.5, 0.13);
  L_gam_omega_pat  ~ lkj_corr_cholesky(4);

  lam_mu_phi_pat   ~ std_normal();
  to_vector(lam_z_pat) ~ std_normal();
  lam_sigma_pat    ~ normal(0.5, 0.13);
  L_lam_omega_pat  ~ lkj_corr_cholesky(4);

  theta_mu_phi_pat   ~ std_normal();
  to_vector(theta_z_pat) ~ std_normal();
  theta_sigma_pat    ~ normal(0.5, 0.13);
  L_theta_omega_pat  ~ lkj_corr_cholesky(4);

  // ---------- control priors ----------
  gam_mu_phi_con   ~ std_normal();
  to_vector(gam_z_con) ~ std_normal();
  gam_sigma_con    ~ normal(0.5, 0.13);
  L_gam_omega_con  ~ lkj_corr_cholesky(4);

  lam_mu_phi_con   ~ std_normal();
  to_vector(lam_z_con) ~ std_normal();
  lam_sigma_con    ~ normal(0.5, 0.13);
  L_lam_omega_con  ~ lkj_corr_cholesky(4);

  theta_mu_phi_con   ~ std_normal();
  to_vector(theta_z_con) ~ std_normal();
  theta_sigma_con    ~ normal(0.5, 0.13);
  L_theta_omega_con  ~ lkj_corr_cholesky(4);

  // ---------- likelihood ----------
  // same convention as original model:
  // col1 = description, col2 = experience
  choices_pat[,1] ~ bernoulli_logit(de_sv_diff_pat);
  choices_pat[,2] ~ bernoulli_logit(ex_sv_diff_pat);

  choices_con[,1] ~ bernoulli_logit(de_sv_diff_con);
  choices_con[,2] ~ bernoulli_logit(ex_sv_diff_con);
}

generated quantities {
  // ---------- patient: population means (0..5) ----------
  real mu_gamma_ex_pat  = Phi(gam_mu_phi_pat[1]) * 5;
  real mu_gamma_de_pat  = Phi(gam_mu_phi_pat[2]) * 5;
  real mu_lambda_ex_pat = Phi(lam_mu_phi_pat[1]) * 5;
  real mu_lambda_de_pat = Phi(lam_mu_phi_pat[2]) * 5;
  real mu_theta_ex_pat  = Phi(theta_mu_phi_pat[1]) * 5;
  real mu_theta_de_pat  = Phi(theta_mu_phi_pat[2]) * 5;

  // ---------- control: population means (0..5) ----------
  real mu_gamma_ex_con  = Phi(gam_mu_phi_con[1]) * 5;
  real mu_gamma_de_con  = Phi(gam_mu_phi_con[2]) * 5;
  real mu_lambda_ex_con = Phi(lam_mu_phi_con[1]) * 5;
  real mu_lambda_de_con = Phi(lam_mu_phi_con[2]) * 5;
  real mu_theta_ex_con  = Phi(theta_mu_phi_con[1]) * 5;
  real mu_theta_de_con  = Phi(theta_mu_phi_con[2]) * 5;

  // ---------- sigmas (for monitoring) ----------
  real sig_gamma_ex_pat = gam_sigma_pat[1];
  real sig_gamma_de_pat = gam_sigma_pat[2];
  real sig_lambda_ex_pat = lam_sigma_pat[1];
  real sig_lambda_de_pat = lam_sigma_pat[2];
  real sig_theta_ex_pat = theta_sigma_pat[1];
  real sig_theta_de_pat = theta_sigma_pat[2];

  real sig_gamma_ex_con = gam_sigma_con[1];
  real sig_gamma_de_con = gam_sigma_con[2];
  real sig_lambda_ex_con = lam_sigma_con[1];
  real sig_lambda_de_con = lam_sigma_con[2];
  real sig_theta_ex_con = theta_sigma_con[1];
  real sig_theta_de_con = theta_sigma_con[2];

  // ---------- correlations between conditions (per group) ----------
  matrix[2,2] gam_omega_pat   = L_gam_omega_pat * L_gam_omega_pat';
  matrix[2,2] lam_omega_pat   = L_lam_omega_pat * L_lam_omega_pat';
  matrix[2,2] theta_omega_pat = L_theta_omega_pat * L_theta_omega_pat';

  real gam_r_pat   = gam_omega_pat[1,2];
  real lam_r_pat   = lam_omega_pat[1,2];
  real theta_r_pat = theta_omega_pat[1,2];

  matrix[2,2] gam_omega_con   = L_gam_omega_con * L_gam_omega_con';
  matrix[2,2] lam_omega_con   = L_lam_omega_con * L_lam_omega_con';
  matrix[2,2] theta_omega_con = L_theta_omega_con * L_theta_omega_con';

  real gam_r_con   = gam_omega_con[1,2];
  real lam_r_con   = lam_omega_con[1,2];
  real theta_r_con = theta_omega_con[1,2];

  // ---------- log-likelihoods for LOO ----------
  real log_lik_de_pat[n_pat];
  real log_lik_ex_pat[n_pat];
  real log_lik_de_con[n_con];
  real log_lik_ex_con[n_con];

  for (i in 1:n_pat) {
    log_lik_de_pat[i] = bernoulli_logit_lpmf(choices_pat[i,1] | de_sv_diff_pat[i]);
    log_lik_ex_pat[i] = bernoulli_logit_lpmf(choices_pat[i,2] | ex_sv_diff_pat[i]);
  }
  for (i in 1:n_con) {
    log_lik_de_con[i] = bernoulli_logit_lpmf(choices_con[i,1] | de_sv_diff_con[i]);
    log_lik_ex_con[i] = bernoulli_logit_lpmf(choices_con[i,2] | ex_sv_diff_con[i]);
  }
  real de_gap_gamma_pat = mu_gamma_de_pat - mu_gamma_ex_pat;
  real de_gap_gamma_con = mu_gamma_de_con - mu_gamma_ex_con;
  real de_gap_lambda_pat = mu_lambda_de_pat - mu_lambda_ex_pat;
  real de_gap_lambda_con = mu_lambda_de_con - mu_lambda_ex_con;
  real de_gap_theta_pat = mu_theta_de_pat - mu_theta_ex_pat;
  real de_gap_theta_con = mu_theta_de_con - mu_theta_ex_con;

  // ---------- Group differences in DE gaps ----------
  real group_diff_gamma_gap = de_gap_gamma_pat - de_gap_gamma_con;
  real group_diff_lambda_gap = de_gap_lambda_pat - de_gap_lambda_con;
  real group_diff_theta_gap = de_gap_theta_pat - de_gap_theta_con;
}