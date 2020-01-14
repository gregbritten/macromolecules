functions {
  real[] macro(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

    real KN    = theta[1];    
    real mu    = theta[2]; 
    real CHsyn = theta[3]; 
    real m_ex  = theta[4];  
    real R_ex  = theta[5];  
    real tau   = theta[6];
    real b     = theta[7];
    
    real CH = fmax(x[1],0.0);
    real PR = fmax(x[2],0.0);
    real Chl = fmax(x[3],0.0);
    real N = fmax(x[4],0.0);

    real PRsynth = mu*N/(KN+ pow(KN*N,0.5) + N);
    real r0      = theta[7]*(PR/CH);
    //real Chl     = Chl*P;
    real Rcell   = x[1]/x[2];
    real excr    = 0.5*theta[4]*(1+tanh(Rcell - theta[5]));
    
    real dCH    = x[2]*(theta[3] - excr);
    real dr     = (1/theta[6])*(r0-x[3]);
    real dPR    = x[2]*PRsynth;
    real dN     = -dPR/(1+exp(-10000*x[4])); 

    return {dCH,dPR,dr,dN};
  }
}
data {
  int<lower = 0> n;           // num obs
  real t_obs[n];              // obs times
  real<lower = 0> y[n,4];     // observed variable at measurement times
}
parameters {
  real<lower = 0> theta[7];   // parameters
  real<lower = 0> x0[4];      // initial population
  real<lower = 1E-15> sigma[4]; // obs error
}
transformed parameters {
  real x[n,4] = integrate_ode_rk45(macro, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0), 1e-6, 1e-5, 1e3) ;
  for(i in 1:n){
    x[i,3] = x[i,3]*x[i,2]*1E6;
  }
}
model {
  x0[1]    ~ normal(0.1,1);
  x0[2]    ~ normal(0.1,1);
  x0[3]    ~ normal(10,10);
  x0[4]    ~ normal(0.1,1);  
  theta[1] ~ normal(0.002,3);
  theta[2] ~ normal(0.3,1);
  theta[3] ~ normal(5,10);
  theta[4] ~ normal(10,10);
  theta[5] ~ normal(13,10);
  theta[6] ~ normal(10,10); 
  theta[7] ~ normal(0.05,5); 
  
  for(i in 1:4){
    y[1:n,i] ~ normal(x[1:n,i], sigma[i]);
  }
}

