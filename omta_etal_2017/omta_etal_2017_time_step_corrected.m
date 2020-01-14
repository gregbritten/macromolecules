%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CNpro=6.6;       %CN ratio of protein
KN   =0.002;      %half-saturation constant  
mu   =0.234;     %max protein synthesis rate
CHsyn=2.19;      %max photosynthetic rate
m_ex =1.43;      %maximum excretion rate
R_ex =13.49;     %CN above which excretion happens
tau  =9.59;      %photoacclimation time
b    =0.05705;   %change in Chl:N relative to N:C

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOGISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01;       %timestep
n=25/dt;         %number of solution points
CH=NaN(n,1);     %container for organic C
PR=NaN(n,1);     %container for organic N
theta=NaN(n,1);  %container for Chl/PR ratio
N=NaN(n,1);      %container for environmental N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIAL CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CH(1)    = 0.112;           %initial organic C
PR(1)    = 0.0169;           %initial organic N
theta(1) = (4.1/1E6)/0.0169; %initial Chl/PR ratio
N(1)     = 0.0911;           %initial environmental N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EULER FORWARD LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:n
    %PRsynth = mu*N(t-1)/(KN + sqrt(KN*N(t-1)) + N(t-1));             %protein synthesis
    PRsynth = mu*N(t-1)/(KN + N(t-1));             %protein synthesis
    theta0  = b*(PR(t-1)/CH(t-1));               %target Chl/PR ratio
    Chl     = theta(t-1)*PR(t-1);                %cellular Chl
    Rcell   = CH(t-1)/PR(t-1);                   %C:N ratio
    excr    = (1/2)*m_ex*(1+tanh(Rcell - R_ex)); %excretion of excess C
    
    dCH_dt    = dt*PR(t-1)*(CHsyn - excr);       %Euler CH increment
    dtheta_dt = dt*(1/tau)*(theta0-theta(t-1));  %theta increment
    dPR_dt    = dt*PR(t-1)*PRsynth;              %PR increment            
    dN_dt     = -dPR_dt/(1+exp(-10000*N(t-1))); %N increment 
    
    CH(t)    = CH(t-1) + dCH_dt;       %Euler CH update
    PR(t)    = PR(t-1) + dPR_dt;       %PR update
    theta(t) = theta(t-1) + dtheta_dt; %theta update
    N(t)     = N(t-1) + dN_dt;         %N update 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ts = linspace(5,25,n);
figure(2)
subplot(2,2,1); plot(ts,CH); title('Organic carbon (mM)');
subplot(2,2,2); plot(ts,PR); title('Organic nitrogen (mM)');
subplot(2,2,3); plot(ts,theta.*PR*1E6); title('Chl (nM)'); 
subplot(2,2,4); plot(ts,N); title('Environmental nitrogen (mM)')
