
clc; clear all; close all;

%% Parse AQWA output files
hydro = struct();
hydro = Read_AQWA(hydro,'AQWA_ANALYSIS.AH1','AQWA_ANALYSIS.LIS');

%% Calculate the IRF for radiation
% hydro = Radiation_IRF(hydro,t_end,n_t,n_w,w_min,w_max), [] for defaults
hydro = Radiation_IRF(hydro ,100 ,1001 ,1001 ,[] ,[]);

%% Calculate the approximate state space IRF for radiation, if using state space
% hydro = Radiation_IRF_SS(hydro,order_max,R2_min), [] for defaults
hydro.order_max = 3;
hydro = Radiation_IRF_SS(hydro ,hydro.order_max ,0.99);

%% Calculate the IRF for excitation
% hydro = Excitation_IRF(hydro,t_end,n_t,n_w,w_min,w_max), [] for defaults
hydro = Excitation_IRF(hydro ,100 ,1001 ,1001 ,[] ,[]);

%% Calculate the IRF for diffraction
hydro = Diffraction_IRF(hydro ,100 ,1001 ,1001 ,[] ,[]);


%% Plot a few key BEM results *optional*
Plot_BEMIO(hydro)

%% Generate file BEM_data.m:
save('../BEM_data.mat', 'hydro');



%% hydro structure 
%(See http://wec-sim.github.io/WEC-Sim/master/user/advanced_features.html#user-advanced-features-bemio)


% A       : [sum(dof),sum(dof),Nf]           : added mass
% Ainf    : [sum(dof),sum(dof)]              : infinite frequency added mass
% B       : [sum(dof),sum(dof),Nf]           : radiation damping
% beta    : [1,Nh]                           : wave headings (deg)
% body    : {1,Nb}                           : body names
% C       : [6,6,Nb]                         : hydrostatic restoring stiffness
% cb      : [3,Nb]                           : center of buoyancy
% cg      : [3,Nb]                           : center of gravity
% code    : string                           : BEM code (WAMIT, AQWA, or NEMOH)
% dof     : [1,Nb]                           : degrees of freedom for each body 
% ex_im   : [sum(dof),Nh,Nf]                 : imaginary component of excitation
% ex_K    : [sum(dof),Nh,length(ex_t)]       : excitation IRF
% ex_ma   : [sum(dof),Nh,Nf]                 : magnitude of excitation force
% ex_ph   : [sum(dof),Nh,Nf]                 : phase of excitation force
% ex_re   : [sum(dof),Nh,Nf]                 : real component of excitation
% ex_t    : [1,length(ex_t)]                 : time steps in the excitation IRF
% ex_w    : [1,length(ex_w)]                 : frequency step in the excitation IRF
% file    : string                           : BEM output filename
% g       : [1,1]                            : gravity
% h       : [1,1]                            : water depth
% Nb      : [1,1]                            : number of bodies
% Nf      : [1,1]                            : number of wave frequencies
% Nh      : [1,1]                            : number of wave headings
% ra_K    : [sum(dof),sum(dof),length(ra_t)] : radiation IRF
% ra_t    : [1,length(ra_t)]                 : time steps in the radiation IRF
% ra_w    : [1,length(ra_w)]                 : frequency steps in the radiation IRF
% rho     : [1,1]                            : density
% ss_A    : [sum(dof),sum(dof),ss_O,ss_O]    : state space A matrix
% ss_B    : [sum(dof),sum(dof),ss_O,1]       : state space B matrix
% ss_C    : [sum(dof),sum(dof),1,ss_O]       : state space C matrix
% ss_conv : [sum(dof),sum(dof)]              : state space convergence flag
% ss_D    : [sum(dof),sum(dof),1]            : state space D matrix
% ss_K    : [sum(dof),sum(dof),length(ra_t)] : state space radiation IRF
% ss_O    : [sum(dof),sum(dof)]              : state space order
% ss_R2   : [sum(dof),sum(dof)]              : state space R2 fit
% T       : [1,Nf]                           : wave periods
% Vo      : [1,Nb]                           : displaced volume
% w       : [1,Nf]                           : wave frequencies
