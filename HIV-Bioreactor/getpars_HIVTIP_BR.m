%% Neha Khetan
% Define parameters for Bioreactor simulations
%% 

function PP = getpars_HIVTIP_BR()
    PP        = {};
    PP.lam    = 0;  %50*10^3; %20*10^3
    PP.h      = 2*10^6;
    PP.dr     = 0;
    PP.drV    = 0;   
    PP.d      = 0.02;
    PP.k      = 2.5*10^-8;
    PP.d2     = 0.7;
    PP.c      = 6;    
    PP.n      = 200;
    PP.D      = 0.8; 
    PP.P      = 30;   
    PP.pr     = 10^-3;
    PP.h0     = 1;
end
