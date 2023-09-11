
function PP = getpars_HIVTIP_BR()
    PP        = {};
    PP.lam    = 0;%50*10^3; %20*10^3;%50*10^3;%40*10^3;


    PP.h      = 2*10^6;
    

    PP.dr     = 0;%0.1;% 0.02;
    PP.drV    = 0; %0.250;%0.025;



    
    PP.d      = 0.02;
    PP.k      = 2.5*10^-8;
    PP.d2     = 0.7;
    PP.c      = 6;
    
    PP.n      = 200;
    PP.D      = 0.8; %0.8       %0.9; %psi
    PP.P      = 30; %100;        %1;   %rho
    
    
    PP.pr     = 10^-3;
    PP.h0     =1;

end
