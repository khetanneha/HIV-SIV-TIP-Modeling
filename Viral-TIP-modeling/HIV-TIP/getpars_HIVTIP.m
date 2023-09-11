%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Parameters for HIV-TIP models
%  Most from Weinberger et al. 2003
%  viral cleaeance : Ramratnam B, Bonhoeffer S, Binley J, Hurley A, Zhang L, Mittler JE, Markowitz M, Moore JP, 
%                    Perelson AS, Ho DD. Rapid production and clearance of HIV-1 and hepatitis C virus assessed by large volume plasma apheresis. Lancet. 1999 
%  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PP = getpars_HIVTIP()
    PP        = {};
    PP.lam    = 30;
    PP.d      = 0.02;
    PP.k      = 1.8*10^-4;
    PP.d2     = 0.7;
    PP.c      = 23;
    
    PP.n      = 200;
    PP.D      = 0.4;
    PP.P      = 25; 
end