clear ; clf; clc;  close all
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neha Khetan, Sep. 2022 - 2023
% Aim: 
% 1. Implement and reproduce: HIV-TIP model ( Weinberger et al 2003 ) 
% 2. Parameter scans and sensitivity analysis
% 3. Modifications: 
%      i. death rate independent of "D"
%     ii. Hill functions for infectivity and interference
%    iii. Heterozygous formation HIV-TIP 
%     iv. Immature virion population and % immature ~ f( infectiivty )
%      v. super-infection [ m =2 , m =3 ]  
%     vi. CD8-mediated immune response
%    vii. Refractory compartment
%   viii. Latent population and secondary reserviours
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




TotalTime = 300;  % in days
% Initial values: cells/uL
T0    =  800;
I0    =  60;
V0    =  1; 
Vt0   =  1; 
It0   =  0;
Id0   =  0;
tval  = [ 0 :0.1:TotalTime ];
ivE   = [ T0 ; I0 ; V0 ; It0  ; Id0 ; Vt0  ];


pvaal     = getpars_HIVTIP( );
tolval    = 1e-6;
options1  = odeset('AbsTol', tolval  , 'RelTol', tolval  ,'Stats','off', 'InitialStep' , 10^-12  );
sol       = ode23s( @( tval , yy2 ) expandedHiv( tval , yy2 , pvaal  )  , tval, ivE );


figure(1),...
    plot( sol.x , sol.y( 1 ,: ) ,'linewidth' , 3 ), hold on,...
    plot( sol.x , sol.y( 4 ,: ) ,'linewidth' , 3 ), hold on,...
    plot( sol.x , ( sol.y( 1 ,: ) + sol.y( 4 ,: )) ,'linewidth' ,2 ), hold on,...
    set( gca , 'fontsize' , 18 ),...
    xlabel('Time (days)'),...
    ylabel('Healthy CD4+ (cells/\mul)'),...
    legend( 'Uninfected', 'TIP-infected' , 'Total healthy T cells')

figure(2),...
    plot( sol.x , sol.y( 2 ,: ) ,'linewidth' , 3 ), hold on,...
    plot( sol.x , sol.y( 5 ,: ) ,'linewidth' , 3 ), hold on,...
    plot( sol.x , ( sol.y( 2 ,: ) + sol.y( 5 ,: )) ,'linewidth' ,2 ), hold on,...
    set( gca , 'fontsize' , 18 ),...
    xlabel('Time (days)'),...
    ylabel('Unhealthy CD4+ (cells/\mul)'),...  
    legend(  'HIV-Infected' , 'TIP + HIV infected' , 'Total infected T cells')

figure(3),...
    plot( sol.x , sol.y( 3 ,: ).*10^3 , 'k' , 'linewidth' , 3 ), hold on,...
    plot( sol.x , sol.y( 6 ,: ).*10^3 , 'r' , 'linewidth' , 3 ), hold on,...
    set( gca , 'yscale' , 'log'),...
    set( gca , 'fontsize' , 18 ),...
    xlabel('Time (days)'),...
    ylabel('Viral load ( virions/mL)'),...  
    legend( 'HIV' , 'HIV-TIP')







