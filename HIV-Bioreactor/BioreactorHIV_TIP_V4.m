%% ========================================================================
% Neha Khetan, June 2023
% HIV-TIP - Bioreactor modeling for design criterion for amplification
% of spontaneously emerged DIPs in an evolutionary experimental scheme
% Design criterion, conditions: This is used for TIP-HIV paper
% Three different Experimental Scenarios hypothesized are simulated in-silico 
% to predict the conditions that favor amplification , and establishment of DIPs
% and the time-scales of these ...........................................
%% Experimenta set-up
% BRsetup: 0 , 1 , 2
% Case 0: Model 1: Serial-passage culturing
% Case 1: Model 2: Non-dilutive reactor culturing of HIV 
% Case 2: Model 2: Target-cell replenishment in reactor
% ---------------------------------------------------------------------
clear global; clc ; clf; close all;

BRsetup = 1;
switch ( BRsetup )
    case 0
        opath     ='./Revised/SP/';
        tp        = 25 ;
        dilFac    = 0.1;
        iEnd      = 100;
    case 1
        opath     ='./Revised/Cont/';
        tp        = 150;
        iEnd      = 150;
    case 2
        opath     ='./Revised/ContCD/';
        tp        = 300;
        iEnd      = 300;
end

timTot  = iEnd;
timstep  = 1;
%Initial condition:    T, I , V , Tt , Td , Vt
ivE0     = [ 1.6*10^6 , 0  , 2*10^4  , 0 , 0 , 0 , 0 ];
ivE      = ivE0;
iStart   = 0;



%% =====================================================================
%
%% =====================================================================
set(0,'defaulttextfontsize',16);
set(0,'defaultaxesfontsize',16);
set( gcf ,'color','w');
tolval    = 1e-12;
abtol     = tolval;
options1  = odeset('AbsTol', tolval  , 'RelTol', tolval );%,'Stats','off', 'InitialStep' , 10^-12  );
Pv          = 0.7;   % [ 0.5 , 0.7 , 1 , 1.5 , 5 ,10  ];%,  5 , 15 , 45  , 100  ];
Dv          = 0.9; %[   0.1 , 0.8 ]; %[ 0.04 , 0.4 , 0.8 ];%  0.08 , 0.2 ,  0.8 , 1  ];
LowCutOff   = 10^-3;%10^-12;


PDpairs = [];
for i = 1:length(Pv)
    for k =1:length(Dv)
        PDpairs = [PDpairs ; Pv(i) , Dv(k) ];
    end
end



dilFac = [  10 , 5 , 3, 2 , 1 ];
for pd = 1:size( dilFac ,2  )
    % Iterating over dilution-factors
    AllStates = [];
    sol = [];
    for iTime = iStart:tp:iEnd-tp
        % tp: passaging time
        Timvals         = [ iTime iTime+tp ];
        timpts          =  [ Timvals(1):timstep:Timvals(2) ];
        pvaal           = getpars_HIVTIP_BR();
        pvaal.D         = Dv;
        pvaal.P         = Pv;
        sol             = ode23s( @( Timvals , yy2 )TipmodelBioreactor( Timvals , yy2 , pvaal  )  , Timvals , ivE  , options1 );
        sol2            = deval( sol ,  timpts  );
        tmpIdx          = sol2(:,:) <LowCutOff ;
        sol2(tmpIdx )   = LowCutOff ;
        AllStates       = [ AllStates , sol2 ];     
        if BRsetup == 0
            ivE            = [ ivE0(1:2)  , sol.y(3,end)/dilFac(pd) , ivE0(4:5) , sol.y(6,end)/dilFac(pd) ,  ivE0(7) ];
        else
            ivE            = [ ivE0(1:7)  ];
        end

    end

    AllStatesIdx   = AllStates(:,:) < LowCutOff;
    AllStates( AllStatesIdx ) = 0;
    TotalTCells               = ( AllStates(1,:) + AllStates( 2,:) + AllStates( 4,:) +AllStates( 5,:));
    TotalInf                  = ((AllStates(2,:) + AllStates(5,:))./TotalTCells).*100;
    Infect                    = ( AllStates(2,:)./TotalTCells).*100 ;
    InfectDual                = ( AllStates(5,:)./TotalTCells).*100 ;

    TotalDIP                  = (AllStates(4,:) + AllStates(5,:));

    %% Ext Fig1 : Viral load, DIP and Total T cells
    % figure(1),...
    %   yyaxis left,...
    %   plot( ((AllStates(4,:)+AllStates(5,:))./TotalTCells).*100 , '-',  'linewidth' ,1.5   ), hold on,...
    %   yyaxis right,...
    %   plot( ((AllStates(4,:)+AllStates(5,:))) , '-',  'linewidth' ,1.5   ), hold on,...
    % figure(2),...
    %     yyaxis left,...
    %     semilogy(  ( AllStates(6,:)./( AllStates(3,:) + AllStates(6,:) )).*100   ,  'linewidth' ,1 ), hold on,...
    %     ylabel(' % DIPs/[DIPs + HIV] '),...
    %     yyaxis right,...
    %     semilogy(  (AllStates(6,:)  ) ,  'linewidth' ,1 ), hold on,...
    %     ylabel('# DIPs/mL'),...

    figure(10),...
        ax = gca();
         yyaxis left,...
         plot(  timpts,     TotalTCells./10^6 , '-'   ,  'linewidth' ,1.5    , 'color' , [0.7 , 0.7 ,0.7 ] ), hold on,...
         plot(  timpts,     TotalDIP./10^6 , '-.'   ,  'linewidth' ,1.5    , 'color' , [0.7 , 0.7 ,0.7 ] ), hold on,...
         yyaxis right,...
         semilogy(  timpts,  AllStates(3,:)  , '-',  'linewidth' ,1.5 ), hold on,...
         semilogy( timpts,  AllStates(6,:)  , '-.', 'linewidth' ,1.5 ), hold on;
         ax.YAxis(1).Color = [ 0.7 , 0.7 , 0.7] ;
         ax.YAxis(2).Color = [ 1 0.25 0.25] ; 
end



