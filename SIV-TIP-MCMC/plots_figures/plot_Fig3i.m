%% Neha Khetan
%  Draft Fig: peak:steady state and abs values from
%  22 May 2023, mathematica --> Matlab
%  Cleaning up: Fig. 4E for submission




clearvars; clc; close all


ipath   = './out-files/'
LHSdataTIP = importdata([  ipath , 'LHS_Values_TIP_M1.out' ]); 
LHSdataBasic = importdata( [ ipath , 'LHS_Values_BASIC.out']);
% Top 100 from lhs fits

%% =======================================================================
%      Experiment vals
% untreated
x0   = [   1:3  ];x1 = [4:7];
y0   = [   1.1,   0.89,  0.87];   
% TIP treated Expt animals
y1 = [ 1.43  , 1.25 , 1.68 , 1.32];
ExpPVL      = [ 6.62  ,  5.61  , 7.01  , 7.30 ];
ExpSSVL     = [ 4.62 , 4.48   , 4.18  , 5.52 ];
ExpCTRSSvl  = [  7.56 ,  8.2 , 7.85]

EXPT_Treated_SS = 1.0e+05 .*[     
                        0.0000    0.2907
                        0.0000    0.2584
                        0.0000    0.1061
                        0.0000    2.0507];
    
EXPT_Treated_MAX =         [   1     4210000
                               2      411000
                               3    10300000
                               4    20300000];

EXPT_CTR_SS       =  1.0e+08 .*[
                            0.0000    0.3510
                            0.0000    1.7902
                            0.0000    0.7126 ];

EXPT_CTR_MAX      =     [  1   144000000
                           2    19400000
                           3    21600000];
%% ========================================================================

cc = 1;
for k = [1,  4 , 3 ,  2 ] 
    
    fname1  = [ 'Model_PVL_A' , sprintf('%d' , k) , '_M0', '.out' ];
    fname1b = [ 'Model_SSL_A' , sprintf('%d' , k) , '_M0', '.out' ];
    fname1c = [ 'Analytical_SSL_A' , sprintf('%d' , k) , '_M0', '.out' ];
  
    fname2  = [ 'Model_PVL_A' , sprintf('%d' , k) , '_M1', '.out' ];
    fname2b = [ 'Model_SSL_A' , sprintf('%d' , k) , '_M1', '.out' ];
    fname2c = [ 'Analytical_SSL_A' , sprintf('%d' , k), '_M1' , '.out' ];
   

    data1    = importdata( [ ipath, fname1 ]);
    data1b   = importdata( [ ipath, fname1b ]);
    data1c   = importdata( [ ipath, fname1c ]);


    data2    = importdata( [ ipath, fname2 ]);
    data2b   = importdata( [ ipath, fname2b ]);
    data2c   = importdata( [ ipath, fname2c ]);
  
 
    d2   = ( data2 )./( data2b(:,1) );
    d1c   = ( data1 )./( (data1c(:,1)) );
 
    % Fits: M1 , SS 
    figure(12),...
        plot( x1(cc) ,  (d1c) , '.' ,'markersize', 18 , 'color' , [0.5, 0.5 0.5 ] ),hold on,...
        %'o' ,'markersize' , 8  , 'color' , [ 0.5 0.5 0.5 ] ), hold on,...
        h1= plot( x1(cc) ,  (d2) , '.' ,'markersize' , 18  , 'color' , [ 0.3922    0.5843    0.902]  ), hold on,...
        h2= plot( x1(cc) ,   log10(EXPT_Treated_MAX(k,2))/log10(EXPT_Treated_SS(k,2)) ,  '*' , 'markersize' , 35 , 'color'  , 'k' ); %'k*' ,'markersize' , 25 )
        h3=plot( x1(cc) ,  LHSdataTIP( cc ,3) , 's','markersize' , 18  , 'color' , [ 0.3922    0.5843    0.902]), hold on,...
        cc = cc + 1;

end




for i = 1:3
    fname0a = [ 'Model_PVL_A' , sprintf('%d' , i ) , '_CTR', '.out' ];
    fname0b = [ 'Model_SSL_A' , sprintf('%d' , i ) , '_CTR', '.out' ];
    fname0c = [ 'Analytical_SSL_A' , sprintf('%d' , i), '_CTR' , '.out' ];

    data1    = importdata( [ ipath, fname0a ]);
    data1b   = importdata( [ ipath, fname0b ]);
    data1c   = importdata( [ ipath, fname0c ]);

    d1   = data1./(data1b);
    figure(12),...
        plot( x0(i) , (d1) , '.' ,'markersize' , 18  , 'color' , [ 0.5 0.5 0.5 ]  ), hold on,...
        plot( x0(i) , log10( EXPT_CTR_MAX(i,2))/log10(EXPT_CTR_SS(i,2) ) , 'k*' ,'markersize' , 35  );
end

figure(12),...
    plot( [ x0,x1] , LHSdataBasic , 's','markersize' , 18  , 'color' , [ 0.4, 0.4 , 0.4 ] , 'MarkerEdgeColor', 'k', 'MarkerEdgeWidth', 2), hold on,...
    set( gca , 'fontsize' , 20 ),...
    xlim( [0.5 , 7.5 ] ),...
    
