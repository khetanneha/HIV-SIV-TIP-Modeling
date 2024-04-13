%% Neha Khetan
%  Draft Fig: peak:steady state and abs values from
%  Modified on 22 May 2023, mathematica -> Matlab
%  Modifying on 12 April 2024 : for rep
%  Cleaning up for submission: ----> S15B


clearvars; clc; close all;

ipath   = './'; 


% untreated
x0   = [   1:3  ];
y0   = [   1.1,   0.89,  0.87]   %[ 1.09 , 1.03 , 1.06]
% TIP treated Expt animals
x1   = [4:7];
y1   = [ 1.43  , 1.25 , 1.68 , 1.32];

%%=========================================================================
EXPT_Treated_TIP =  1.0e+05 .*[     0.0000    0.2907
                                    0.0000    2.0507
                                    0.0000    0.1061
                                    0.0000    0.2584];

EXPT_CTR_SS     =    1.0e+08 .*[    0.0000    0.3293
                                    0.0000    1.5722
                                    0.0000    0.6186 ];


cc = 1;
for k = [ 1, 4 , 3, 2 ] 
    %fname1  = [ 'Model_PVL_A' , sprintf('%d' , k) , '_M0', '.out' ];
    %fname1b = [ 'Model_SSL_A' , sprintf('%d' , k) , '_M0', '.out' ];
    %fname1c = [ 'Analytical_SSL_A' , sprintf('%d' , k) , '_M0', '.out' ];

    %fname2  = [ 'Model_PVL_A' , sprintf('%d' , k) , '_M1', '.out' ];
    %fname2b = [ 'Model_SSL_A' , sprintf('%d' , k) , '_M1', '.out' ];
    fname2c = [ 'Analytical_SSL_A' , sprintf('%d' , k), '_M1' , '.out' ];



    %data1    = importdata( [ ipath, fname1 ]);
    %data1b   = importdata( [ ipath, fname1b ]);
    %data1c   = importdata( [ ipath, fname1c ]);

    %data2    = importdata( [ ipath, fname2 ]);
    %data2b   = importdata( [ ipath, fname2b ]);
    data2c   = importdata( [ ipath, fname2c ]);


    % From S.S. | M1 - V and TIP
    figure(1),...
        plot( x1(cc) , data2c(:,2) , '.' ,'markersize', 18 , 'color' , [0.5, 0.5 0.5 ] ), hold on,...
        plot( x1(cc) , log10(EXPT_Treated_TIP(k,2)) , '*' , 'markersize' , 30 , 'color'  , 'b' ),hold on,...
    cc = cc + 1;

end


xpre = [ 1:3];
x0   = [ 1,2,3];

for i = 1:3
    %fname0a = [ 'Model_PVL_A' , sprintf('%d' , i ) , '_CTR', '.out' ];
    %fname0b = [ 'Model_SSL_A' , sprintf('%d' , i ) , '_CTR', '.out' ];
    fname0c = [ 'Analytical_SSL_A' , sprintf('%d' , i), '_CTR' , '.out' ];

    %data1    = importdata( [ ipath, fname0a ]);
    %data1b   = importdata( [ ipath, fname0b ]);
    data1c   = importdata( [ ipath, fname0c ]);

    figure(1),...
        plot( x0(i)    , data1c(:,1) , '.' ,'markersize', 18    , 'color' , [0.5, 0.5 0.5 ] ), hold on,...
         plot( x0(i) , log10(EXPT_CTR_SS(i,2)) , 'r*'  ,  'markersize' , 30  ),hold on,...
end
    set( gca , 'fontsize' , 24),...
        ylabel('Log_{10} Total SIV-Gag RNA copies/mL', 'fontsize', 20 )


