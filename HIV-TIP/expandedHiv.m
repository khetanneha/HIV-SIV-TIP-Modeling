
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neha Khetan, Basic TIP model- Weinberger et al 2003
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dY = expandedHiv( T , Y , PP )

% T  : Y1 -Healthy
% I  : Y2 -HIV infected
% V  : Y3 -HIV-virus
% It : Y4 - Tip infected
% Id : Y5 - HIV + TIP: dually infectd
% Vt : Y6 - TIP-infected



lam = PP.lam;
d   = PP.d;
k   = PP.k;
d2  = PP.d2;
c   = PP.c;
n   = PP.n;
P   = PP.P;
D   = PP.D;
d3  = D*d2;

dY = zeros( 6 , 1 );
dY(1) = lam - Y(1)*( d + ( k*Y(3) )  + ( k*Y(6) ) );
dY(2) = (k*Y(3)*Y(1)) - ( d2 * Y(2) );
dY(3) = ( n*d2 * Y(2) ) - ( c*Y(3) ) + ( D * n * d3 * Y(5) );
dY(4) = (k*Y(1)*Y(6)) - ( k*Y(3)*Y(4) ) - (d*Y(4));
dY(5) = (k*Y(3)*Y(4)) - ( d3 * Y(5) );
dY(6) = ( P^2 * D*n*d3*Y(5) ) - ( c*Y(6) );
end