%% ========================================================================
% Neha Khetan, June 2023
% HIV-TIP - Bioreactor modeling for design criterion for amplification
% of spontaneously emerged DIPs in an evolutionary experimental scheme
% Design criterion, conditions; THIS is used for TIP-HIV paper
%% ========================================================================

function dY = TipmodelBioreactor( T , Y , PP )
    % 
    lam = PP.lam;
    d   = PP.d;
    k   = PP.k;
    d2  = PP.d2;
    c   = PP.c;
    n   = PP.n;
    P   = PP.P;
    D   = PP.D;
    d3  = D*d2;
    h   = PP.h;
    dr  = PP.dr;
    Pr  = PP.pr;
    drV = PP.drV;
    h0  = PP.h0;
    
    dY = zeros( 7 , 1 );
    dY(1) = lam   + Y(1)*h0*( ( h -( Y(1)+Y(2) + Y(4) + Y(5) ) )/h )  - (d*Y(1)) - k*Y(1)*( Y(3) + Y(6)  )- (dr*Y(1));
    dY(2) = (k*Y(3)*Y(1)) - ( d2 * Y(2) ) - (dr*Y(2));
    dY(3) = ( n*d2 * Y(2) ) - ( c*Y(3) ) + ( D * n * d3 * Y(5) ) - (drV*Y(3));
    dY(4) = (k*Y(1)*Y(6) )  + Y(4)*h0*(  (h - (Y(1)+Y(2) + Y(4) + Y(5) ))/h ) - ( k*Y(3)*Y(4) ) + (Pr*k*Y(1)*Y(3))- (d*Y(4))- (dr*Y(4));
    dY(5) = (k*Y(3)*Y(4)) - ( d3 * Y(5) ) - (dr*Y(5));
    dY(6) = ( P^2 * D*n*d3*Y(5) ) - ( c*Y(6) ) -  (drV*Y(6)) ;    
    dY(7) = (dr*Y(1)) + (dr*Y(2)) + (dr*Y(3)) + (dr*Y(4)) + (dr*Y(5)) + (dr*Y(6));
    %dY(7) =  (drV*Y(3))+ (drV*Y(6));
end


