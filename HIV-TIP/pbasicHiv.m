% Neha Khetan, Basic Viral dynamics model- Novak 1999
% M. A. Nowak, R. M. May, Virus dynamics : mathematical principles of immunology and
% virology. (Oxford University Press, Oxford ; New York, 2000).

function dY = pbasicHiv( T , Y , PP )
    lam = PP.lam;
    d   = PP.d;
    k   = PP.k;
    d2  = PP.d2;
    c   = PP.c;
    n   = PP.n;
    dY = zeros( 3 , 1 );
    dY(1) = lam - Y(1)*( d ) - k*Y(3)*Y(1) ;
    dY(2) = (k*Y(3)*Y(1)) - ( d2 * Y(2) );
    dY(3) = ( n*d2 * Y(2) ) - ( c*Y(3) ) ;
end