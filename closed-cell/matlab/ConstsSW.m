function [ K1, K2, KB, K1P, K2P, K3P, KSI, KS, KF, KW ] = ConstsSW( S, T )

TK = 273.15 + T;
    IS = 19.924*S/(1000 - 1.005*S);

    K1 = exp ( -2307.1266/TK + 2.83655 - 1.5529413*log(TK) ...
   +    (-4.0484/TK - 0.20760841)*sqrt(S) + 0.08468345*S ...
   -   0.00654208*S^1.5 + log(1 - 0.001005*S)    );

    K2 = exp (   -3351.6106/TK - 9.226508 - 0.2005743*log(TK) ...
   +    (-23.9722/TK - 0.106901773)*sqrt(S) + 0.1130822*S - 0.00846934*S^1.5 ... 
   + log(1 - 0.001005*S) );

    KB = exp (  (-8966.90 - 2890.53*sqrt(S) - 77.942*S + 1.728*S^1.5 -0.0996*S^2)/TK ... 
    + (148.0248 + 137.1942*sqrt(S) + 1.62142*S) ...
   +   (-24.4344 - 25.085*sqrt(S) -0.2474*S)*log(TK) ...
   +    0.053105*sqrt(S)*TK     );

    K1P = exp (   -4576.752/TK + 115.525 - 18.453*log(TK) ...
   +    (-106.736/TK + 0.69171)*sqrt(S) + (-0.65643/TK - 0.01844)*S     );


    K2P = exp ( -8814.715/TK + 172.0883 - 27.927*log(TK) ...
   +    (-160.340/TK + 1.3566)*sqrt(S) + (0.37335/TK - 0.05778)*S );


    K3P = exp (  -3070.75/TK - 18.141 ...
   +    (17.27039/TK + 2.81197)*sqrt(S) + (-44.99486/TK - 0.09984)*S );

    KSI = exp (  -8904.2/TK + 117.385 - 19.334*log(TK) ...
   +    (-458.79/TK + 3.5913)*sqrt(IS) ...
   +    (188.74/TK - 1.5998)*IS + (-12.1652/TK + 0.07871)*(IS)^2 ...
   +    log(1 - 0.001005*S));

    KW = exp (  -13847.26/TK + 148.9652 - 23.6521*log(TK) ...
   +    (118.67/TK - 5.977 + 1.0495*log(TK))*sqrt(S) - 0.01615*S );

    KS = exp ( -4276.1/TK + 141.328 - 23.093*log(TK) ...
   +    (-13856/TK + 324.57 - 47.986*log(TK)) * sqrt(IS) ...
   +    (35474/TK - 771.54 + 114.723*log(TK)) * IS ...
   -     2698*IS^1.5/TK + 1776*IS^2/TK + log(1 - 0.001005*S)  );

      KF = 1590.2/TK - 12.641 + 1.525 * sqrt(IS) + log(1 - 0.001005*S) ;
      KF = KF + log(1 + (0.1400/96.062)*(S/1.80655)/KS);
      KF = exp(KF);


end

