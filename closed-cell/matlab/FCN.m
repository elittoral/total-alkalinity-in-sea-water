function [ FVEC ] = FCN( NPts, Npar, X )  %M, N

F=X(1);
AT=X(2);
CT=X(3);
K1=X(4);

for i=1:NPts  %100
    FVEC(i)=AT - CT*((K1*F*H(I) + 2*K1*K2) / ((F*H(I))^2 + K1*F*H(I) + K1*K2)) ...
       - BT/ (1 + F*H(I)/KB) - PT*((K1P*K2P*F*H(I) ...
 + 2*K1P*K2P*K3P - (F*H(I))^3) / ((F*H(I))^3 + K1P*(F*H(I))^2 ...
 + K1P*K2P*F*H(I) + K1P*K2P*K3P)) ...
      - SiT/(1 + F*H(I)/KSi) ...
       + ST/(1 + KS*Z/(F*H(I))) ...
       + FT/(1 + KF/(F*H(I))) ...
       + (W0 + W(I))/W0 * (F*H(I)/Z - KW/(F*H(I)));

    
end

end