function [ E0, AT ] = EstimE0( W0, WA, WB, EA, EB, C, KNernst )
% Calculate Gran function (W0+W)exp(E/K) and fit to y=a0+a1*x
% A1 = (yB - yA) / (xB - xA)
% A0 = yA - A1*xA
  A1 = ((W0 + WB) * exp(EB/KNernst) - (W0 + WA) * exp(EA/KNernst))/ (WB - WA);
     A0 = (W0 + WA) * exp(EA/KNernst) - A1 * WA;

% Calculate estimate of AT
     AT = (-A1/A0) * C / W0;

% Calculate [H] at those 2 points and hence an average E0
     HA = (WA * C - W0 * AT) / (W0 + WA) ;
     HB = (WB * C - W0 * AT) / (W0 + WB) ;
     E0 = (EA - KNernst * log(HA) + EB - KNernst * log(HB)) / 2;


end

