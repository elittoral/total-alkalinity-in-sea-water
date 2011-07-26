function [ E0 ] = EstimE0( W0, WA, WB, EA, EB, C, KNernst )

  A1 = ((W0 + WB) * exp(EB/KNernst) - (W0 + WA) * exp(EA/KNernst))/ (WB - WA);
     A0 = (W0 + WA) * exp(EA/KNernst) - A1 * WA;


     AT = (-A1/A0) * C / W0;

     HA = (WA * C - W0 * AT) / (W0 + WA) ;
     HB = (WB * C - W0 * AT) / (W0 + WB) ;
     E0 = (EA - KNernst * log(HA) + EB - KNernst * log(HB)) / 2;


end

