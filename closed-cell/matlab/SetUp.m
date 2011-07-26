function [ KNernst, E0, K2, BT, KB, K1P, K2P, K3P, KSI, ST, KS, Z, FT, KF, W0, W, KW, H ] = SetUp( S, T, V0, DAcid, NPts, V, E, C )
% INPUT
% S - salinity of titrated sample
% T - titration temperature (deg C)
% V0 - volume of sample titrated (cm3)
% DAcid - density of the acid titrant (g/cm3)
% NPts - number of titration points
% V - array of volumes of titrant used (cm3)
% E - Array of corresponding e.m.f.s (V)
% C - concentration of titrant acid (mol/kg)
%
% OUTPUT
% KNernst - Nernst parameter (appropriate to titration)
% E0 - Estimate of E0 of pH cell
% ----------
% W0 - mass of sample titrated (g)
% W - array containing titrant amounts (g)
% H - array containing [H'] = 10**((E0 - E)/K) based on the initial E0 estimate
% BT - total boron(mol/kg-soln)
% ST - total sulfate(mol/kg-soln)
% FT - total fluoride (mol/kg-soln)
% PT - total phosphate (mol/kg-soln)
% SiT - total silicate (mol/kg-soln)
%
% K1 - [H][HCO3]/[H2CO3]
% K2 - [H][CO3]/[HCO3]
% KB - [H][BO2]/[HBO2]
% K1P - [H][H2PO4]/[H3PO4]
% K2P - [H][HPO4]/[H2PO4]
% K3P - [H][PO4]/[HPO4]
% KSI - [H][SiO(OH)3]/[Si(OH)4]
% KS - [h][SO4]/[HSO4]
% KF - [H][F]/[HF]
% KW - [H][OH]
% Z - pH scale conversion factor [H] = [h](1 + ST/KS)

if (S >= 5)
    W0 = V0 * DensSW(S, T);
    [ BT, ST, FT ] = ConcnsSW(S);
    [ K1, K2, KB, K1P, K2P, K3P, KSI, KS, KF, KW ] = ConstsSW(S, T);
    Z = 1 + ST/KS;
else
    W0 = V0*DensNaCl(S, T);
    BT=0;ST=0;FT=0;
    [ K1, K2, KW ]= ConstsNaCl(S, T);
    Z = 1;
end

for i=1:NPts
    W(i) = V(i) * DAcid;
end

KNernst = 8.31451 * (273.15 + T) / 96485.309;
if (E(1) > E(NPts))
    KNernst = -KNernst;
end

E0 = EstimE0(W0, W(NPts-1), W(NPts), E(NPts-1), E(NPts), C, KNernst);

for i=1:NPts
    H(i)=exp( (E(i)-E0) / KNernst );
end
end

