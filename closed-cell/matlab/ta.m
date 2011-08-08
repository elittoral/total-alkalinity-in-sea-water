% Matlab version of fortran script writed by Andrew G. Dickson
% changelist:
% 1. removed MaxPts limitation
% 2. input/output in standard csv format (data separated by comma)
% 3. updated http://cdiac.ornl.gov/oceans/Handbook_2007.html

%% init
% csv: V/cm3, E/V, 
% S - salinity of sample, 
% PT[umol/kg] - total phosphate, 
% SiT[umol/kg] -total silicate, 
% T[oC] - temperature of sample when titrated, 
% W0[cm3] - volume of sample titrated, 
% C[mol/kg] - concentration of acid titrant, 
% DACID[g/cm3] - density of acid titrant 

clear all
IN=importdata('dane.in.csv');
%matrix cheatsheet (columns,rows)  ;)
meta=IN.data(1,3:9);
datos=IN.data(:,1:2);
S=meta(1);PT=meta(2);SiT=meta(3);T=meta(4);V0=meta(5);C=meta(6);DAcid=meta(7);
V=datos(:,1);E=datos(:,2);NPts=size(datos,1);
%% CALC
%MaxPts=100;
%NPar = 4; LWA = MaxPts*NPar + 5*NPar + MaxPts;

[ KNernst, E0, K2, BT, KB, K1P, K2P, K3P, KSI, ST, KS, Z, FT, KF, W0, W, KW, H ] = SetUp( S, T, V0, DAcid, NPts, V, E, C );

%% TODO this:
x0=[1 2 2 1];

%x(1)=F;
%x(2)=AT;
%x(3)=CT;
%x(4)=K1;

      %%%%%%% que pasa con W(I) y H(I) ????   NONLINEAR DATA-FITTING
% F = @(x,xdata)x(2) - x(3)*((x(4)*x(1)*H(I) + 2*x(4)*K2) / ((x(1)*H(I))^2 + x(4)*x(1)*H(I) + x(4)*K2)) ...
%        - BT/ (1 + x(1)*H(I)/KB) - PT*((K1P*K2P*x(1)*H(I) ...
%  + 2*K1P*K2P*K3P - (x(1)*H(I))^3) / ((x(1)*H(I))^3 + K1P*(x(1)*H(I))^2 ...
%  + K1P*K2P*x(1)*H(I) + K1P*K2P*K3P)) ...
%       - SiT/(1 + x(1)*H(I)/KSI) ...
%        + ST/(1 + KS*Z/(x(1)*H(I))) ...
%        + FT/(1 + KF/(x(1)*H(I))) ...
%        + (W0 + W(I))/W0 * (x(1)*H(I)/Z - KW/(x(1)*H(I)));

% plot(W,H,'ro')
% plot(V,E,'ro')

M = @(x,H)(W0*(x(1)*H).^2 - KW*W0*Z + W0*Z*x(1)*H.*  (   x(2) - x(3)*((x(4)*x(1).*H + 2*x(4)*K2) ./ ((x(1)*H).^2 + x(4)*x(1).*H + x(4)*K2)) ...
       - BT./ (1 + x(1)*H./KB) - PT*((K1P*K2P*x(1)*H ...
 + 2*K1P*K2P*K3P - (x(1)*H).^3) ./ ((x(1)*H).^3 + K1P*(x(1)*H).^2 ...
 + K1P*K2P*x(1)*H + K1P*K2P*K3P)) ...
      - SiT./(1 + x(1)*H/KSI) ...
       + ST./(1 + KS*Z./(x(1)*H)) ...
       + FT./(1 + KF./(x(1)*H))   ))./(KW*Z - (x(1)*H).^2);

[x,resnorm,~,exitflag,output] = lsqcurvefit(M,x0,H,W)

% hold on
% plot(M(x,H),H)
% hold off

%TOL = SQRT(DPMPAR(1)) TOL  - tolerance for fitting

%CALL LMDIF1(FCN, NPts, NPar, X, FVEC, TOL, INFO, IWA, WA, LWA) 


%% out
% LMDIF1 exit parameter 1
% Function calls 9
% Jacobian calls 28
% E0 = 0.393609 V
% AT = 2320.21 umol/kg
% CT = 2344.26 umol/kg
% pK1 = 5.9090
% s = 1.209 umol/kg
% V/cm3, E/V, -log[H], dH/[umol/kg]
% csv: -log[H], dH/(umol/kg), E0[V], AT[umol/kg], CT[umol/kg], pK1
OUT(:,1)=-log10(H*x(1));
OUT(1,2)=E0-KNernst*log(x(1));
OUT(1,3)=x(2);
OUT(1,4)=x(3);
OUT(1,5)=-log10(1-6*x(4));

csvwrite('dane.out.csv',OUT);
