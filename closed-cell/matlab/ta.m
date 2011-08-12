%% closed-cell total alkanity
% (2011) <http://www.elittoral.es elittoral S.L.N.E.> and <http://www.bioges.org Bioges>
%
% Matlab version of fortran script writed by Andrew G. Dickson (1994); 
% changelist:
%
% # removed MaxPts limitation
% # input/output in csv format 
% # updated to http://cdiac.ornl.gov/oceans/Handbook_2007.html
% # nonlinear data fitting using _levenberg-marquardt_ (default) or
% _trust-region-reflective_ method provided by lsqcurvefit Matlab function
%
%% init
% 
% * csv: V/cm3, E/V, 
% * S - salinity of sample, 
% * PT[umol/kg] - total phosphate, 
% * SiT[umol/kg] -total silicate, 
% * T[oC] - temperature of sample when titrated, 
% * W0[cm3] - volume of sample titrated, 
% * C[mol/kg] - concentration of acid titrant, 
% * DACID[g/cm3] - density of acid titrant 
clear all
IN=importdata('dane.in.csv');
%matrix cheatsheet (columns,rows)  ;)
meta=IN.data(1,3:9);
datos=IN.data(:,1:2);
S=meta(1);PT=1e-6*meta(2);SiT=1e-6*meta(3);T=meta(4);V0=meta(5);
C=meta(6);DAcid=meta(7);
V=datos(:,1);E=datos(:,2)/1000;NPts=size(datos,1);
%% CALC
%MaxPts=100;
%NPar = 4; LWA = MaxPts*NPar + 5*NPar + MaxPts;

[ KNernst, E0, K2, BT, KB, K1P, K2P, K3P, KSI, ST, KS, Z, FT, KF, ...
    W0, W, KW, H, AT ] = SetUp( S, T, V0, DAcid, NPts, V, E, C );
 
x0=[1 2e-3 2e-3 1e-6];

%x(1)=F;
%x(2)=AT;
%x(3)=CT;
%x(4)=K1;
% M = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);

M = @(x,H)(W0*(x(1)*H).^2 - KW*W0*Z + W0*Z*x(1)*H.*  ( x(2) - ...
x(3)*((x(4)*x(1).*H + 2*x(4)*K2) ./ ...
((x(1)*H).^2 + x(4)*x(1).*H + x(4)*K2)) ...
       - BT./ (1 + x(1)*H./KB) - PT*((K1P*K2P*x(1)*H ...
 + 2*K1P*K2P*K3P - (x(1)*H).^3) ./ ((x(1)*H).^3 + K1P*(x(1)*H).^2 ...
 + K1P*K2P*x(1)*H + K1P*K2P*K3P)) ...
      - SiT./(1 + x(1)*H/KSI) ...
       + ST./(1 + KS*Z./(x(1)*H)) ...
       + FT./(1 + KF./(x(1)*H))   ))./(KW*Z - (x(1)*H).^2);
   
%%
%
% $$ \textbf{W}=\frac{W_{0}(x_{1}\textbf{H})^{2} - K_{W}W_{0}Z+W_{0}Z\cdot
% x_{1}\textbf{H}\cdot \textit{fun}}{K_{W}Z-(x_{1}\textbf{H})^{2}} $$,
%
% $$ \textit{fun} = x_{2}-x_{3}\frac{x_{4}x_{1}\textbf{H}+2x_{4}K_{2}}{(x_{1}\textbf{H})^{2}+x_{4}x_{1}\textbf{H}+x_{4}K_{2}}-\frac{B_{T}}{1+\frac{x_{1}\textbf{H}}{K_{B}}}-P_{T}\frac{K_{1p}K_{2p}x_{1}\textbf{H}+2K_{1p}K_{2p}K_{3p}-(x_{1}\textbf{H})^{3}}{(x_{1}\textbf{H})^{3}+K_{1p}(x_{1}\textbf{H})^{2}+K_{1p}K_{2p}x_{1}\textbf{H}+K_{1p}K_{2p}K_{3p}}-\frac{S_{iT}}{1+\frac{x_{1}\textbf{H}}{K_{SI}}}+\frac{S_{T}}{1+\frac{K_{S}Z}{x_{1}\textbf{H}}}+\frac{F_{T}}{1+\frac{K_{F}}{x_{1}\textbf{H}}} $$
% (matlab latex bug)
% 

[x,resnorm,residual,exitflag,output,lambda,jacobian] = ...
    lsqcurvefit(M,x0,H,W,[],[],optimset('Algorithm','levenberg-marquardt'))


 plot(W,H,'ro')
 hold on
 plot(M(x,H),H)
 hold off

%TOL = SQRT(DPMPAR(1)) TOL  - tolerance for fitting
%CALL LMDIF1(FCN, NPts, NPar, X, FVEC, TOL, INFO, IWA, WA, LWA) 

%% out
% 
% * LMDIF1 exit parameter 1
% * Function calls 9
% * Jacobian calls 28
% * E0 = 0.393609 V
% * AT = 2320.21 umol/kg
% * CT = 2344.26 umol/kg
% * pK1 = 5.9090
% * s = 1.209 umol/kg
% * V/cm3, E/V, -log[H], dH/[umol/kg]
% * csv: -log[H], dH/(umol/kg), E0[V], AT[umol/kg], CT[umol/kg], pK1
F=x(1);AT=x(2);CT=x(3);K1=x(4);
OUT(:,1)=-log10(H*x(1));
OUT(1,2)=E0-KNernst*log(x(1));
OUT(1,3)=1e+3*x(2);
OUT(1,4)=1e+3*x(3);
OUT(1,5)=-log10(1e-6*x(4));

csvwrite('dane.out.csv',OUT);

%% TODO list
% 
% # real data test
% # CALC revision
% # fix syntax bugs
