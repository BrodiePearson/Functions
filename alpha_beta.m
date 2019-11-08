function [alpha,beta]=alpha_beta(T,S,rho0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%computes alpha and beta for the linear approximation of the EOS
% t and s are the local salinity values
%  [N,M,L]=size(T);

sqrtS=sqrt(S);
Q01=6.793952e-2; Q02=-9.095290e-3;
Q03=+1.001685e-4; Q04=-1.120083e-6; Q05=+6.536332e-9;
U00=+0.824493; U01=-4.08990e-3 ; U02=+7.64380e-5;
U03=-8.24670e-7; U04=+5.38750e-9; V00=-5.72466e-3;
V01=+1.02270e-4; V02=-1.65460e-6; W00=+4.8314e-4;

cff=1./rho0;
 
	alpha=-cff*( Q01+T.*( 2.*Q02+T.*( 3.*Q03+T.*(4.*Q04 +T.*5.*Q05 )))+S.*( U01+T.*( 2.*U02+T.*(3.*U03 +T.*4.*U04 ))+sqrtS.*( V01+T*2.*V02)));
        beta= cff*( U00+T.*(U01+T.*(U02+T.*(U03+T.*U04)))+1.5*(V00+T.*(V01+T.*V02)).*sqrtS+2.*W00*S );
       


