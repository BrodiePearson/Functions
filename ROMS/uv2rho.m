function [varu_rho,varv_rho]=uv2rho(u,v)

% Converts u,v points to rho points from an Arakawa Grid in ROMS.

% Adapted from Pierrick Penven/ROMSTOOLS

[Mp,L,d,t]=size(u);
Lp=L+1;
Lm=L-1;
varu_rho=zeros(Mp,Lp,d,t);
varu_rho(:,2:L,:,:)=0.5*(u(:,1:Lm,:,:)+u(:,2:L,:,:));
varu_rho(:,1,:,:)=varu_rho(:,2,:,:);
varu_rho(:,Lp,:,:)=varu_rho(:,L,:,:);

[M,Lp,d,t]=size(v);
Mp=M+1;
Mm=M-1;
varv_rho=zeros(Mp,Lp,d,t);
varv_rho(2:M,:,:,:)=0.5*(v(1:Mm,:,:,:)+v(2:M,:,:,:));
varv_rho(1,:,:,:)=varv_rho(2,:,:,:);
varv_rho(Mp,:,:,:)=varv_rho(M,:,:,:); 

  

