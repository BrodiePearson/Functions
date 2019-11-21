function [u_rho,v_rho, w_rho]=uvw2rho(u,v,w)

% Converts u,v points to rho points from an Arakawa Grid in ROMS.

% Adapted from Pierrick Penven/ROMSTOOLS

[Mp,L,N,t]=size(u);
Lp=L+1;
Lm=L-1;
u_rho=zeros(Mp,Lp,N,t);
u_rho(:,2:L,:,:)=0.5*(u(:,1:Lm,:,:)+u(:,2:L,:,:));
u_rho(:,1,:,:)=u_rho(:,2,:,:);
u_rho(:,Lp,:,:)=u_rho(:,L,:,:);

[M,Lp,N,t]=size(v);
Mp=M+1;
Mm=M-1;
v_rho=zeros(Mp,Lp,N,t);
v_rho(2:M,:,:,:)=0.5*(v(1:Mm,:,:,:)+v(2:M,:,:,:));
v_rho(1,:,:,:)=v_rho(2,:,:,:);
v_rho(Mp,:,:,:)=v_rho(M,:,:,:); 

[M,L,N,t]=size(w);
Np=N+1;
Nm=N-1;
w_rho=zeros(M,L,Nm,t);
w_rho(:,:,1:Nm,:)=0.5*(w(:,:,1:Nm,:)+w(:,:,2:N,:));
% w_rho(:,:,1,:)=w_rho(:,:,2,:);
% w_rho(:,:,N,:)=w_rho(:,:,N,:); 

  

