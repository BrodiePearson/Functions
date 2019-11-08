function [ droxi,droeta ] = gradient_rho_terms(rho,pm,pn )
%computes gradient of any field located at rho points
droxi=zeros(size(pm));
droeta=zeros(size(pm));
droxi(:,2:end)=(rho(:,2:end)-rho(:,1:end-1)).*0.5.*(pm(:,2:end)+pm(:,1:end-1));
droeta(2:end,:)=(rho(2:end,:)-rho(1:end-1,:)).*0.5.*(pn(2:end,:)+pn(1:end-1,:));

%grad_rho_2=droxi.^2+droeta.^2;
end