function b = BarycentricMap(u,v,w)
% calculates coordinates and plots barycentric map

% Inputs
% u,v,w = velocity compents in in column vector form nd x nt form



umean = nanmean(u(:));
vmean = nanmean(v(:));
wmean = nanmean(w(:));

up = u-umean;
vp = v-vmean;
wp = w-wmean;

% Find correlations
% for ii = 1:size(up,1)
%     upsq(ii,:) = up(ii,:).*up(ii,:);
%     vpsq(ii,:) = vp(ii,:).*vp(ii,:);
%     upvp(ii,:) = up(ii,:).*vp(ii,:);
%     upwp(ii,:) = up(ii,:).*wp(ii,:);
%     wpvp(ii,:) = wp(ii,:).*vp(ii,:);
%     wpsq(ii,:) = wp(ii,:).*wp(ii,:);
% end

    upsq = up.*up;
    vpsq = vp.*vp;
    upvp = up.*vp;
    upwp = up.*wp;
    wpvp = wp.*vp;
    wpsq = wp.*wp;

% %Set nans to zero
% upsq(isnan(upsq)) = 0;
% vpsq(isnan(vpsq)) = 0;
% upvp(isnan(upvp)) = 0;
% upwp(isnan(upwp)) = 0;
% wpvp(isnan(wpvp)) = 0;
% wpsq(isnan(wpsq)) = 0;

% Find means
mupsq = nanmean(upsq(:));
mvpsq = nanmean(vpsq(:));
mupvp = nanmean(upvp(:));
mupwp = nanmean(upwp(:));
mwpvp = nanmean(wpvp(:));
mwpsq = nanmean(wpsq(:));


% Computes TKE for each time step with dimensions 1 x time
% TKE = 1/2*(mupsq + mvpsq); 
TKE = 1/2*(mupsq + mvpsq + mwpsq); 

% Initialize anisotropy tensor
a = zeros(3,3,length(TKE));

% Compute components of anisotropy tensor
a(1,1,:) = mupsq./(2.*TKE)-1/3;
a(1,2,:) = mupvp./(2.*TKE);
a(1,3,:) = mupwp./(2.*TKE);
a(2,1,:) = mupvp./(2.*TKE);
a(2,2,:) = mvpsq./(2.*TKE)-1/3;
a(2,3,:) = mwpvp./(2.*TKE);
a(3,1,:) = mupwp./(2.*TKE);
a(3,2,:) = mwpvp./(2.*TKE);
a(3,3,:) = mwpsq./(2.*TKE)-1/3;


%GLAD_ANISOTROPY.a = a;   % Save matrix to structure


% Find eigenvalues and eigenvectors
for kk = 1:length(TKE)
    [V,D] = eig(a(:,:,kk));
    l = [D(1,1), D(2,2), D(3,3)];
    lambda1(kk) = max(l);
    lambda2(kk) = median(l);
    lambda3(kk) = min(l);

% Define coordinate system
c(kk,1) = lambda1(kk)-lambda2(kk);
c(kk,2) = 2*(lambda2(kk)-lambda3(kk));
c(kk,3) = 3*lambda3(kk) + 1;

end

% Plot Barycentric map
plotAnisotropicBarycentricMap(c,'bo');

end