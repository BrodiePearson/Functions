function HD = HD_Calc(lsf,tsf,r)
%% Info:

% USE:

%   To calculate the isotropic, homogenous, stationary rotational and 
% divergenct structure functions of found using VSF_Calc/Bin with
% par.direction = 'lt'

% INPUT: nd x nt = matrix of dimensions number drifters x number of time steps
%       l = nd x nt longitudinal structure function values (m/s)^2
%       t = nd x nt transverse structure function values (m/s)^2
%       r = nd x nt separation distances (km)


% OUTPUT: nd x nt = matrix with dimensions of number of velocity
% differences by number of time steps
%       l = nd x nt cublic spline interpolated longitudinal structure function values (m/s)^2
%       t = nd x nt cublic spline interpolated transverse structure function values (m/s)^2
%       rot = nd x nt cublic spline interpolated rotational structure function values (m/s)^2
%       div = nd x nt cublic spline interpolated divergence structure function values (m/s)^2
%       r = nd x nt cublic spline interpolated separation distances (km)

[lsfc,gof,output] = fit(r',lsf(:,1),'smoothingspline');
[tsfc,gof,output] = fit(r',tsf(:,1),'smoothingspline');


x = (min(r):max(r))';
yl = feval(lsfc,x);
yt = feval(tsfc,x);

yl(1) = lsf(1);
yt(1) = tsf(1);


for ii = 2:length(yt)
    
    int = trapz((1./x(1:ii)).*(yt(1:ii)-yl(1:ii)));   % Numerically integrate using trapezoid rule
    
    
    sfrot(ii-1) = yt(ii-1)+int;   % Rotational Component
    sfdiv(ii-1) = yl(ii-1)-int;   % Divergent Component
    
end

x = x(1:end-1);
yl = yl(1:end-1);
yt = yt(1:end-1);


% Save data to structure
HD.rot = sfrot;
HD.div = sfdiv;
HD.l = yl;
HD.t = yt;
HD.r = x;


end
