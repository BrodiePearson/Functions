function K = ripley_k_function(xy_pos,xK,box,method)
% KFUNCTION calculates Ripleys K function 
% K = kfunction(dataXY,xK,box,method) - returns vector K containing value
% of Ripley's K-function of dataXY in the distances in xK.
% dataXY - N-by-2 vector where N is number of datapoints. Each row
% corresponds to x and y coordinates of each datapoint
% xK - corresponds to the distances where K function should be computed.
% K is the same size as xK...
% box - rectangular boudnary of the data: box = [xlim1, xlim2, ylim1,
% ylim2]
% method - switch between edge correction. If method=0, no edge correction
% is applied. If method=1, datapoint is used for estimation of K(h) only if
% it is at least h units away from the box


if nargin<4 method=0; end
[N,k] = size(xy_pos);
if k~=2 error('dataXY must have two columns'); end

rbox = min([    xy_pos(:,1)'-box(1);
                box(2)-xy_pos(:,1)';
                xy_pos(:,2)'-box(3); 
                box(4)-xy_pos(:,2)']);
% rbox is the nearest distance of each datapoint to the box

DIST = squareform(pdist(xy_pos,'euclidean'));
DIST = sort(DIST);

if method==1 % edge correction...
K = zeros(length(xK),1);
Nk = length(K);
%wb = waitbar(0,'Computing Ripley''s K-function...');
for k=1:Nk
    %waitbar(k/Nk,wb);    
    I = find(rbox>xK(k));
    if ~isempty(I)
        K(k) = sum(sum(DIST(2:end,I)<xK(k)))/length(I);
    end
end
%close (wb);

elseif method==0 % no correction
    K = zeros(length(xK),1);
    for k=1:length(K)
        K(k) = sum(sum(DIST(2:end,:)<xK(k)))/N;
    end
end

% fprintf ('\b'); fprintf ('\b'); fprintf ('\b'); fprintf ('\b');
% fprintf ('100%%\n');

lambda = N/((box(2)-box(1))*(box(4)-box(3)));
K = K/lambda;