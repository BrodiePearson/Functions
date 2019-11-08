function  [phi,psi] = flowfun(u,v,flag)
%  Uses command CUMSIMP (Simpson rule summation).

%  Kirill K. Pankratov, March 7, 1994.

% Check input arguments .............................................
issu=0; issv=0; isflag=0;    % For input checking
isphi = 1; ispsi = 1;        % Is phi and psi to be computed
if nargin==1, issu = isstr(u); end
if nargin==2, issv = isstr(v); end
if nargin==1&~issu, v=imag(u); end
if issv, flag = v; v = imag(u); isflag = 1; end 
if nargin==0|issu            % Not enough input arguments
  disp([10 '  Error: function must have input arguments:'...
  10 '  matrivces  U and V  (or complex matrix W = U+iV)' 10 ])
  return
end
if any(size(u)~=size(v))     % Disparate sizes
  disp([10 '  Error: matrices U and V must be of equal size' 10])
  return
end
if nargin==3, isflag=1; end
u = real(u);

 % Check the flag string . . . . . . . .
Dcn = str2mat('+','potential','phi');
Dcn = str2mat(Dcn,'-','streamfunction','psi');
if isflag
  lmin = min(size(flag,2),size(Dcn,2));
  flag = flag(1,1:lmin);
  A = flag(ones(size(Dcn,1),1),1:lmin)==Dcn(:,1:lmin);
  if lmin>1, coinc = sum(A'); else, coinc = A'; end
  fnd = find(coinc==lmin);
  if fnd~=[], if fnd<4, ispsi=0; else, isphi=0; end, end
end

phi = [];        % Create output
psi = [];

lx = size(u,2);  % Size of the velocity matrices
ly = size(u,1);

% Now the main computations .........................................
% Integrate velocity fields to get potential and streamfunction
% Use Simpson rule summation (function CUMSIMP)

 % Compute potential PHI (potential, non-rotating part)
if isphi
  cx = cumsimp(u(1,:));  % Compute x-integration constant
  cy = cumsimp(v(:,1));  % Compute y-integration constant
  phi = cumsimp(v)+cx(ones(ly,1),:);
  phi = (phi+cumsimp(u')'+cy(:,ones(1,lx)))/2;
end

 % Compute streamfunction PSI (solenoidal part)
if ispsi
  cx = cumsimp(v(1,:));  % Compute x-integration constant
  cy = cumsimp(u(:,1));  % Compute y-integration constant
  psi = -cumsimp(u)+cx(ones(ly,1),:);
  psi = (psi+cumsimp(v')'-cy(:,ones(1,lx)))/2;
end

 % Rename output if need only PSI
if ~isphi&ispsi&nargout==1, phi = psi; end