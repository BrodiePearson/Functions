function [mpsd,ky,kx]=my2dpsd_r(data,dy,dx,win)
  % this function returns a the four quadrant 2D psd and the corresponding
  % wavenumber vectors kx and ky, for data of size (Ny,Nx)
  
  if nargin<4
     disp('--hann windown is applied');
     win=1;
  end
  
  
  sd=size(data);
  %ny=sd(1); nx=sd(2);
  ky=([0:sd(1)-1])/sd(1)/dy;        % get the frequencies in cp
  kx=([0:sd(2)-1])/sd(2)/dx;        % get wave number is cp
  data=data-mean(mean(data));      % always remove the mean!
  
%   win=hamming(sd(1));
%   data=repmat(win,[1 sd(2)]).*data;
  if win==1
    factor=sqrt(8/3);
    
    wc=factor*window(@hann,sd(1));
    wr=factor*window(@hann,sd(2));
    [maskr,maskc]=meshgrid(wr,wc);
    clear wc wr
   %  w=maskr.*maskc;
    data=maskr.*maskc.*data;
    clear maskr maskc
  end
  
  A=fft2(data);                    % get the coefficients
  A=A/sd(1)/sd(2);                 % normalize sanely
% get the spectral coefficients
  mpsd=A.*conj(A);                 % the sum gives the variance right now
  if mod(sd(1),2)==0               
     ky=[-sd(1)/2+1:sd(1)/2]/sd(1)/dy;  % in wave number
    iflip1=sd(1)/2+2;
    nl1=sd(1)/2-1;
    nr1=nl1+2;             
  else
    ky=[-(sd(1)-1)/2:(sd(1)-1)/2]/sd(1)/dy; 
    iflip1=(sd(1)+3)/2;
    nl1=(sd(1)-1)/2;
    nr1=nl1+1;
  end
  if mod(sd(2),2)==0                  % stop and flip in the right place 
    kx=[-sd(2)/2+1:sd(2)/2]/sd(2)/dx;  % in wave number
    iflip2=sd(2)/2+2;
    nl2=sd(2)/2-1;
    nr2=nl2+2;
  else
    kx=[-(sd(2)-1)/2:(sd(2)-1)/2]/sd(2)/dx; 
    iflip2=(sd(2)+3)/2;
    nl2=(sd(2)-1)/2;
    nr2=nl2+1;
  end
  ls(:,1:nl2)=mpsd(:,iflip2:end);    % and flip them around in wave number
  rs(:,1:nr2)=mpsd(:,1:iflip2-1);
  clear mpsd
  mpsd=[ls rs];
  us(1:nl1,:)=mpsd(iflip1:end,:);
  bs(1:nr1,:)=mpsd(1:iflip1-1,:);
  clear mpsd
  mpsd=[us; bs];

  mpsd=mpsd*dy*dx*sd(1)*sd(2);         % normalize for spectral density

