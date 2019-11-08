function [ mpsd,f,kh ] = my_f_kh_spec_win_all(data,dt,dx,win)
%This function returns the 2D f,kh spectra
%It expects data to be size(Nt,Ny,Nx) and for Nx=Ny

if nargin<4
     disp('--hann windown is applied');
     %disp('--hamming windown is applied');
     win=1;
end
[nt,ny,nx]=size(data);
 data=data-mean(data(:));      % always remove the mean!
Kh_max = round(sqrt( (nx/2)^2 + (ny/2)^2 ));

kx = zeros(nx,1); 
  for i=1:nx 
    kx(i) = i-1;  
    if (kx(i)>nx/2)  
     kx(i) = nx+1-i; 
    end 
  end 
  
  ky = zeros(ny,1); 
  for j=1:ny 
    ky(j) = j-1;  
    if (ky(j)>ny/2)  
     ky(j) = ny+1-j; 
    end 
  end 
  if win==1
   disp('applying hann window');
  
   factor=sqrt(8/3);
    wt=factor*window(@hann,nt);
   wy=factor*window(@hann,ny); %ny=nx so wy=wx
  %more memory efficient
  [masky2,maskt2]=meshgrid(wy,wt);
  clear wy wt
  %data=maskx.*masky.*maskt.*data
  data=(permute(repmat(masky2,[1 1 nx]),[1 3 2])).*(repmat(masky2,[1 1 nx])).*(repmat(maskt2,[1 1 nx])).*data;
  clear masky2 maskt2
  end

  A=fftn(data);                % get the coefficients
  clear data;
  A=A/nt/ny/nx;       % normalize sanely
  mpsd3=A.*conj(A);             % the sum gives the variance right now
  clear A;
   
  if mod(nt,2)==0               % only half the data is good  
        stop=nt/2+1;                % we are making a onesided spectrum
    else
        stop=(nt+1)/2;              % stop in the right place in frequency
  end
  
  mpsd3=mpsd3(1:stop,:,:);
  if mod(nt,2)==0
     mpsd3(2:end-1,:,:)=2*mpsd3(2:end-1,:,:);  % multiply the right number by 2
  else
     mpsd3(2:end,:,:)=2*mpsd3(2:end,:,:);
  end
  
  mpsd=zeros(length(mpsd3(:,1,1)),Kh_max+1);
  
  
  Kh=zeros(ny,nx);
   for j=1:ny 
     for i=1:nx 
       Kh(j,i) = sqrt(kx(i)*kx(i) + ky(j)*ky(j)); 
     end 
   end 
   clear kx ky

  
  for j=1:ny
      for i=1:nx
          ik=round(Kh(j,i))+1;
          mpsd(:,ik)=mpsd(:,ik)+mpsd3(:,j,i);
      end
  end
  clear Kh;
  clear mpsd3  
  mpsd=mpsd*nt*dt*nx*dx;
  
  
   kh = [0:Kh_max]';
 kh=kh./(nx)./dx;
f=([0:nt-1])/nt/dt; 
  
  f=f(1:stop);
  
end

