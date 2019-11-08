% function out=make2dredspec(in,stdev)
clear all; clc;close all;
in = NaN([100,100]); stdev = 1;
% First, generate random white (gaussian) spectrum
out=fft2(randn(size(in)));

% construct frequency vectors and multiply
[nx,ny]=size(in);

dk=1/nx;
dl=1/ny;

[k,l]=meshgrid((-ny/2:ny/2-1)*dl,(-nx/2:nx/2-1)*dk);

filt=1./sqrt(k.^2+l.^2);
filt(~isfinite(filt))=0;

out=real(ifft2(out.*fftshift(filt)));

out= (out-mean(out(:)))./std(out(:));

% [psd kxx] = pwelch(out);
% psdxx = nanmean(psd,2);
% 
% loglog(kxx,psdxx)
% hold on
% loglog(kxx(10:end),kxx(10:end).^-2)
% 
% title('Spectra -2')
% xlabel('k')



