% function out=make2dspec(in,stdev,n)
clear all; close all; clc;
in = 100; stdev = 1; n = 100;
% First, generate random white spectrum
out=fft2(randn(size(in)));

% construct frequency vectors and multiply
[nx,ny]=size(in);

dk=1/nx;
dl=1/ny;

[k,l]=meshgrid((-ny/2:ny/2-1)*dl,(-nx/2:nx/2-1)*dk);

filt=(k.^2+l.^2).^(n/2);
filt(~isfinite(filt))=0;

out=real(ifft2(out.*fftshift(filt)));

out=stdev.*out./std(out(:));
