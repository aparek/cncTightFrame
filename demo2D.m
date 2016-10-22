%% 2D Image denoising using undecimated wavelet transform
% 
% Ankit Parekh (ankit.parekh@nyu.edu), NYU School of Engineering
% Reference:
% Convex denoising using non-convex tight frame regularization
% Ankit Parekh and Ivan W. Selesnick
% IEEE Signal Process. Lett., 2015

%% Initializations
clear; clc; close all
rmse = @(y,x) sqrt( sum( (y(:)-x(:)).^2) / numel(y) );
psnr = @(x,y) 10*log10 ( 1 / rmse(im2double(x),im2double(y))^2);
addpath Functions
addpath Functions/WaveletFunctions

%% Load Image and generate Noisy version of the image
S = rgb2gray(imread('peppers.png'));
S = double(S);

rng('default')
sigma = 50;
Y = S + sigma*randn(size(S));
J = 5;
[AH, A, normA] = MakeTransforms('2D-CDWT',size(Y),J);

figure(1), clf
subplot(1,2,1)
imshow(uint8(S));
title('Original Image')

subplot(1,2,2)
imshow(uint8(Y));
title(sprintf('Noisy Image, PSNR = %2.2f dB',psnr(uint8(S),uint8(Y))))

%% Denoise the image
Lam_l1 =  0.9*sigma;
Lam_ncvx = 1.1*sigma;
a = 1/(Lam_ncvx);   
mu = 2;
Nit = 15;
nit = 1:Nit;

x = bp_ncvx2DCWT(Y,A,AH,J,Lam_ncvx,a,mu,Nit,'atan');
xL1 = bp_ncvx2DCWT(Y,A,AH,J,Lam_l1,1,mu,Nit,'l1');

%% Plot the denoised images

figure(2), clf

subplot(2,2,1)
imshow(uint8(S))
title('Original Image')
box off


subplot(2,2,2)
imshow(uint8(Y))
title(sprintf('Noisy Image \nPSNR = %2.1f dB',psnr(uint8(S),uint8(Y))))
box off


subplot(2,2,3)
imshow(uint8(xL1))
title(sprintf('L1 \nPSNR = %2.1f dB',psnr(uint8(S),uint8(xL1))))
box off


subplot(2,2,4)
imshow(uint8(x))
title(sprintf('Non-convex \nPSNR = %2.1f dB',psnr(uint8(S),uint8(x))),...
    'horizontalAlignment','center')
box off
