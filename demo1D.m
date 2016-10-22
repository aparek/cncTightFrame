%% 1D signal denoising using undecimated wavelet transform
% 
% Ankit Parekh (ankit.parekh@nyu.edu), NYU School of Engineering
%
% Reference:
% Convex denoising using non-convex tight frame regularization
% Ankit Parekh and Ivan W. Selesnick
% IEEE Signal Process. Lett., 2015
%
% Code available at: 
%             https://sites.google.com/a/nyu.edu/ankit-parekh/research

%% Initializations
clear; clc; close all
rmse = @(y,x) sqrt( sum( (y(:)-x(:)).^2) / numel(y) );
SNR = @(x,y) 20*log10( norm(x(:)) / norm(x(:)-y(:)) );
addpath Functions
addpath Functions/WaveletFunctions

%% Generate noisy signal
N = 512;
n = 0:N-1;
rng('default')
s = MakeSignal('Piece-Regular', N);
sigma = 4.0;
noise = sigma * randn(size(s));
y = s + noise;

figure(1), clf
subplot(2,1,1)
plot(s,'k')
title('Clean Signal')
box off
xlim([0 N])

subplot(2,1,2)
plot(y,'k')
title('Noisy signal')
box off
xlim([0 N])

%% Perform wavelet denoising

K = 3;                                                                      
J = 4;                                                                      
[AH,A,normA] = MakeTransforms('DWT',N,[J K]);                               
lam = cell(J,1);
lamL1 = cell(J,1);
a = cell(J,1);
beta = 2.3;
betaL1 = 1.8;
w = A(y);

for i = 1:J                                                                 
    lam{i} = beta * sigma ./ sqrt(2).^i * ones(size(w{i}));
    lamL1{i} = betaL1 * sigma ./ sqrt(2).^i * ones(size(w{i}));
    a{i} = 1./lam{i};
end

mu = 2;
Nit = 25;

[x,cost] = bp_ncvxUDWT(y,A,AH,J,lam,a,mu,Nit,'atan');
disp('Solved using non-convex regularization')
[xL1,costL1] = bp_ncvxUDWT(y,A,AH,J,lamL1,a,mu,Nit,'l1');
disp('Solved using convex regularization (L1 norm)')
%% Plot the cost function

figure(2), clf
plot(1:Nit, cost,'.-k')
ylabel('Objective Function')
xlabel('Number of iterations')
box off
xlim([0 Nit])
title('L1 regularization')

figure(3),clf
plot(1:Nit, costL1,'k');
ylabel('Objective Function')
xlabel('Number of iterations')
box off
xlim([0 Nit])
title('Nonconvex regularization')
%% Plot the denoised signals
figure(4), clf
gap = 100;
plot(n,y,'k',n,x-gap,'k',n,xL1-2*gap,'k');
text(256, 50, 'Noisy Signal','horizontalAlignment','center')
text(256, -50,...
    sprintf('Non-convex regularization (\\beta = 2.3, RMSE = %1.2f)',...
    rmse(s,x)),'horizontalAlignment','center')
text(256, -150, ...
    sprintf('L1 regularization (\\beta = 1.8, RMSE = %1.2f)',...
    rmse(s,xL1)), 'horizontalAlignment','center')
box off
xlim([0 N])
set(gca,'YTick',[],'YTickLabel',[])
xlabel('Time (n)')

