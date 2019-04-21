clear all;
clc;
clf;

% Generating a noisy test image
randn('seed',316);
Img = double(imread('cameraman.tif'))/255;
noise_level = 0.05;
d = Img + noise_level*randn(size(Img));

% Running with default values
[x,f_val] = dbc_tv(d);
subplot(1,3,1), imagesc(Img), colormap gray, title('Original image'), ...
    axis square, axis off;
subplot(1,3,2), imagesc(d), colormap gray, title('Noisy image'), ...
    axis square, axis off;
subplot(1,3,3), imagesc(x), colormap gray, title('Reconstructed image'),...
    axis square, axis off;


% Running 50 iterations of DAM-r for TV_l1 and theta=0.05
[x,f_val] = dbc_tv(d,...
                   'type',   'random',...                   
                   'theta',   0.05,...
                   'TV',     'l1',...
                   'maxiter', 50);
               
% Running 10 iterations with output
[x,f_val] = dbc_tv(d,'maxiter',10,'output', 1);               

% Restricting the number of iterations of the root finding Newton method to
% one.
t_inexact = tic;
[x,f_val] = dbc_tv(d,'maxiter',20,'maxiter_subproblem',1);               
tElap = toc(t_inexact);
disp([tElap,f_val(end)])

t_exact = tic;
[x,f_val] = dbc_tv(d,'maxiter',20);               
tElap = toc(t_exact);
disp([tElap,f_val(end)])
