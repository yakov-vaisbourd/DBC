%% Dual Block Coordinate (DBC) Algorithms for Solving the Total Variation 
% (TV) Denoising Problem.
% -----------------------------------------------------------------------
%% Copyright (2016): A.Beck, L.Tetruashvili, Yakov Vaisbourd and A.
%% Shemtov
% 
% dbc_tv is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%% Syntax
% function [x, H_all] = dbc_tv(d,varargin)
%% Description
% This function implements DAM method for solving the Total Variation (TV) 
% denoising problem:    
% 
% $$ \min_{\mathbf{x}\in \mathbf{R}^{m\times n}} H(\mathbf{x})\equiv 
% \frac{1}{2}\|\mathbf{x}-\mathbf{d}\|_F^2 + \theta\cdot TV(\mathbf{x}) $$
%                                                                               
% where $$ \mathbf{d} $$ is an 1D signal or a 2D image contaminated by noise.
% TV is the total-variation regularizer, which can either be the isotropic 
% (iso) or the anisotropic (l1) TV.
% For a further description of the model and the algorithm please refer to
% the DBC manual.

%% Input Arguments
%                                                                        
% * d -  The observed signal (1D or 2D) which is noisy                        
%                                                                        
% OPTIONAL INPUTS                                                        
%                                                                        
% * 'theta', scalar - The regularization parameter (default: 0.1).
%                                                                        
% * 'TV', 'iso' - Isotropic TV (default).
% * 'TV', 'l1' - Anisotropic TV .
%                                                                        
% * 'type', 'cyclic' - Cyclic version of the algorithm (default).
% * 'type', 'random' - Random version of the algorithm.
%                                                                                       
% * 'max_iter', scalar - Maximum number of iteration (default: 10,000).
%                                                                        
% * 'epsilon', scalar - Tolerance for relative error used in the stopping 
%                       criteria for the main problem (default: 0).
%                                                                        
% * 'epsilon_subproblem', scalar - Tolerance for relative error used in the
%                                  stopping criteria for the subproblems 
%                                  (default: 1e-4).
%
% * 'safeguard_subproblem', scalar - set to 1 in order to activate 
%                                    safeguard for the Newton method 
%                                    (default: 0).
%                                                                           
% * 'output', scalar - 1 if a report on the iterations is given, 0 if the  
% report is silenced (default).
%                                                                        
%                                                                        
%                                                                        
%% Output Arguments                                                              
%                                                                        
% * x       - The solution of the denoising problem.
%                                                                        
% * fun_all - Array containing the sequence the primal objective function 
%             values.
%                                                                        
%% Code Documentation
function [x,fun_all] = dbc_tv(d,varargin)

%% Assigning parameters according to arguments and/or default values

if nargin == 0
    error('Error: please specify a signal to denoise\n');
    
elseif nargin > 1+2*9
    error('Error: too many input arguments\n');
    
elseif rem(length(varargin),2)~=0
        error('Wrong list of arguments. Please refer to the function manual');    
end
    
%% Set default values and initialize parameters
tic;
[m, n] = size(d);
method = 'DAM';
theta = 0.1 ;
TV = 'iso';
type = 'cyclic';
maxiter = 200;
epsilon = 0;
sp_maxiter = 3;
sp_epsilon = 0;
sp_safeguard =0;
output = 0;
num = 0;



%% Arguments values
for i=1:2:length(varargin)-1
    switch lower(varargin{i})      
        case 'theta' ,              theta               = varargin{i+1};
        case 'tv' ,                 TV                  = varargin{i+1};
        case 'type' ,               type                = varargin{i+1};
        case 'maxiter' ,            maxiter             = varargin{i+1};
        case 'epsilon' ,            epsilon             = varargin{i+1};
        case 'epsilon_subproblem' , sp_epsilon          = varargin{i+1};
        case 'maxiter_subproblem' , sp_maxiter          = varargin{i+1};
        case 'safeguard_subproblem',sp_safeguard        = varargin{i+1};
        case 'output' ,             output              = varargin{i+1};
        case 'num' ,                num                 = varargin{i+1};        
        otherwise , ...
     error('Wrong list of arguments. Please refer to the function manual');
    end
end
   
%% Input Validation
if ~isscalar(maxiter) || ~isscalar(theta) || ~isscalar(epsilon)...
                                                       || ~isscalar(output)
    error('Please enter scalars as arguments for numerical values.');
end

if rem(maxiter,1)~=0 || maxiter<1 || rem(sp_maxiter,1)~=0 || sp_maxiter<1
    error...
 ('Please choose the maximum number of iterations as a positive integer.');
end

if epsilon<0 || sp_epsilon<0
    error('Please choose a positive value for epsilon.');
end

if output~=0 && output~=1
    error('Please choose 0 or 1 for "output".');
end
    
if epsilon==0
    epsilon=-1 ;
end 

if sp_epsilon==0
    sp_epsilon=[] ;
end 

%% Compute the decomposition transformations

[Trans, TransSizes,g_num] = decompose(m,n,TV,num);

% Problem variables and functions
sigma = 1;

%% Define functions 
% Objective function 
if n==1
    if strcmpi(TV, 'l1')
        ObjValue = @(x) 0.5*norm(x-d)^2+theta*norm(x(1:end-1)-x(2:end),1);
    elseif strcmpi(TV, 'iso')
        ObjValue = @(x) 0.5*norm(x-d)^2+theta*...
            sum(sqrt((x(1:end-2)-x(2:end-1)).^2+(x(2:end-1)-x(3:end)).^2));
    else
        error('Choose between isotropic (iso) and anisotropic (l1) TV.')
    end
else
    if strcmpi(TV, 'l1')
        ObjValue = @(x) 0.5*norm(x-d,'fro')^2+theta*...
                        (sum(sum(abs(x(1:end-1,:)-x(2:end,:))))+...
                         sum(sum(abs(x(:,1:end-1)-x(:,2:end)))));
    elseif strcmpi(TV, 'iso')
        ObjValue = @(x) 0.5*norm(x-d,'fro')^2+theta*...
                 (sum(sum(sqrt((x(1:end-1,1:end-1)-x(2:end,1:end-1)).^2+...
                            (x(1:end-1,1:end-1)-x(1:end-1,2:end)).^2)))+...
                                  sum(abs(x(1:end-1,end)-x(2:end,end)))+...
                                    sum(abs(x(end,1:end-1)-x(end,2:end))));
    else
        error('Choose between isotropic (iso) and anisotropic (l1) TV.')
    end
end 


% Opt. sub-problem solvers: in each case, we define the 2 (or more)  
%                           functions giving the solution for the 2 
%                           (or more) subproblems at each iteration

SUBPROB1 = {@(v,alpha) kron(min(alpha,abs(v(1:2:end-1,:)-...
             v(2:2:end,:))/2).*sign(v(1:2:end-1,:)-v(2:2:end,:)),[1 ; -1]);
         @(v,alpha) sp_newton(v,alpha,sp_maxiter,sp_epsilon,sp_safeguard)};

SUBPROB2 = @(x,y) -x+d-y;
        
% Set the indices of each sub-problem to be solved according to the
% problem dimenssions, the TV type and the decomposition.
alg_type = strcat(TV,'_',int2str(min(n,2)),'D');
switch alg_type
    case 'l1_1D'                
        SUBPROB1_INDX = ones(1,g_num);
    case 'l1_2D'
        SUBPROB1_INDX = [1,1,1,1];
    case 'iso_1D'
        SUBPROB1_INDX = 2*ones(1,g_num);
    case 'iso_2D'        
        SUBPROB1_INDX = [2,2,2;1,1,1];
end

    
%% Initialization
k=1 ;
x = d;
fun_all = zeros(maxiter+1,1);
H = ObjValue(x);
fun_all(k) = H;
diff = Inf ;
y = zeros(m,n,g_num);

order_last = ceil(g_num*rand);

if output
    fprintf('Solving with %s-%c     \n',method,type(1));
    fprintf('#iteration  function-value  relative-difference\n');
    fprintf('-----------------------------------------------\n');
    fprintf('%5d    %10.12f   \n', k-1, H);
end

while k<=maxiter && abs(diff)>epsilon
    if strcmp(type,'random')
        order = rem(order_last+cumsum(ceil((g_num-1)*rand(1,g_num)))-1,...
                                                                  g_num)+1;
        order_last = order(g_num);
    else
        order = 1:g_num;
    end  
    for i=order        
        temp_sum = sum(y(:,:,[1:i-1,i+1:end]),3);            
        x = d-temp_sum; % Auxiliary Expression
        x(Trans(1:TransSizes(i,1),i,1)) = ...
        x(Trans(1:TransSizes(i,1),i,1))-SUBPROB1{SUBPROB1_INDX(1,i)}...
                            (x(Trans(1:TransSizes(i,1),i,1)),theta);
        if size(SUBPROB1_INDX,1)>1
            x(Trans(1:TransSizes(i,2),i,2)) = ...
        x(Trans(1:TransSizes(i,2),i,2))-SUBPROB1{SUBPROB1_INDX(2,i)}...
                            (x(Trans(1:TransSizes(i,2),i,2)),theta);
        end
        y(:,:,i) = SUBPROB2(x,temp_sum);
    end
    k = k+1;
    Hprev = ObjValue(x) ;
    diff = H-Hprev;
    H = Hprev ;    
    fun_all(k) = H;

    if output
        fprintf('%5d    %10.12f        %10.12f', k-1, H, abs(diff));
        fprintf('\n') ;
    end
end
fun_all = fun_all(1:k);
end


