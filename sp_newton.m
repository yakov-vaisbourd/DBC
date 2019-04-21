%% TV-isotropic
% Solve the subproblem emerged in the denoising version of the dual block 
% coordinate descent methods with the isotropic total variation regulizer.
% 
%% -----------------------------------------------------------------------
% Copyright (2016): A.Beck, L.Tetruashvili, Yakov Vaisbourd and A.
% Shemtov
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
%   function res = sp_newton(v,theta)
%   function res = sp_newton(v,theta,k_max)
%   function res = sp_newton(v,theta,[],epsilon)
%% Description
% This function solves simlulatniusly N 3-dimensional optimization 
% problems of the following form:   
%
% $$ \min_{\mathbf{x}\in \mathbf{R}^3} \frac{1}{2}\|\mathbf{x}
% -\mathbf{p}\|^2+\theta\sqrt{(x_1-x_2)^2+(x_2-x_3)^2}+
% \langle\mathbf{c},\mathbf{x}\rangle $$ 
%
% The output is equal to: $$ \mathbf{x}-\mathbf{v}. $$ 
% If some parameter values are not specified then the following default
% values are assumed:
%
% * epsilon = 0.
% * k_max = 3.
%
% If you wish to specify some parameter value without setting some of the
% parameters preceded him you should set [ ] for such parameters. For
% example, running the function without specifying the value of 
% _epsilon_ should be done as follows:
%
% x = sp_newton(v,theta,[ ],5)
%% Input Arguments
% * v - parameter vector (3N dimensional vector). This vector satisfies
%       v = p-c.
% * theta ($$ \theta $$) - parameter (scalar).
% * k_max - maximal number of iterations, a stopping criteria for the 
% Newton method.
% * epsilon ($$ \epsilon $$) - tollerance, a stopping criteria for the 
% Newton method.
% safeguard = activate safeguard for the Newton method.

%% Code Documentation
function res = sp_newton(v,theta,varargin)

% Assigning parameters according to arguments and/or default values
if nargin < 2
    error('Error: please specify the parameters v and theta.\n') ;    
elseif nargin > 2+2*3
    error('Error: too many input arguments\n') ;    
end

% Default values
k_max = 3;
epsilon = [];
safeguard = 0;

%%% Arguments values
for i=1:length(varargin)
    switch i        
        case 1 , k_max      = varargin{i};
        case 2 , epsilon    = varargin{i};
        case 3 , safeguard  = varargin{i};
    end
end

%%% Input validation
if size(v,2)~=1
    error('v should be a column vectors.');
end

if rem(size(v,1),3)~=0 || isempty(v)
    error('The dimension of the vector v should be divisible by 3');
end

if ~isempty(theta) && (~isscalar(theta) || theta<0)
    error('theta should be set to some positive scalar.');
end

if ~isempty(k_max) && (rem(k_max,1)~=0 || k_max<=0)
    error...
 ('Please choose the maximum number of iterations as a positive integer.');
end

if ~isempty(epsilon) && ~(epsilon>0)
    error('Please ch    oose epsilon to be positive number.')
end

% Set default values
if isempty(epsilon)
    epsilon = 0;
end
if isempty(k_max)
    if epsilon == 0
        k_max = 3;
    else
        k_max = Inf;
    end
end

%%% Initialization
n = size(v,1);
N = size(v,1)/3; % The quantity of the 3-dimensional optimization problems 
                 % to be solved simultaniuslly.
res = zeros(n,1);

%%% Main procedure

% Indicate the problems for which the unconstrained solution of the dual
% problem does not satisfy the constraints.
IND = 5*v(1:3:n-2).^2+2*v(2:3:n-1).^2+5*v(3:3:n).^2-2*v(1:3:n-2).*...
       v(2:3:n-1)-8*v(1:3:n-2).*v(3:3:n)-2*v(2:3:n-1).*v(3:3:n)>9*theta^2 ;
% Define auxiliary parameters.
t = zeros(N,2);
t(:,1) = v(3:3:n)-v(1:3:n-2);
t(:,2) = 2*v(2:3:n-1)-v(1:3:n-2)-v(3:3:n);


% Initialize the dual variables and additional parameters.
diff = zeros(N,1);
lambda = zeros(N,1);


lambda(IND) = max(sqrt((t(IND,1).^2+t(IND,2).^2)/2)/theta-1.5,0) ;                            
A = zeros(N,1) ;

% Safeguard - initialize lower and upper bounds for lambda
if safeguard
    phi_val = zeros(N,1);
    phi_val_l = zeros(N,1);
    lambda_l = zeros(N,1);%lambda;
    lambda_u = Inf*ones(N,1);
end

for k=1:k_max
    A(IND) = t(IND,1).^2.*(3+lambda(IND)).^2 + t(IND,2).^2.*...
                                                       (1+lambda(IND)).^2 ;   
    diff(IND) = A(IND).*(sqrt(A(IND)/2)/theta -(1+lambda(IND)).*...
                (3+lambda(IND)))./(t(IND,1).^2.*(3+lambda(IND)).^3+...
                                          t(IND,2).^2.*(1+lambda(IND)).^3);
    lambda(IND) = lambda(IND) + diff(IND) ;
        
    % Safeguard
    if safeguard        
        phi_val_l(IND) = 1/theta-1./sqrt(((t(IND,1)./(1+lambda_l(IND)))...
                                  .^2+(t(IND,2)./(3+lambda_l(IND))).^2)/2);        
        diff(phi_val_l>0) = 0;              
        phi_val(IND) = 1/theta-1./sqrt(((t(IND,1)./(1+lambda(IND))).^2+...
                                        (t(IND,2)./(3+lambda(IND))).^2)/2);
        temp = false(N,1);
        temp(IND) = phi_val(IND)>0 & lambda(IND)>lambda_l(IND);
        lambda_l(temp) = lambda(temp);
        temp = false(N,1);
        temp(IND) = phi_val(IND)<0 & lambda(IND)<lambda_u(IND);
        lambda_u(temp) = lambda(temp);
        temp = false(N,1);
        temp(IND) = lambda(IND)<lambda_l(IND);
        lambda(temp) = lambda_l(temp);
        temp = false(N,1);
        temp(IND) = lambda(IND)>lambda_u(IND);
        lambda(temp) = lambda_u(temp);    
    end
    
    IND = abs(diff)>epsilon ;
    if max(IND)==0
        break;
    end               
end

% Auxiliary variables and the "hard case".
p1 = t(:,1)./(2*(lambda+1));
p2 = (t(:,1)~=0).*t(:,2)./(2*(lambda+3))+(t(:,1)==0).*((abs(t(:,2))>3*...
               theta*sqrt(2))*theta.*sign(t(:,2))/sqrt(2)+(abs(t(:,2))...
                                            <=3*theta*sqrt(2)).*t(:,2)/6);

% A'*eta
res(1:3:n-2) = -p1-p2; 
res(2:3:n-1) = 2*p2;
res(3:3:n) = p1-p2;

end