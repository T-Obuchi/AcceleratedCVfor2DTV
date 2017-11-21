function [LOOE,ERR]=LOOEapprox_2DTV(x,y,A,Nx,Ny,lambda_T,delta,theta)
%--------------------------------------------------------------------------
% LOOEapprox_2DTV.m: An approximate leave-one-out error estimator for
% linear regression with l1 and total variation regularizations
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Compute and return an approximate leave-one-out error (LOOE) and its
%    standard error for linear regression penalized by l1 norm and
%    two-dimensional total variation (TV). 
%
% USAGE:
%    [LOOE,ERR] = LOOEapprox_2DTV(x,y,A,Nx,Ny,lambda_T,delta,theta)
%
% INPUT ARGUMENTS:
%    x           Explanatory variables (N=Nx*Ny dimensional vector). 
%                A two-dimensional (2D) image is expected in common cases.
%
%    y           Measurement result (M dimensional vector)
%
%    A           Measurement matrix (M*N dimensional matrix)
%
%    Nx          One side length of x in 2D.
%
%    Ny          Another side length of x in 2D.
%
%    lambda_T    Regularization weight of TV
%
%    delta       Softening constant of TV. Default value is 10^(-4).
%
%    theta       Threshold to determine clusters induced by TV. 
%                Default value is 10^(-12).
%
% OUTPUT ARGUMENTS:
%    LOOE        Approximate value of the leave-one-out error 
%
%    ERR         Approximate standard error of the leave-one-out error 
%
% DETAILS:
%    The following linear regression penalized by the l1 norm and total 
%    variation is considered:
%
%                \hat{x}=argmin_{x}( 1/2 ||y-A*x||_2^2 
%                       + lambda_l1*||x||_1 + lambda_T*T(x) ),
%
%    where T(x) is the two-dimensional total variation
%
%                T(x)=\sum_i\sqrt{ \sum_{j \partial i}(x_j-x_i)^2 },
%
%    and \partial i represents the set of neighbors indices of x_i. 
%    The input x is supposed to be the solution of this regression. 
%    The LOOE is defined by
%
%                LOOE=(1/(2M))*\sum_{mu=1}^{M}( y_{\mu}-\sum_{i=1}^{N}
%                A_{\mu i}\hat{x}_i^{\backslash \mu} )^2,
%
%    where \hat{x}^{\backslash \mu} is the solution of the above 
%    regression without the mu-th observation. 
%    This LOO solution \hat{x}^{\backslash \mu} is approximated 
%    from the full solution \hat{x}, yielding an approximate LOOE.
%    For the details of the method, see arXiv:1611.07197.
%
% REFERENCES:
%    Tomoyuki Obuchi, Shiro Ikeda, Kazunori Akiyama, Yoshiyuki Kabashima 
%    Accelerating cross-validation with total variation and its application
%    to super-resolution imaging
%    arXiv:1611.07197
%
% DEVELOPMENT:
%    18 Jul 2017: Original version was written.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default values of delta and theta
if nargin < 6
    error('more input arguments needed.');
end
if nargin < 7 || isempty(delta)
    delta = 10^(-4);
end
if nargin < 8 || isempty(theta)
    theta = 10^(-12);
end

% Geometric information in 2D
AROUNDs=AROUNDs_TV(Nx,Ny);

% Residual squares
RSS=(y-A*x).^2;                                         

% Total variation (TV)
TVterms=Allterms_TV_delta(x,Nx,Ny,delta);                          % All TV terms
[vec_ord,ind_ord]=sort(TVterms);                                   % Sort by the ascending order

% Compute TV and find its active components 
Nlckd=sum( TVterms < delta+theta );                                % Number of locked TV terms 
Active_TV=ones(Nx*Ny,1);                                           % Active terms in TV
Active_TV(ind_ord(1:Nlckd))=0;                                     % Smallest Nlckd terms are killed    

% Find clusters and Active and Inactive variables
Locked_cddt=[ind_ord(1:Nlckd),AROUNDs(ind_ord(1:Nlckd),:)];        
Locked_all=Find_Locked_all(Locked_cddt,Nx,Ny);                     % Set of clusters
[Active,Inactive]=Find_Active(x,Locked_all,Nx,Ny);                 % Active and inactive variables are enumerated

% Construct and merge Hessians
H_TV_full=Hessian_TV_full_delta(x,Nx,Ny,TVterms,Active_TV,delta);  % Full TV Hessian
H_TV_mgd=Merge_RC(x,H_TV_full,Active,Locked_all);                  % Merged TV Hessian using Active set info.
A_mgd=Merge_C(x,A,Active,Locked_all);                              % Merged measurement matrix
H_mgd=A_mgd'*A_mgd + lambda_T*H_TV_mgd;                            % Merged cost-function Hessian

% Approximate LOOE 
LOOEfactor=(1-diag(A_mgd(:,:)*( H_mgd\A_mgd(:,:)' ) )).^(-2);      
LOOE=mean(LOOEfactor.*RSS);
ERR=std(LOOEfactor.*RSS)/sqrt(length(y));
end