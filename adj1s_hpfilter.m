function [ ycycle_adj, lm_1, kappa ] = adj1s_hpfilter(y,lm_2,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for "On adjusting the one-sided Hodrick-Prescott filter" (2020) 
% by Elias Wolf (FU-Berlin), Frieder Mokinski (Deutsche Bundesbank), and
% Yves Schüler (Deutsche Bundesbank)
%
% Version date 2022/01/14
% If you encounter any bug, please mail Yves Schüler at yves.schueler (at) bundesbank.de
%
% This program estimates the cyclical component of the adjusted one-sided HP
% filter (HP-1s*) by estimating the two-sided HP filter (HP-2s) on an expanding sample and keeping the last observation.
%
% For the exact references that are cited in this program, please see
% the paper.
%
%       Input:    
%                   y       - a Tx1 vector (y_1,...,y_T)', where T is the
%                   number of observations. 
%                       
%                   lm_2  - a scalar. This is the value of the smoothing
%                   parameter usually employed for the two-sided HP filter
%                   (HP-2s). This input is optional. The default is lm_2=1,600.
%
%                   opt    - 0/1. If opt=1, the adjustment parameters are
%                   obtained through the optimization routine (computationally more intensive). This input is optional. The default is
%                   opt=0 (computationally less intensive).
%
%      Output:  
%                   ycycle_adj  - a Tx1 vector with the adjusted extracted cyclical component
%                  
%                   lm_1          - a scalar. This is the adjusted value of the smoothing parameter used as an input to the one-sided HP filter
%
%                   kappa       - a scalar. This is the mulitplicative scaling factor, with wich the cyclical component of the
%                   one-sided HP filter is rescaled
%
%       Example:
%                   [ycycle_adj, lm_1, kappa]=adj1s_hpfilter(y,1600)
%                   yields a Tx1 vector of the extracted cyclical component
%                   using a smoothing parameter of 650 and scaling
%                   parameter of size 1.1513 for the one-sided HP filter
%                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1 || isempty(y), error('Time series is missing or empty'); end
if nargin < 2, lm_2= 1600; opt = 0; warning('Smoothing parameter not privided. Smoothing parameter is set to 1,600');  end
if nargin < 3,  opt = 0; end
if sum(isnan(y))>0, error('y contains missing values: please remove'); end
if size(y,2)>1, error('y contains more than one time series: the program is designed for one series only'); end
if size(y,1)<3, error('y contains less than 3 observations: please check your data input'); end

if opt == 1
    % finds adjustment parameters lm_1 and kappa by harmonizing PTF of HP-2s
    % and HP-1s
    T_ptf= 1000; %size of vector of filter weight for observation t=T
    
    grid = 2*pi*(0:1:floor(T_ptf/2))/T_ptf; %grid of spectrum
    
    ptf_hp2s = ((4*lm_2*(1-cos(grid)).^2./(1+4*lm_2*(1-cos(grid)).^2)).^2)'; % PTF of HP-2s in large samples (see King & Rebelo (1993))
    
    
    x = zeros(1,2); %starting value set to polynomials (see below)
    x(1) =  -5.071743 +  0.4080780*lm_2  -0.749537*lm_2^(1/2) + 14.50001*lm_2^(1/3) -50.36123*lm_2^(1/4) + 41.44369*lm_2^(1/5);
    x(2) =  0.984928 +  0.547209*(1/lm_2)  - 0.703636*(1/sqrt(lm_2)) + 3.291657*(1/(lm_2)^(1/3)) - 3.165415*(1/(lm_2)^(1/4)) + 1.759527 *(1/(lm_2)^(1/5));
    [x_hat,~] = fminsearch(@(x)fn_harmonize(x,T_ptf,ptf_hp2s,grid),x); %harmonization carried out here

    lm_1 = x_hat(1);
    kappa = x_hat(2);
else
    % polynomials providing lm_1 and kappa designed for values of lm_2 in
    % the range of 6.25 up to 1000000. Polynomials fitted on the basis of
    % values given in Table 5 in Appendix A
    lm_1 = -5.071743 +  0.4080780*lm_2  -0.749537*lm_2^(1/2) + 14.50001*lm_2^(1/3) -50.36123*lm_2^(1/4) + 41.44369*lm_2^(1/5);
    kappa = 0.984928 +  0.547209*(1/lm_2)  - 0.703636*(1/sqrt(lm_2)) + 3.291657*(1/(lm_2)^(1/3)) - 3.165415*(1/(lm_2)^(1/4)) + 1.759527 *(1/(lm_2)^(1/5));
end

T = length(y);
ycycle = zeros(T,1); %set first two observations to zero

for ii = 3:T

I = eye(ii); 
Q = diff(I,2);
A_inv = (I + lm_1*(Q')*Q)\I;

ycycle(ii,1) = y(ii,1) - A_inv(end,:)*y(1:ii,1);
end

ycycle_adj = kappa*ycycle;

 
end

%% function used in the harmonization rountine if opt=1
function [value] = fn_harmonize(x,T,ptf_hp2s,grid)

I = eye(T);
Q = diff(I,2);
A_inv = (I + x(1)*(Q')*Q)\I;

    k=[-(T-1):0]';
    for ii =1:length(grid)
           tf_hp1s(1,ii) = 1- ones(1,length(k))*(A_inv(end,:)'.*exp(1i*grid(ii)*k)); %transfer function of HP-1s given t=T observations
    end

    ptf_hp1s = x(2)^2*abs(tf_hp1s').^2;
    value =  100*sum((ptf_hp1s-ptf_hp2s).^2);
end


