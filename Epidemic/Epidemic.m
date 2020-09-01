% @ copyright
% Authors:
%   Ricardo Aguas
%   Rodrigo M Corder
%   Jessica G King
%   Guilherme Goncalves
%   Marcelo U Ferreira
%   M Gabriela M Gomes
%
% This work is protected under the @Attribution-NonCommercial 4.0 International intellectual property license.
% You are free to:
%   Share - copy and redistribute the material in any medium or format
%   Adapt - remix, transform, and build upon the material Under the following terms:
%   Attribution - You must give appropriate credit to the authors, and indicate if any changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.
%   NonCommercial - You may not use the material for commercial purposes.
%   ShareAlike - If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

clear
close all
clc

global pop ga de rh N xvec qvec ndays drampup inidist tspan var qvec p ...
    k1 k2 model drampdown dmax Rc R0_init si ph

addpath('ODE')
addpath('Distributions')

%% Choose model

model = input('\n*** Choose a model: ***\n\n 1. Variable susceptibility \n 2. Variable connectivity \n 3. Variable connectivity (reducing CV during social distancing) \n\n Model: ');
switch model
	case 1
        model_name = 'Variable susceptibility';
    case 2
        model_name = 'Variable connectivity';
    case 3
        model_name = 'Variable connectivity (reducing CV during social distancing)';
    otherwise
        disp('Error: Invalid number')
        return
end

%% Define epidemiological parameters

fprintf('\n*** Insert parameters: ***\n\n')

R0_init = input(' Initial R0 = ');
if(R0_init<0) 
    fprintf('\n Error: invalid number \n\n')
    return 
end
CV = input(' Coefficient of variation in susceptibility or exposure (> or = 0): CV = ');
var = CV^2;
if(var<0)
    fprintf('\n Error: invalid number \n\n')
    return 
end
d = input(' Social distancing (0 - 1) = ');
if(d<0 || d>1) 
    fprintf('\n Error: invalid number \n\n')
    return 
end
inidist = round(input(' Time (in days) to initial social distancing measures = '));
if(inidist<0) 
    fprintf('\n Error: invalid number \n\n')
    return 
end

% default parameters
ga    = 1/4;            % rate of recovery or death
de    = 1/4;            % rate of progression from exposed to infectious
rh    = 0.5;            % transmissibility of exposed individuals (0 to 1)
rep   = 0.1;            % rate of reporting
ph 	  = 0;              % death fraction
si 	  = 0;              % reinfection factor
ndays = 366;            % total number of days
drampup = 21;           % length (in days) of the ramp-up of distancing measures
drampdown = 120;        % length (in days) of the ramp-down of distancing measures
dmax = 30;              % length (in days) of maximum distancing
dt = 1;                 % unit of time (days)
tspan = (0:dt:ndays-1); % time span
pop = 60e6;             % population size
p = 1-d;                % social proximity (= 1 - social distancing)

%% Define distribution of susceptibility for heterogeneous models

Dist_gamm;  % Susceptibility distribution

k1 = xvec*qvec;     % first central moment
k2 = xvec.^2*qvec;  % second central moment

%% Initial conditions

inf = pop/5e6;
% proportion of infections at the start of the epidemic 

s0= zeros(4*N,1);
s0(1:N)       = (pop)*qvec - (inf/rep)/(1-exp(-de))*qvec - (inf/rep)*qvec;
s0(N+1:2*N)   = (inf/rep)/(1-exp(-de))*qvec;
s0(2*N+1:3*N) = (inf/rep)*qvec;
s0(3*N+1:4*N) = (0/rep)*pop*qvec;
% Initial number of inidividuals within each class (S, I or R), 
% by susceptibility factor

%% Forward model

% Incidence function
incidence=@(m)forwardmodel(s0,m,N);

% Model input
fit = [p, inidist, R0_init];
res_interv = incidence(fit);

% Model output
inc = res_interv(1,:)*rep;
Reff = res_interv(2,:);


%% Plot

figure(1); 
clf('reset')

subplot(2,1,1)
plot(tspan,inc/pop*100,'-','color','blue','linewidth',1.5)
ylim([0 1.3*max(inc/pop*100)])
xTicks = [0:60:length(tspan)];
set(gca,'xtick',xTicks,'xticklabel',{},'fontsize',12)
ylabel('Daily new cases (%)','fontsize',12)
title(model_name,'FontWeight','Bold','fontsize',14)

subplot(2,1,2)
plot([0 ndays],[R0_init R0_init],'--','color', 'blue','linewidth',1.5)
hold on
plot([0 ndays],[1 1],'--','color', 'k','linewidth',1.5)
hold on
plot(tspan,Rc,'-','color','blue','linewidth',1.5)
hold on
plot(tspan,Reff,'-','color','m','linewidth',1.5)
ylim([0 R0_init*2])
xTicks = [0:30:length(tspan)];
set(gca,'xtick',xTicks,'fontsize',12)
xtickangle(45)
legend('R0','R0 = 1','Rc','Reff')
