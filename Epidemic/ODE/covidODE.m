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

function dsdt = f(t,s)
 
    global model ga de rh pop N xvec qvec drampup p inidist R0_init k1 k2 dmax drampdown K si ph be

    % inidist          % length of the period prior to social distancing
    % drampup          % length of the rampup period of social distancing
    % dmax             % length of the maximum social distancing period
    % drampdown        % length of the rampdown period of social distancing
    % prox             % social proximity (= 1 - social distancing)
    % social_prox      % social proximity over time
    
    % ga               % rate of recovery/death
    % de               % rate of progression from exposed to infectious
    % rh               % transmissibility of exposed individuals (0 to 1)
    % si               % reinfection factor
    % sus              % vector of susceptibility/connectivity factors
    % be               % transmission coefficient
    % lam              % force of infection
    
    dsdt = zeros(4*N,1);
    
    % Social proximity over time
    if(t<=inidist)
        prox = 1; 
    elseif(t>inidist && t<=inidist+drampup)
        prox = 1-(t-inidist)*(1-p)/drampup;
    elseif(t>inidist+drampup && t<=inidist+drampup+dmax)
        prox = p;
    elseif (t>inidist+drampup+dmax && t<=inidist+drampup+dmax+drampdown)
        prox = 1 - (inidist+drampup+dmax+drampdown-t)*(1-p)/drampdown;
    else
        prox = 1;
    end
    
    
    if(model == 1)
        % Variable susceptibility
        sus = xvec;
        be = R0_init/(rh/de+1/ga)/k1;
        lam = be*prox/pop*ones(size(sus))*(rh*s(N+1:2*N)+s(2*N+1:3*N));
    elseif(model == 2)
        % Variable connectivity
        sus = xvec;
        be = R0_init/(rh/de+1/ga)*k1/k2;
        lam = be*prox/pop*sus*(rh*s(N+1:2*N)+s(2*N+1:3*N))/k1;
    elseif(model == 3) 
        % Variable connectivity (reducing CV during social distancing)
        sus = 1 + prox*(xvec-1);
        be = R0_init/(rh/de+1/ga)*k1/k2;
        lam = be*prox/pop*sus*(rh*s(N+1:2*N)+s(2*N+1:3*N))/(sus*qvec);
    end
   
    
    dsdt(1:N)       = - lam.*sus'.*s(1:N);                                      % Susceptible
    dsdt(N+1:2*N)   =   lam.*sus'.*(s(1:N)+si*s(3*N+1:4*N)) - de*s(N+1:2*N);    % Exposed
    dsdt(2*N+1:3*N) =   de*s(N+1:2*N)    - ga*s(2*N+1:3*N);                     % Infectious
    dsdt(3*N+1:4*N) =   (1-ph)*ga*s(2*N+1:3*N) - si*lam.*sus'.*s(3*N+1:4*N);           % Recovered

end