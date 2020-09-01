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

function inc_Reff = forwardmodel(s0,m,N)

    global ga rh de tspan model Rc xvec qvec pop k1 k2 social_prox si be drampup dmax drampdown
    
    % inidist          % number of days prior to social distancing
    % R0_init          % initial R0 (prior to social distancing)
    % p                % social proximity (= 1 - social distancing)
    % social_prox      % social proximity over time
    % Rc               % R0 over time (as a function of social distancing)
    % sus              % vector of susceptibility/connectivity factors
    % inc              % incidence over time
    % Reff             % effective R0 over time
    % incidence        % output of 'forwardmodel' which includes inc and Reff
    
    % NOTE: in model 3, sus changes as a function of social distancing,
    % whereas in the remaining models it is unchanged.
    
    p       = m(1);
    inidist = m(2);
    R0_init = m(3);
    
    [t,soln] = ode45('covidODE',tspan,s0);
    % ODE defining the SEIR model of disease transmission
    
    social_prox = ones(size(tspan));
    Rc = R0_init*ones(size(tspan));
    for i = inidist+1:inidist+drampup
        social_prox(i) = 1-(i-inidist)*(1-p)/drampup;
        Rc(i) = social_prox(i)*Rc(i);
    end
    for i = inidist+drampup+1:inidist+drampup+dmax
        social_prox(i) = p;
        Rc(i) = social_prox(i)*Rc(i);
    end
    for i = inidist+drampup+dmax+1:inidist+drampup+dmax+drampdown
        social_prox(i) = 1 - (inidist+drampup+dmax+drampdown-i)*(1-p)/drampdown;
        Rc(i) = social_prox(i)*Rc(i);
    end

    if(model == 1)
        % Variable susceptibility
        sus = xvec;
        inc = de*sum(soln(:,N+1:2*N),2)';
        D = sum((rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'),1)/(rh/de+1/ga);
        lam = be*social_prox/pop.*(ones(size(sus))*(rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'));
        Reff = lam./D.*(sus*(soln(:,1:N)'+si*soln(:,3*N+1:4*N)'));
        inc_Reff = cat(1,inc,Reff);
    elseif(model == 2)
        % Variable connectivity
        sus = xvec;
        inc = de*sum(soln(:,N+1:2*N),2)';
        D = sum((rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'),1)/(rh/de+1/ga);
        lam = be*social_prox/k1/pop.*(sus*(rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'));
        Reff = lam./D.*(sus*(soln(:,1:N)'+si*soln(:,3*N+1:4*N)'));
        inc_Reff = cat(1,inc,Reff);
    elseif(model == 3) 
        % Variable connectivity (reducing CV during social distancing)
        sus = 1 + social_prox'*(xvec-1);
        inc = de*sum(soln(:,N+1:2*N),2)';
        D = sum((rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'),1)/(rh/de+1/ga);
        lam = be*social_prox/k1/pop.*diag(sus*(rh*soln(:,N+1:2*N)'+soln(:,2*N+1:3*N)'))'./(sus*qvec)';
        Reff = lam./D.*diag(sus*(soln(:,1:N)'+si*soln(:,3*N+1:4*N)'))';
        inc_Reff = cat(1,inc,Reff);
        Rc = Rc.*((sus.^2*qvec)./(sus*qvec))'/(k2/k1);
    end
end

