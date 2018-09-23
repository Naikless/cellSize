function [ tau, sigma_max ] = reactor_isoV( P, T, X, mech )
%% constructs ode system to solve zero dimensional constant volume reactor 
% combustion for given initial conditions 

% INPUT:
% P = pressure in (Pa)
% T =  temperature in (K)
% X = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'
 
% OUTPUTS:
% tau = induction time in (s)
% sigma_max = maximum normalized heat release in (1/s)

try 
    gas = Solution([mech '.xml']);
catch
    try
        gas = Solution([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end
set(gas,'P',P,'T',T,'X',X);

% calculate adiabatic flame temperature
equilibrate(gas,'UV');
T_eq = temperature(gas);
set(gas,'P',P,'T',T,'X',X);

% initial values
Y0=massFractions(gas);
T0=temperature(gas);
rho=density(gas);
sigma0 = 0;

% ode solver options
options=odeset('RelTol',1.e-6,'AbsTol',1e-12, 'Events', @(t,z) terminate(z,T_eq) );

% set initial and final time
t0=0; tf=10; % you may need to increase tf

%% solve ode system
[t,z]=ode15s( @(t,z) isoVsigma_ode(z, gas, rho) , [t0 tf], [T0; Y0; sigma0], options);


sigma_max = max(gradient(z(:,end),t));

%% determine ignition delay from maximum temperature gradient
[~,imaxgradT]=max(gradient(z(:,1),t)); % max grad T
tau = t(imaxgradT);
disp(['ignition delay: ' num2str(tau*1e6) ' microseconds'])

%% plot temperature in K
if length(t) < imaxgradT + 200
    endplot = length(t);
else
    endplot = imaxgradT + 200;
end
 figure(1), plot(t(1:endplot)*1e3,z(1:endplot,1),'k'), xlabel('time (ms)'), ylabel('T (K)')

end

function dzdt = isoVsigma_ode(z,g,rho)
% RHS for solving constant volume reactor
T=z(1); % temperature
y=z(2:end-1); % mass fractions

set(g,'T',T,'Density',rho,'Y',y); %update gas

Ms  = molecularWeights(g);
dy = ydot(g); % = netProdRates(g).*Ms/rho; change in mass fraction
ui  = (enthalpies_RT(g) - 1) *gasconstant * T ./ Ms; % specific internal energies
cv  = cv_mass(g); % specific heat

dT  = - sum(ui .* dy)/cv; % change in temperature

dsigma = sum ((meanMolecularWeight(g)./ Ms - enthalpies_RT(g) * gasconstant /  cp_mole(g)) .* dy);


dzdt=[dT;dy;dsigma]; % return

end

function [value, isterminal, direction] = terminate(z, T_eq)
value      = 0.95*T_eq - z(1);
isterminal = 1;   % Stop the integration
direction  = 0;
end