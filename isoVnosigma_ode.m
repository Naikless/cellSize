function dzdt = isoV_ode(t,z,g,rho)
% ode rhs for constant pressure reactor
T=z(1); % temperature
y=z(2:end); % mass fractions

set(g,'T',T,'Density',rho,'Y',y); %update gas

Ms  = molecularWeights(g);
dy = ydot(g); % = netProdRates(g).*Ms/rho; change in mass fraction
ui  = (enthalpies_RT(g) - 1) *gasconstant * T ./ Ms; % specific internal energies
cv  = cv_mass(g); % specific heats

dT  = - sum(ui .* dy)/cv; % change in temperature

dzdt=[dT;dy]; % return
