function [ P_vN, T_vN ] = PostShockCantera(P_0, T_0, X, mech, shockSpeed)
% this function calculates the fluid properties of the van Neumann - post 
% shock state using the Hugoniot and Rayleigh equations. 

% adapted for cantera based on a function written by Philipp Berndt 

% INPUT:
% P_0 = pressure in (Pa)
% T_0 =  temperature in (K)
% X = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'
% shockSpeed = CJ detonation velocity in (m/s)

% OUTPUT:
% P_vN = von Neumann pressure in (Pa)
% T_vN = von Neumann temperature in (K)
%
% ATTENTION: succesful solving for the Hugoniot energy functional strongly
% depends on an appropriate choice of the initial value for post-shock
% density RHO! RHO_init must be greater than the expected RHO, so 
% RHO_init = 15 has been chosen as default value. However, this value might
% have to be decreased for very weak shocks (due to dillution, very rich 
% mixtures, etc.)

try 
    gas = importPhase([mech '.xml']);
catch
    try
        gas = importPhase([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end
set(gas,'P',P_0,'T',T_0,'X',X);
RHO_0 = density(gas);
U_0 = intEnergy_mass(gas);

RHO_init = 15; % choose wisely, see above!
RHO = fzero(@(RHO) hugoniotEnergyFunctional(gas, RHO_0, P_0, U_0, shockSpeed, RHO), RHO_init, optimset('Display', 'off', 'TolFun', eps, 'TolX', eps));

rayleighPressure = P_0 - (shockSpeed*RHO_0)^2 * (1/RHO - 1/RHO_0);
energy = (rayleighPressure + P_0) * (1/RHO_0 - 1/RHO) / 2 + U_0;

setState_UV(gas,[energy 1/RHO]);
P_vN = pressure(gas);
T_vN = temperature(gas);


function [ err ] = hugoniotEnergyFunctional(gas, RHO_0, P_0, U_0, shockSpeed, RHO)
	% RHS for the Post-Shock state calculation
	%
	% Tries to find the density for which the Rayleigh line and the Hugoniot
	% intersect in the P/V diagram. Calculates the pressure for both functions
	% and returns their difference.

	% Calculate pressure using the Rayleigh line
	rayleighPressure = P_0 - (shockSpeed*RHO_0)^2 * (1/RHO - 1/RHO_0);

	% Calculate energy using the Hugoniot
	energy = (rayleighPressure + P_0) * (1/RHO_0 - 1/RHO) / 2 + U_0;

	% Determine the pressure according to this energy
	setState_UV(gas,[energy 1/RHO]);
    P_test = pressure(gas);


	% We aim to find the intersection of the Rayleigh line and the Hugoniot,
	% which both yielded a pressure, which has to be equal..
	err = rayleighPressure - P_test;

	% Penalty for values close to the initial value (A convergence error is
	% better than a useless result..)
	if abs(RHO - RHO_0) < 2e-1
		err = err + sign(err) * (2e-1 - abs(RHO - RHO_0))^2;
	end
