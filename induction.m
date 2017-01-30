function [ tau, sigma_max, shockSpeed ] = induction( P, T, X, mech )
% This function calculates the induction time of a CJ detonation for given 
% initial conditions. 

% INPUT:
% P = pressure in (Pa)
% T =  temperature in (K)
% X = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'

% OUTPUT:
% tau = induction time in (s)
% sigma_max = maximum normalized heat release in (1/s)
% shockSpeed = CJ detonation velocity in (m/s)

try 
    gas = importPhase([mech '.xml']);
catch
    try
        gas = importPhase([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end
%% get shockSpeed via CJ calculation
set(gas,'T',T,'P',P,'X',X);
X = moleFractions(gas);
shockSpeed = CJvelocity(P,T,X,mech);
disp(['Shock speed is ' num2str(shockSpeed) ' m/s']);

%% get von Neumann state parameters
set(gas,'T',T,'P',P,'X',X);
[P_vN, T_vN] = PostShockCantera(P, T, X, mech, shockSpeed);

%% solve ode for reaction
[tau , sigma_max] = reactor_isoV(P_vN, T_vN, X, mech);