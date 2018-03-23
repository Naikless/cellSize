function [ lambda ] = cellSizeGavrikov( P, T, X, mech )
% This function calculates the approximate detonation cell size lambda for
% given initial conditions. 
% It uses the NASA-CEA code to acquire the CJ - detonation properties and 
% the Rankine-Hugoniot relations to calculate the post shock (von Neumann) state.
% Induction time is calculated solving ODEs for a zero dimensional, constant volume reactor model.
% The relation between induction time and detonation cell size is established by a correlation model,
% described in the work of [Gavrikov, Efimenko and Dorofeev, 2000].

% INPUT:
% P = pressure in (Pa)
% T = temperature in (K)
% X = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'

% OUTPUT:
% lambda = detonation cell size in (m)
 
try 
    gas = Solution([mech '.xml']);
catch
    try
        gas = Solution([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end
%% Get required values from induction time calculation
[ tau, ~, shockSpeed ] = induction( P, T, X, mech );
[~, T_vN] = PostShockCantera(P, T, X, mech, shockSpeed);

%% Calculate effectivel activation energy for D1 = shockSpeed and D2 = 1.6* shockSpeed
shockSpeed_plus = 1.6 * shockSpeed;
shockSpeed_mean = (1.6+1)/2 * shockSpeed;

set(gas,'T',T,'P',P,'X',X);
[P_vN_plus , T_vN_plus ] = PostShockCantera(P, T, X, mech, shockSpeed_plus);
set(gas,'T',T,'P',P,'X',X);
[P_ps, T_ps] = PostShockCantera(P, T, X, mech, shockSpeed_mean);

tau_plus  = reactor_isoV(P_vN_plus,  T_vN_plus,  X, mech);
tau_mean = reactor_isoV(P_ps, T_ps, X, mech);

E_a = gasconstant * log(tau/tau_plus) / (T_vN^-1 - T_vN_plus^-1);

%% Calculate induction length
set(gas,'T',T,'P',P,'X',X);
rho0 = density(gas);
set(gas,'P',P_ps, 'T', T_ps, 'X', X);
rho1 = density(gas);
u_ps = shockSpeed_mean * rho0 / rho1;
delta = u_ps * tau_mean;


%% Approximate cell size through polynom [Dick, Ju, Lee, 2007]
lambda = GavrikovCorr( delta, E_a, T_ps, T_vN, T );
end

