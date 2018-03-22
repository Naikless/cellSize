function [ lambda ] = cellSize( P, T, X, mech )
% This function calculates the approximate detonation cell size lambda for
% given initial conditions. 
% It uses the NASA-CEA code to acquire the CJ - detonation properties and 
% the Rankine-Hugoniot relations to calculate the post shock (von Neumann) state.
% Induction time is calculated solving ODEs for a zero dimensional, constant volume reactor model.
% The relation between induction time and detonation cell size is established by a correlation model,
% described in the work of [Ng, Ju, Lee, 2007].

% INPUT:
% P = pressure in (Pa)
% T = temperature in (K)
% X = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'

% OUTPUT:
% lambda = detonation cell size in (m)
 
try 
    gas = importPhase([mech '.xml']);
catch
    try
        gas = importPhase([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end
%% Get required values from induction time calculation
[ tau, sigma_max, shockSpeed ] = induction( P, T, X, mech );
[P_vN, T_vN] = PostShockCantera(P, T, X, mech, shockSpeed);

%% Calculate non dimensional global activation energy epsilon [Schultz, Shepherd, 2000]
shockSpeed_plus  = 1.01 * shockSpeed;
shockSpeed_minus = 0.99 * shockSpeed;
set(gas,'T',T,'P',P,'X',X);
[P_vN_plus , T_vN_plus ] = PostShockCantera(P, T, X, mech, shockSpeed_plus);
set(gas,'T',T,'P',P,'X',X);
[P_vN_minus, T_vN_minus] = PostShockCantera(P, T, X, mech, shockSpeed_minus);

tau_plus  = reactor_isoV(P_vN_plus,  T_vN_plus,  X, mech);
tau_minus = reactor_isoV(P_vN_minus, T_vN_minus, X, mech);

epsilon = T_vN^-1 * log(tau_plus/tau_minus) / (T_vN_plus^-1 - T_vN_minus^-1);

%% Calculate induction length
set(gas,'T',T,'P',P,'X',X);
rho0 = density(gas);
set(gas,'P',P_vN, 'T', T_vN, 'X', X);
rho1 = density(gas);
u_vN = shockSpeed * rho0 / rho1;
delta = u_vN * tau;

%% Calculate CJ-plane particle velocity
[ ~, CJstate ]  = CJvelocity(P,T,X,mech);
u_CJ = CJstate.uCJ;

%% Approximate cell size through polynom [Dick, Ju, Lee, 2007]
lambda = cellSizeCorr(epsilon, delta, sigma_max, u_CJ);
end

