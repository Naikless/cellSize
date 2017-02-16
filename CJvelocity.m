function [ V_CJ, CJstate ]  = CJvelocity(P1,T1,X1,mech)
% This function calculates the CJ-detonation velocity for given initial
% conditions. 

% It is based on (i.e. mostly copied from) the CJ-detonation part of the 
% NASA CEA Fortran95 code, available at
% http://www.grc.nasa.gov/WWW/CEAWeb/ceaHome.htm,
% which itself is based on "Calculation of Detonation Properties and Effect
% of Independent Parameters on Gaseous Detonations" by Frank J. Zeleznik 
% and Sanford Gordon, ARS Journal, April 1962.

% INPUT:
% P1 = pressure in (Pa)
% T1 = temperature in (K)
% X1 = mole fraction column vector
% mech = name of kinetic mechanism, e.g. 'gri30'

% OUTPUT:
% V_CJ = CJ detonation velocity in (m/s)
% CJstate = fluid properties at CJ-plane (i.e. where Ma=1)
%           consisting of: CJstate.P = pressure in (Pa)
%                          CJstate.T = temperature in (K)
%             CJstate.mole_fractions = mole fraction column vector
%                        CJstate.uCJ = particle velocity in (m/s)

try 
    gas = importPhase([mech '.xml']);
catch
    try
        gas = importPhase([mech '.cti']);
    catch
        error('Selected mechanism not found!');
    end
end

set(gas,'P',P1,'T',T1,'X',X1);

M1 = meanMolecularWeight(gas);
h1 = enthalpy_mass(gas);

P = 15 * P1;

% Acquire initial guess for T2 by finding flame temperature according to
% chapter 8.3 of CEA documentation
h_init = h1 + .75*gasconstant*T1/M1 * P/P1;
setState_HP(gas, [h_init P]);
equilibrate(gas,'HP');
T_init = temperature(gas);
T = T_init;

% Further improving initial guesses for T and P by recursion formulas
for i=1:3
    
    set(gas, 'T',T, 'P',P, 'X',X1);
    equilibrate(gas,'TP');
    M2 = meanMolecularWeight(gas);
    alpha = T1/T * M2/M1;
    [gamma_s , cp_eq] = partDeriv( T, P, X1, gas );  
    P = P1 * (1+gamma_s)/(2*gamma_s*alpha) * ...
        (1+(1-(4*gamma_s*alpha)/(1+gamma_s)^2)^0.5);
    r = alpha * P/P1;
    T = T1 * T_init/T1 - 0.75 * gasconstant/(M1*cp_eq) * 15 ...
        + (gasconstant*gamma_s)/(2*M1*cp_eq) * (r^2-1)/r * P/P1; 
end

% Main iteration loop
I=0;
while 1
    I = I+1;
    set(gas, 'T',T, 'P',P, 'X',X1);
    equilibrate(gas,'TP');
    M2 = meanMolecularWeight(gas);
    [ gamma_s , cp_eq , Dlvtp, Dlvpt ] = partDeriv( T, P, X1,gas ); 
    alpha = T1/T * M2/M1;
    r = alpha * P/P1;
    a11 = P1/P + gamma_s*r*Dlvpt;
    a12 = gamma_s*r*Dlvtp;
    a21 = 0.5*gamma_s*(r^2 - 1 - Dlvpt*(1+r^2)) + Dlvtp - 1;
    a22 = -0.5*gamma_s*Dlvtp*(r^2+1) - M2*cp_eq/gasconstant;
    b1 = P1/P - 1 + gamma_s*(r-1);
    b2 = M2*(enthalpy_mass(gas)-h1)/(gasconstant*T) - 0.5*gamma_s*(r^2-1);
    d = a11*a22 - a12*a21;
    x1 = (a22*b1-a12*b2)/d;
    x2 = (a11*b2-a21*b1)/d;
    
    delta = 1;
    temp = max([abs(x1) abs(x2)]);
    if temp > 0.4054652
        delta = 0.4054652/temp;
    end
    P = P*exp(x1*delta);
    T = T*exp(x2*delta);
    soundspeed_equil = sqrt(gasconstant*gamma_s*T/M2);
    V_CJ = r*soundspeed_equil;
    
    % Convergence test
    if I > 8
        display('Could not converge in 8 iterations')
        break
    elseif I <= 8 && temp < 0.5e-4
        break
    end
end
CJstate = struct('P',P,'T',T,'mole_fractions',moleFractions(gas), ...
    'uCJ', soundspeed_equil);
end

function [ gamma_s , cp_eq, Dlvtp, Dlvpt ] = partDeriv( T, P, X, gas )
% Calculates all needed variables related to partial derivatives, i.e 
% gamma_s (isentropic exponent), cp_eq (equilibrium specific heat at
% constant pressure), Dlvtp (d ln(v) / d ln(t) at constant pressure),
% Dlvpt (d ln(v) / d ln(p) at constant temperature). Derivatives are
% approximated via central differences.


deltaT = T*1e-6;
deltaP = P*1e-8;

v_P = zeros(1,2);
v_T = zeros(1,2);
h = zeros(1,2);

for i = 1:2
    
    sign = [-1 1];
    
    set(gas,'P',P+sign(i)*deltaP,'T',T,'X',X);
    equilibrate(gas,'TP');
    v_P(i) = density(gas)^-1;
    
    set(gas,'P',P,'T',T+sign(i)*deltaT,'X',X);
    equilibrate(gas,'TP');
    v_T(i) = density(gas)^-1;
    h(i) = enthalpy_mass(gas);
end

Dlvtp = log(v_T(2)/v_T(1))/log((T+deltaT)/(T-deltaT));
Dlvpt = log(v_P(2)/v_P(1))/log((P+deltaP)/(P-deltaP));
cp_eq = (h(2)-h(1))/(2*deltaT);

set(gas,'T',T, 'P',P, 'X',X);
equilibrate(gas,'TP');

gamma_s = -(Dlvpt + ...
    gasconstant * Dlvtp^2 / (cp_eq * meanMolecularWeight(gas)))^-1;
end