function [ upstream, downstream ] = PostShock(config, downstream, shockSpeed)
% POSTSHOCK Calculate the post shock state of a shock moving with velocity
% shockSpeed into the downstream state described by downstream.
%
% Internal calculation assumes the shock moves to the left, so shockSpeed
% should be negative.
%
% The returned state is in the wave reference frame and must be converted
% back to the lab frame to actually put it to use, using
%  ps = SetVelocity(ps, ps(2)/ps(1) + shockSpeed)
%

% Transform the downstream state into the wave coordinate system
downstreamEnergy = (downstream(3) - downstream(2)^2/downstream(1)/2) / downstream(1);
downstream = SetVelocity(downstream, downstream(2)/downstream(1) - shockSpeed);

% Calculate the post-shock state using the conservation equations
% TODO 3 is kind of an arbitrary initial guess for ρ. Better do an initial guess!
rho = fsolve(@(rho) hugoniotEnergyFunctional(config, downstreamEnergy, downstream, rho), 3, optimset('Display', 'off', 'TolFun', eps, 'TolX', eps));
rayleighPressure = downstream(end-config.i.P) - downstream(2)^2 * (1/rho - 1/downstream(1));
energy = (rayleighPressure + downstream(end-config.i.P)) * (1/downstream(1) - 1/rho) / 2 + downstreamEnergy;

upstream = downstream;
upstream(1) = rho;
upstream(3) = energy * rho + upstream(2)^2 / rho / 2;
upstream(4:end-config.i.COUNT) = upstream(4:end-config.i.COUNT) * rho / downstream(1);
upstream = config.kineticsif.advance(upstream, 0);

function [ err ] = hugoniotEnergyFunctional(config, downstreamEnergy, downstream, rho)
	% RHS for the Post-Shock state calculation
	%
	% Tries to find the density for which the Rayleigh line and the Hugoniot
	% intersect in the P/V diagram. Calculates the pressure for both functions
	% and returns their difference.

	% Calculate pressure using the Rayleigh line
	rayleighPressure = downstream(end-config.i.P) - downstream(2)^2 * (1/rho - 1/downstream(1));

	% Calculate energy using the Hugoniot
	energy = (rayleighPressure + downstream(end-config.i.P)) * (1/downstream(1) - 1/rho) / 2 + downstreamEnergy;

	% Determine the pressure according to this energy
	test = downstream;
	test(1) = rho;
	test(2) = 0;
	test(3) = energy * rho;
	test(4:end-config.i.COUNT) = test(4:end-config.i.COUNT) * rho / downstream(1);
	test = config.kineticsif.advance(test, 0);

	% We aim to find the intersection of the Rayleigh line and the Hugoniot,
	% which both yielded a pressure, which has to be equal..
	err = rayleighPressure - test(end-config.i.P);

	% Penalty for values close to the initial value (A convergence error is
	% better than a useless result..)
	if abs(rho - downstream(1)) < 2e-1
		err = err + sign(err) * (2e-1 - abs(rho - downstream(1)))^2;
	end
