function [phi, phi_sum, delta] = forwardhb(Pi, obslike_X, u, loc, M, T)

% do forward filtering
phi = zeros(M, T);
phi_sum = zeros(1, T);
delta = zeros(M, T);

for i = 1 : M
	xi = loc(u(1), i);
	phi(i, 1) = obslike_X(xi);
end
phi_sum(1) = sum(phi(:, 1));
phi(:, 1) = phi(:, 1)/phi_sum(1);

tmp = zeros(M, 1);
for t = 2 : T
	for i = 1 : M
		xi = loc(u(t), i);
		for j = 1 : M
			xj = loc(u(t-1), j);
			tmp(j) = phi(j, t-1)*Pi(xj, xi)*obslike_X(t, xi);
		end	
		[ phi(i, t), delta(i, t) ] = max(tmp, [], 1);
	end
	phi_sum(t) = sum(phi(:, t));
	phi(:, t) = phi(:, t)/phi_sum(t);
end
