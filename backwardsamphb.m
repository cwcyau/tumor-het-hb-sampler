function x = backwardsamphb(Pi, phi, u, loc, M, T)

% do backward sampling
v = zeros(1, T);
x = zeros(1, T);
v(T) = randsample(1:M, 1, 1, phi(:, T) );
x(T) = loc(u(T), v(T));

for t = (T-1):-1:1
	p = zeros(M, 1);
	xi = loc(u(t+1), v(t+1));
	for j = 1 : M
		xj = loc(u(t), j);
		p(j) = phi(j, t)*Pi(xj, xi);
	end
	v(t) = randsample(1:M, 1, 1, p);
	x(t) = loc(u(t), v(t));
end
