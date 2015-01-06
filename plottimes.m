clear all;

load m1;

X(:, 1) = timeTaken;
Y(:, 1) = loglik_all;

load m2;

X(:, 2) = timeTaken;
Y(:, 2) = loglik_all;

load m3;

X(:, 3) = timeTaken;
Y(:, 3) = loglik_all;


figure(1); clf;
plot(cumsum(X, 1), Y, '-');
legend('m=1','m=2','m=3');

print -dpsc2 -r300 timings.ps;
