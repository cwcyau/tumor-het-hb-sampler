clear all;
close all;

load zare-0.mat;

outfile = 'plots/tumor-mixing.ps';
outfile_png = '~/Dropbox/external_papers/hammingball/diagrams/tumor-mixing.png';


K = 8;
s = 1;

for m = 1 : K
    timeTaken(m) = sum( out{m}.timeloop );
end

timeTaken(K+1) = sum(out_gibbs.timeloop);

%
% full posterior
%

figure(1); clf;

set(gcf, 'Renderer', 'Painters');
orient landscape;

mplot = 3;
nplot = K/2+2;
fontSz = 6;
markSz = 4;

subplot(mplot, nplot, 1);
hold on;
plot(1:K, timeTaken(1:K), 'ko-', 'markersize', markSz);
line([0 K], timeTaken(K+1)*[1 1], 'color', 'k', 'linestyle', '--');
%text(K, timeTaken(K+1), 'Gibbs');
set(gca, 'FontSize', fontSz, 'Box', 'On');
xlim([-0.5 K+0.5]);
xlabel('Hamming Ball Size');
ylabel('Time / s');


ind = opt.start_it:100:opt.nits;
n_ind = length(ind);

w_s_all = reshape(out_gibbs.w_all(:, s, :), [params.K opt.nits]);
w_s_all = sort(w_s_all(:, ind), 1, 'descend');

subplot(mplot, nplot, 2);
bar(w_s_all', 'stacked');
xlim([0 n_ind]);
ylim([0 1]);
set(gca, 'FontSize', fontSz);
ylabel('Component weights, w');
xlabel('Sample');
title(['Gibbs'], 'FontSize', fontSz+2);

j = 1;
for m = 2 : 2 : K

    w_s_all = reshape(out{m}.w_all(:, s, :), [params.K opt.nits]);
    w_s_all = sort(w_s_all(:, ind), 1, 'descend');

    subplot(mplot, nplot, j+2);
    bar(w_s_all', 'stacked');
    xlim([0 n_ind]);
    ylim([0 1]);
    title(['m = ' num2str(m)], 'FontSize', fontSz+2);
    set(gca, 'FontSize', fontSz);
    xlabel('Sample');
    colormap hsv;
    
%    subplot(mplot, nplot, j+1+(K/2+1));
%    plot(w_s_all(1, :), 'k.');
%    xlim([0 n_ind]);
%    ylim([0 1]);
%    set(gca, 'FontSize', fontSz);
%    title(['m = ' num2str(m)]);
%    xlabel('Sample');  
    
    j = j + 1;  

end

print(outfile, '-dpsc2', '-r600');



%
% max. component only
%
figure(2); clf;

set(gcf, 'Renderer', 'Painters');
orient landscape;

mplot = 2;
nplot = 4;
fontSz = 8;
markSz = 1;
lineSz = 0.5;

subplot(mplot, nplot, 1);
hold on;
hnd = bar(1:K, 100*(timeTaken(1:K)-timeTaken(K+1))/timeTaken(K+1), 0.25);
set(hnd, 'FaceColor', 'k');
%line([0 K+1], timeTaken(K+1)*[1 1], 'color', 'k', 'linestyle', '--');
%text(K, timeTaken(K+1), 'Gibbs');
set(gca, 'FontSize', fontSz, 'Box', 'On', 'XTick', 1:K);
xlim([0.5 K+0.5]);
xlabel('Hamming Ball Size, m');
ylabel('Relative Comp. Speed / %');
title('(a)', 'FontSize', fontSz+4);

ind = opt.start_it:100:opt.nits;
n_ind = length(ind);


cmap = hsv(K/2+1);

j = 1;
for m = 2 : 2 : K

    w_s_all = reshape(out{m}.w_all(:, s, :), [params.K opt.nits]);
    w_s_all = sort(w_s_all(:, ind), 1, 'descend');

    subplot(mplot, nplot, [2 3 4]);
    hold on;
    plot(0.8*w_s_all(1, :)+j+0.1, '.', 'color', 'k', 'markersize', markSz);
    if m < K
    	line([0 n_ind], [1 1]+j, 'color', 'k', 'LineWidth', lineSz);
    end
    xlim([0 n_ind]);
    ylim([0 5]);
    set(gca, 'FontSize', fontSz);
    xlabel('Sample');  
    ax = axis;
    text(-0.05*n_ind, j+0.5, ['m = ' num2str(m)], 'FontSize', fontSz, 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    
    j = j + 1;  

end

w_s_all = reshape(out_gibbs.w_all(:, s, :), [params.K opt.nits]);
w_s_all = sort(w_s_all(:, ind), 1, 'descend');

subplot(mplot, nplot, [2 3 4]);
hold on;
plot(0.8*w_s_all(1, :)+0.1, '.', 'color', 'k', 'markersize', markSz);
line([0 n_ind], [1 1], 'color', 'k', 'LineWidth', lineSz);
xlim([0 n_ind]);
ylim([0 5]);
set(gca, 'FontSize', fontSz, 'Box', 'On', 'YTick', [0.1 0.5 0.9 1.1 1.5 1.9 2.1 2.5 2.9 3.1 3.5 3.9 4.1 4.5 4.9], 'YTickLabel', {'0', '0.5', '1', '0', '0.5', '1', '0', '0.5', '1', '0', '0.5', '1', '0', '0.5', '1' }, 'TickLength', [0.005 0.005]);
xlabel('Sample');  
text(-0.05*n_ind, 0.5, 'Gibbs', 'FontSize', fontSz, 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
title('(b)', 'FontSize', fontSz+4);
    
print(outfile, '-dpsc2', '-r600');

system(['convert -density 600 -rotate 90 ' outfile ' ' outfile_png]);

gzip(outfile);
delete(outfile);

