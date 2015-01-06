doPrint = 1;

opt = out.opt;

% plot
myfigure(1);  clf;

subplot(2, 3, [1:3]);
hold on;
hist(out.X_all(:), 1:opt.nConfig);
title('Column counts');

if dat.N < 18

	% convert to row format
	nCellConfig = 2^dat.N; % total number of cell configurations
	X_count = zeros(1, nCellConfig);
	for it = opt.start_it : opt.nits
		X = out.X_all(:, it);
		Xmat = opt.Xtable(X, :);
		for k = 1 : params.K
			ind = bi2de(Xmat(:, k)') + 1;
			X_count(ind) = X_count(ind) + sum(out.w_all(k, :, it), 2);	
			X_row_all(k, it) = ind;
		end
	end


	subplot(2, 3, [4:6]);
	hold on;
	hnd = bar(0:nCellConfig-1, X_count(:)/(opt.nits-opt.start_it+1), 1);
	set(hnd, 'EdgeColor', 'r', 'FaceColor', 'r');
	ax = axis;
	xlim([-0.5 nCellConfig+0.5]);
	title('Row counts');

end

if doPrint
    myprint(outfile);
end

myfigure(2); clf;

subplot(2, 2, 1);
plot(out.loglik_all(100:10:end), 'b-');
%ylim([-400 0]);

subplot(2, 2, 2);
hist(out.loglik_all(100:10:end), 30);

%subplot(2, 2, 3);
%plot(out.res_all(100:end)', '.');
%ylim([0 1]);

%subplot(2, 2, 4);
%hist(out.res_all(100:end)', 30);
%xlim([0 1]);

if doPrint
	myprint(outfile, 1);
end

myfigure(3); clf;

colormap(hot);

for s = 1 : dat.nSamples
	subplot(1, dat.nSamples, s);
	imagesc(reshape(out.phi_all(s, :, opt.start_it:end), [dat.N opt.nits-opt.start_it+1]), [ 0 0.5 ]);
	set(gca, 'FontSize', 4);
end

if doPrint
    myprint(outfile, 1);
end

myfigure(6); clf;

for s = 1 : dat.nSamples

	for k = 1 : params.K

		subplot(params.K, dat.nSamples, dat.nSamples*(k-1) + s);
		hold on;
		plot( squeeze(out.w_all(k, s, opt.start_it:5:end)), 'b.', 'markersize', 1);
		set(gca, 'FontSize', 4);
		ylim([0 1]);
		title(['k: ' num2str(k) ', s: ' num2str(s)]);
	
	end
	
end

if doPrint
    myprint(outfile, 1); 
end


for s = 1 : dat.nSamples

    w_s_all = reshape(out.w_all(:, s, :), [params.K opt.nits]);
    w_s_all = sort(w_s_all(:, opt.start_it:2:end), 1, 'descend');

    myfigure(7+s); clf;

    subplot(2, 2, [1 2]);
    bar(w_s_all', 'stacked');
    title(['Sample: ' num2str(s)]);

    subplot(2, 2, 3);
    plot(w_s_all(1, :), w_s_all(2, :), 'k.');

    subplot(2, 2, 4);
    [h, xc, yc] = hist2d(w_s_all(1, :), w_s_all(2, :), 20);
    %contour(xc, yc, h, 30);
    surf(xc, yc, flipud(h));
    shading interp;
    view(2);

    if doPrint
        myprint(outfile, 1); 
    end

end
