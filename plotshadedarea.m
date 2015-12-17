function varargout = plotshadedarea(x, y, options)

% just plot one line
if size(y, 1) == 1
    plot(x, y, options);
end

% plot shaded area
if size(y, 1) == 2 
    px=[x, fliplr(x)]; % make closed patch
    py=[y(1,:), fliplr(y(2,:))];
    patch(px, py, 1, 'FaceColor', options, 'EdgeColor', 'none');
end
 
if size(y, 1) == 3 % also draw mean
    px=[x, fliplr(x)];
    py=[y(1,:), fliplr(y(3,:))];
    patch(px, py, 1, 'FaceColor', options, 'EdgeColor', 'none');
    plot(x, y(2,:), options);
end
 
alpha(0.2); % make patch transparent