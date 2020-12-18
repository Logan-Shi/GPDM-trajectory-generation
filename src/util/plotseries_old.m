function plotseries(X, segments, colour, missing)

if (~exist('missing', 'var'))
	missing = [];
end

nmissing = setdiff(1:size(X,1), missing);

hold on;
q = size(X,2);

mask = segments-1;
mask(find(mask == 0)) = [];

if (q == 2)
    dX = X(2:end,:) - X(1:end-1,:);
    dX(mask,:) = zeros(length(mask),q);
    dX = [dX; [0 0]];
    quiver(X(nmissing,1), X(nmissing,2), dX(nmissing,1), ...
	dX(nmissing,2), 0, colour);
    plot(X(nmissing,1), X(nmissing,2), [colour 'o']);
    
	quiver(X(missing,1), X(missing,2), dX(missing,1), ...
	dX(missing,2), 0, 'g');
    plot(X(missing,1), X(missing,2), ['go']);
elseif (q >= 3)
    dX = X(2:end,:) - X(1:end-1,:);
    dX(mask,:) = zeros(length(mask),q);
    dX = [dX; zeros(1,q)];
	
    quiver3(X(nmissing,1), X(nmissing,2), X(nmissing,3), ...
	dX(nmissing,1), dX(nmissing,2), dX(nmissing,3), 0, colour);
    plot3(X(nmissing,1), X(nmissing,2), X(nmissing,3), [colour 'o']);
    
	quiver3(X(missing,1), X(missing,2), X(missing,3), ...
	dX(missing,1), dX(missing,2), dX(missing,3), 0, 'g');
    plot3(X(missing,1), X(missing,2), X(missing,3), ['go']);
end

hold off;
set(gcf, 'Renderer', 'OpenGL');
axis equal;
