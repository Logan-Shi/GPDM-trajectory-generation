function plotseries(X, segments, colour, missing)
%

hold on

q = size(X,2);
symb = 'o';
if nargin<=3
    missing = []; 
end

nmissing = setdiff(1:size(X,1), missing);
mask = segments-1;
mask(find(mask == 0)) = [];

if (q == 2)
    S = X;
    dS = S(2:end,:) - S(1:end-1,:);
    
    
    dS(mask,:) = zeros(length(mask),q);
    
    dS = [dS; [0 0]];
    quiver(S(:,1), S(:,2), dS(:,1), dS(:,2), 0, colour);
    plot(S(:,1), S(:,2), [colour symb]);
elseif (q >= 3)
    S = X;
    if size(X,1) > 1
        dS = S(2:end,:) - S(1:end-1,:);

        dS(mask,:) = zeros(length(mask),q);

        dS = [dS; zeros(1,q)];
       
        
        quiver3(S(nmissing,1), S(nmissing,2), S(nmissing,3), ...
            dS(nmissing,1), dS(nmissing,2), dS(nmissing,3), 0, colour);
        if (size(missing,2) > 0) 
            quiver3(S(missing,1), S(missing,2), S(missing,3), ...
                dS(missing,1), dS(missing,2), dS(missing,3), 0, 'g');
        end
    end
    plot3(S(nmissing,1), S(nmissing,2), S(nmissing,3), [colour symb]);
    plot3(S(missing,1), S(missing,2), S(missing,3), 'go');
end

if (q == 1)
    S = X;
    dS = S(2:end,:) - S(1:end-1,:);


    dS(mask,:) = zeros(length(mask),q);

    dS = [dS; 0];
    quiver(S(:,1), 0, dS(:,1), 0, 0, colour);
    plot(S(:,1), 0, [colour symb]);
end

hold off;
%set(gcf, 'Renderer', 'OpenGL');
side_max = max(max(X)); 
side_min = min(min(X));
axis([side_min side_max side_min side_max side_min side_max side_min side_max]); 
grid on;
axis equal;
