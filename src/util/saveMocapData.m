function saveMocapData(file, Y, initY, varY, frames) 

if (~exist('frames', 'var'))
    frames = 1:size(Y,1);
end

mask = [32:33 34 35:36 44:45 46 47:48 55 62];

% % % choose one of 2 rescaling schemes (must be consistent with load function
% % 
% if (varY(1) < varY(2))
%     Y(:,4:6) = Y(:,4:6)*sqrt(varY(2)/varY(1));
% end
% if (varY(1) < varY(3))
%     Y(:,1:3) = Y(:,1:3)*sqrt(varY(3)/varY(1));
% end
% % 
% Y(:,4:6) = Y(:,4:6)*sqrt(varY(2)/varY(1));
% Y(:,1:3) = Y(:,1:3)*sqrt(varY(3)/varY(1));

Y = Y(:,1:size(initY,2));

globY = Y(:,1:3);
globY(1,:) = initY(1:3) + globY(1,:);
Y(1,1:3) = initY(1:3);
Y(2,1:3) = globY(1,1:3);

for n=3:size(Y,1)
	Y(n,1:3) = sum(globY(1:n-1,1:3));
end

%Y(:,1:3) = zeros(size(Y,1),3);
tmpY = zeros(size(Y,1), 62); 
count = 1; 
for n=1:62
    if (~ismember(mask, n))
        tmpY(:,n) = Y(:,count); 
        count = count + 1; 
    end
end
Y = tmpY; 
matrix_to_amc(file, Y(frames,:));
