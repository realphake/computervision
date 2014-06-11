function [M,S, pointViewMatrix] = ChainFME( )
    [~,matches,f1,f2] = FME(sprintf('frame00000%03d.png',1), sprintf('frame00000%03d.png',2), 'ransac');
%     f1(1,:) = f1(1,:) - mean(f1(1,:));
%     f1(2,:) = f1(2,:) - mean(f1(2,:));
%     f2(1,:) = f2(1,:) - mean(f2(1,:));
%     f2(2,:) = f2(2,:) - mean(f2(2,:));
    
    matchIndexes = matches;
    pointViewMatrix(1,:) = f1(1,matchIndexes(1,:));
    pointViewMatrix(2,:) = f1(2,matchIndexes(1,:));
    pointViewMatrix(3,:) = f2(1,matchIndexes(2,:));
    pointViewMatrix(4,:) = f2(2,matchIndexes(2,:));
    
    for i = 2:30
        [~,matches,~,f2] = FME(sprintf('frame00000%03d.png',i), sprintf('frame00000%03d.png',i+1), 'ransac');
%         f2(1,:) = f2(1,:) - mean(f2(1,:));
%         f2(2,:) = f2(2,:) - mean(f2(2,:));
        newLine = zeros(1,size(matchIndexes,2));
        newPointViewLine = zeros(2,size(matchIndexes,2));
        for j = 1:length (matches)
            column = find(matchIndexes(end, :)==matches(1, j));
			if ~isempty(column)
                newLine(column(1)) = matches(2,j);
				newPointViewLine(:,column(1)) = f2(1:2,matches(2,j));
			end
        end
        matchIndexes = [matchIndexes, zeros(size(matchIndexes,1),size(newLine,2)-size(matchIndexes,2)) ; newLine];
        pointViewMatrix = [pointViewMatrix ; newPointViewLine];
    end
    
    pointViewMatrix = pointViewMatrix(:,min(pointViewMatrix)~=0);
    [U, W, V] = svd(pointViewMatrix);
    U3 = U(:,1:3);
    W3 = W(1:3,1:3);
    V3 = V(:,1:3);
    
    M = U3* (W3.^0.5);
    S = (W3.^0.5) * (V3');
    
    fscatter3( S(1,:), S(2,:), S(3,:), [0,10] );
    
end
