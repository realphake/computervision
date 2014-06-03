function [M, S, pointViewMatrix, pvm] = ChainFME( )

    [~,matches,f1,f2] = FME(sprintf('obj02_%03d.png',16), sprintf('obj02_%03d.png',1), 'ransac');
	f1m(1,:) = f1(1,:) - mean(f1(1,:));
    f1m(2,:) = f1(2,:) - mean(f1(2,:));
    f2m(1,:) = f2(1,:) - mean(f2(1,:));
    f2m(2,:) = f2(2,:) - mean(f2(2,:));
    
    matchIndexes = matches;
    pointViewMatrix(1,:) = f1m(1,matchIndexes(1,:));
    pointViewMatrix(2,:) = f1m(2,matchIndexes(1,:));
    pointViewMatrix(3,:) = f2m(1,matchIndexes(2,:));
    pointViewMatrix(4,:) = f2m(2,matchIndexes(2,:));
    
    pvm(1,:) = f1(1,matchIndexes(1,:));
    pvm(2,:) = f1(2,matchIndexes(1,:));
    pvm(3,:) = f2(1,matchIndexes(2,:));
    pvm(4,:) = f2(2,matchIndexes(2,:));
    
    clear f2m f1m
    
    for i = 1:15
        [~,matches,~,f2] = FME(sprintf('obj02_%03d.png',i), sprintf('obj02_%03d.png',i+1), 'ransac');
        f2m(1,:) = f2(1,:) - mean(f2(1,:));
        f2m(2,:) = f2(2,:) - mean(f2(2,:));
        
        newLine = zeros(1,size(matchIndexes,2));
        newPointViewLine = zeros(2,size(matchIndexes,2));
        pvL = zeros(2,size(matchIndexes,2));
        for j = 1:length (matches)
            column = find(matchIndexes(end, :)==matches(1, j));
			if isempty(column)
				newLine = [newLine,matches(2,j)];
                newPointViewLine = [newPointViewLine,[f2m(1,matches(2,j));f2m(2,matches(2,j))]];
                pvL = [pvL,[f2(1,matches(2,j));f2(2,matches(2,j))]];
			else
				newLine(column(1)) = matches(2,j);
				newPointViewLine(:,column(1)) = f2m(1:2,matches(2,j));
                pvL(:,column(1)) = f2(1:2,matches(2,j));
			end
        end
        matchIndexes = [matchIndexes, zeros(size(matchIndexes,1),size(newLine,2)-size(matchIndexes,2)) ; newLine];
        pointViewMatrix = [pointViewMatrix, zeros(size(pointViewMatrix,1),size(newPointViewLine,2)-size(pointViewMatrix,2)) ; newPointViewLine];
        pvm = [pvm, zeros(size(pvm,1),size(pvL,2)-size(pvm,2)) ; pvL];
        
        clear f2m f1m
    end
    
    [U, W, V] = svd(pointViewMatrix);
    U3 = U(:,1:3);
    W3 = W(1:3,1:3);
    V3 = V(:,1:3);
    
    M = U3* (W3.^0.5);
    S = (W3.^0.5) * (V3');
    
    fscatter3( S(1,:), S(2,:), S(3,:), [0,10] );
    
end
