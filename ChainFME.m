function ChainFME( )

    [~,matches,f1,f2] = FME(sprintf('obj02_%03d.png',16), sprintf('obj02_%03d.png',1), '');
	pointViewMatrix = matches;
    
    for i = 1:15
        [~,matches,f1,f2] = FME(sprintf('obj02_%03d.png',i), sprintf('obj02_%03d.png',i+1), '');
		currentLength = length (pointViewMatrix) ;
		pointViewMatrix = [pointViewMatrix;zeros(size(pointViewMatrix,1))];
		for j = 1:length (matches)
			[contains,column] = ismemberf (pointViewMatrix(end, :), matches(1, i));
			if ~contains
				column = currentLength+1;
				currentLength += 1;
				pointViewMatrix = [pointViewMatrix,zeros(size(pointViewMatrix,2))];
				pointViewMatrix(end, column) = matches(2,j);
			else
				pointViewMatrix(end, column) = matches(2,j);
			end
			
		end
    end
    
    [U, W, V] = svd(pointViewMatrix);
    U3 = U(:,1:3);
    V3 = V(:,1:3);
    W3 = W(1:3,1:3);
    
    M = U3* (W3.^0.5);
    S = (W3.^0.5) * (V3');
    % OR:
    % M = U3
    % S = W3* (V3') 
    
    fscatter3( S(1,:), S(2,:), S(3,:), [0,10] );
    
end
