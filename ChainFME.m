function ChainFME( )

    [~,matches,f1,f2] = FME(sprintf('obj02_%03d.png',16), sprintf('obj02_%03d.png',1), '');
    pointViewMatrix = zeros( 2, matches(1,length(matches)) );
    pointViewMatrix(1,matches(1,:)) = f2 ( 1, matches(2,:) ) ;
    pointViewMatrix(2,matches(1,:)) = f2 ( 2, matches(2,:) ) ;
    
    for i = 1:15
        [~,matches,f1,f2] = FME(sprintf('obj02_%03d.png',i), sprintf('obj02_%03d.png',i+1), '');
        if matches(1,length(matches)) > length(pointViewMatrix)
            temp = zeros( i*2+2, matches(1,length(matches)) );
        else
            temp = zeros( i*2+2, length(pointViewMatrix) );
        end
        temp(1:size(pointViewMatrix, 1), 1:size(pointViewMatrix, 2)) = pointViewMatrix;
        temp(i*2+1, matches(1,:)) = f2 ( 1, matches(2,:) ) ;
        temp(i*2+2, matches(1,:)) = f2 ( 2, matches(2,:) ) ;
        pointViewMatrix = temp;
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
