function pointViewMatrix = ChainFME( )

    [~,matches] = FME(sprintf('obj02_%03d.png',16), sprintf('obj02_%03d.png',1), '');
    pointViewMatrix = zeros( 1, matches(1,length(matches)) );
    pointViewMatrix(matches(1,:)) = matches(2,:);
    
    for i = 1:15
        [~,matches] = FME(sprintf('obj02_%03d.png',i), sprintf('obj02_%03d.png',i+1), '');
        if matches(1,length(matches)) > length(pointViewMatrix)
            temp = zeros( i+1, matches(1,length(matches)) );
        else
            temp = zeros( i+1, length(pointViewMatrix) );
        end
        temp(1:size(pointViewMatrix, 1), 1:size(pointViewMatrix, 2)) = pointViewMatrix;
        temp(i+1, matches(1,:)) = matches(2,:);
        pointViewMatrix = temp;
    end

end

