% This function makes a 3d pointcloud out of a string of images.
% Specifically, this code uses the House dataset. This code can be run by
% placing it in the same directory as the images or by adding the image
% folder to the MATLAB working directory. Also make sure FME.m is in the
% same folder, and the vlfeat-0.9.18 folder. Then just run the command as
% below, [M,S, pointViewMatrix] = ChainFME( ).
function [M,S, pointViewMatrix] = ChainFME( )
    % We start off by initializing the pointViewMatrix we'll be making,
    % using real data because otherwise we'd need to jump through weird
    % hoops to include data from the first image.
    [~,matches,f1,f2] = FME(sprintf('frame00000%03d.png',1), sprintf('frame00000%03d.png',2), 'ransac');
    % Matchindexes will be filled with the contents of "matches". These are
    % meaningless numbers without the feature matrices whose column indices
    % they contain, but they have a convenient quality: they can be checked
    % for equality easily. The matrix will be discarded after the loop
    % terminates.
    matchIndexes = matches;
    % The columns of pointViewMatrix will contain all the locations of a
    % single point across the images. The rows will alternately be x and y
    % coordinates in the images in sequence. 
    pointViewMatrix(1,:) = f1(1,matchIndexes(1,:));
    pointViewMatrix(2,:) = f1(2,matchIndexes(1,:));
    pointViewMatrix(3,:) = f2(1,matchIndexes(2,:));
    pointViewMatrix(4,:) = f2(2,matchIndexes(2,:));
    
    % Stop one short of the full number (49) of images, because each
    % iteration already looks ahead by one.
    for i = 2:48
        % Same deal as before, but only look at the second of these images.
        % The first on is there for comparison, but we've already put it in
        % the matrix in the previous iteration.
        [~,matches,~,f2] = FME(sprintf('frame00000%03d.png',i), sprintf('frame00000%03d.png',i+1), 'ransac');
        % Originally, we also centered the features, but it gave strange
        % results. It really doesn't add too much anyway.
        newLine = zeros(1,size(matchIndexes,2)); % Create a new line.
        newPointViewLine = zeros(2,size(matchIndexes,2)); % In both matrices.
        for j = 1:length (matches)
            % Now we're going to add the current feature to the matrix, but
            % only if we've seen it before and of course in the same column
            % as where it was added before. That's where we use the
            % matchIndexes matrix.
            column = find(matchIndexes(end, :)==matches(1, j));
			if ~isempty(column)
                newLine(column(1)) = matches(2,j);
				newPointViewLine(:,column(1)) = f2(1:2,matches(2,j));
            end
            % You might notice a lot of found matches are being thrown
            % away. That's to make sure they aren't duplicates. I'll go
            % into this in a moment.
        end
        % The lines that were just made are now added below their matrices.
        matchIndexes = [matchIndexes, zeros(size(matchIndexes,1),size(newLine,2)-size(matchIndexes,2)) ; newLine];
        pointViewMatrix = [pointViewMatrix ; newPointViewLine];
    end
    
    % A single feature can't be presented in mutliple columns, because that
    % will ruin the pointcloud presentation. Because of that, I didn't add
    % any new features that weren't seen in the first image. In this line,
    % I remove any feature that was lost at some point. That is pretty
    % drastic, but for some reason the result still isn't correct if I
    % don't do this.
    pointViewMatrix = pointViewMatrix(:,min(abs(pointViewMatrix))~=0);
    
    % To create the pointcloud, decompose the pointViewMatrix and
    % concatenate their elements so they reconstruct to a 3xN matrix. Then
    % reconstruct only half of the original pointViewMatrix for the motion
    % data M and the other half for the shape data S.
    [U, W, V] = svd(pointViewMatrix);
    U3 = U(:,1:3);
    W3 = W(1:3,1:3);
    V3 = V(:,1:3);
    M = U3* (W3.^0.5);
    S = (W3.^0.5) * (V3');
    
    % Then print the structure to the screen!
    fscatter3( S(1,:), S(2,:), S(3,:), [0,10] );
    
end
