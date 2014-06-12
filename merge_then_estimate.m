function[ out_pointcloud, R, T, iterations_stats ] = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech, maxIterations)


    % ** Preferably use merge_pointclouds as the main function **
    % useTheseNumbers: a list of image numbers to use example: [0:1:99]
    % sampleSize: total amount of points sampled, except for normal space,
    % which indicates the samples per bin
    % sampleTech: 'none', 'uniform', 'random', 'normals'
    % maxIterations: set to 0 to not use the maximum iterations, otherwise
    % the ICP algorithm will only use the maxIterations for merging the two
    % pointclouds
   
    % keeping track of how many iterations were necessary for each frame
    iterations_stats = [];

    % initialize output variables

    T = zeros(1, 3, length(useTheseNumbers));
    R = zeros(3, 3, length(useTheseNumbers));
    R(:, :, 1) = eye(3);
    cutoffpoint_noise = 2;
    out_pointcloud = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    % remove noise from the read point cloud

    if ( strcmpi( sampleTech, 'normals' ) )
        % if the technique involves normal space, then read in the normals
        % as well and remove the noise that is located in both data sets
        normals = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '_normal.pcd']);
        [out_pointcloud, normals] = removeNoise(out_pointcloud, cutoffpoint_noise, normals); %base
    else
        out_pointcloud = removeNoise(out_pointcloud, cutoffpoint_noise); %base
    end
    
    
    merges_done = 1;
    for i = useTheseNumbers(2:length(useTheseNumbers))
		disp(['opening pointcloud: ', num2str(i)]);
        % read in the target point cloud
        Base = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        if ( strcmpi( sampleTech, 'normals' ) )
            % if the technique involves normal space, then read in the normals
            % as well and remove the noise that is located in both data sets
            next_normals = readPcd(['00000000', sprintf('%02d', i), '_normal.pcd']);
            [ Base, next_normals ] = removeNoise(Base, cutoffpoint_noise, next_normals);
            % get the rotation and translation of the base from ICP with
            % while using normal space uniform sampling
            [output, R_obtained, T_obtained, iterations] = ICP(Base, out_pointcloud, sampleSize, sampleTech, maxIterations, normals);
            normals = [normals;(R_obtained * next_normals')'];
        else
            % get the rotation and translation of the base from ICP
            Base = removeNoise(Base, cutoffpoint_noise);
            [output, R_obtained, T_obtained, iterations] = ICP(Base, out_pointcloud, sampleSize, sampleTech, maxIterations);
        end
        
        iterations_stats = [iterations_stats; i iterations];

        % update the output variables
        out_pointcloud = [Base;output];

        merges_done = merges_done+1;
        T(:, :, merges_done) = T_obtained;
        R(:, :, merges_done) = R_obtained;
    end
end
