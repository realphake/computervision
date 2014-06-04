function [ out_pointcloud, R, T, iterations_stats ] = estimate_then_merge(useTheseNumbers, sampleSize, sampleTech, maxIterations)
    % ** Preferably use merge_pointclouds as the main function **
    % useTheseNumbers: a list of image numbers to use example: [0:1:99]
    % sampleSize: total amount of points sampled, except for normal space,
    % which indicates the samples per bin
    % sampleTech: 'none', 'uniform', 'random', 'normals'
    % maxIterations: set to 0 to not use the maximum iterations, otherwise
    % the ICP algorithm will only use the maxIterations for merging the two
    % pointclouds
    
    % initialize output variables

    % keeping track of how many iterations were necessary for each frame
    iterations_stats = [];

    T = zeros(1, 3, length(useTheseNumbers));
    R = zeros(3, 3, length(useTheseNumbers));
    R(:, :, 1) = eye(3);
    last_T = [0, 0, 0];
    last_R = eye(3);
    cutoffpoint_noise = 2;
    Base = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    % remove noise from the read point cloud
    if ( strcmpi( sampleTech, 'normals' ) )
        % if the technique involves normal space, then read in the normals
        % as well and remove the noise that is located in both data sets
        normals = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '_normal.pcd']);
        [ Base, normals ] = removeNoise(Base, cutoffpoint_noise, normals); 
    else
        Base = removeNoise(Base, cutoffpoint_noise);
    end
    Base = removeNoise(Base, cutoffpoint_noise); 
    out_pointcloud = Base;
    
    merges_done = 1;
    for i = useTheseNumbers(2:length(useTheseNumbers))

        disp('Opening pointcloud ', num2str(i));
        % read in the target point cloud
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        if ( strcmpi( sampleTech, 'normals' ) )
            % if the technique involves normal space, then read in the normals
            % as well and remove the noise that is located in both data sets
            next_normals = readPcd(['00000000', sprintf('%02d', i), '_normal.pcd']);
            [ Target, next_normals ] = removeNoise(Target, cutoffpoint_noise, next_normals);
            % get the rotation and translation of the base from ICP with
            % while using normal space uniform sampling
            [output, R_obtained, T_obtained, iterations] = ICP(Base, Target, sampleSize, sampleTech, maxIterations, normals);

            normals = next_normals;
        else
            % get the rotation and translation of the base from ICP
            Target = removeNoise(Target, cutoffpoint_noise); %target
            
            [output, R_obtained, T_obtained, iterations] = ICP(Base, Target, sampleSize, sampleTech, maxIterations);
        end
        % update the output variables

        iterations_stats = [iterations_stats; i iterations];

        out_pointcloud = [out_pointcloud;(last_R * output')' + repmat(last_T,length(output),1)];
        last_T = last_T + T_obtained * last_R';
        last_R = last_R * R_obtained;
        Base = Target;
        merges_done = merges_done+1;
        T(:, :, merges_done) = last_T;
        R(:, :, merges_done) = last_R;
    end
end