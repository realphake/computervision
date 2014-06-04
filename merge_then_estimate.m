function[ out_pointcloud, R, T, iterations_stats ] = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech)

    % keeping track of how many iterations were necessary for each frame
    iterations_stats = [];

    T = zeros(1, 3, length(useTheseNumbers));
    R = zeros(3, 3, length(useTheseNumbers));
    R(:, :, 1) = eye(3);
    cutoffpoint_noise = 2;
    out_pointcloud = readPcd(['data/00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    if ( strcmpi( sampleTech, 'normals' ) )
        normals = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '_normal.pcd']);
        [out_pointcloud, normals] = removeNoise(out_pointcloud, cutoffpoint_noise, normals); %base
    else
        out_pointcloud = removeNoise(out_pointcloud, cutoffpoint_noise); %base
    end
    merges_done = 1;
    for i = useTheseNumbers(2:length(useTheseNumbers))
        disp(num2str(i));
        Target = readPcd(['data/00000000', sprintf('%02d', i), '.pcd']);
        if ( strcmpi( sampleTech, 'normals' ) )
            next_normals = readPcd(['data/00000000', sprintf('%02d', i), '_normal.pcd']);
            [ Target, next_normals ] = removeNoise(Target, cutoffpoint_noise, next_normals); %target
            [output, R_obtained, T_obtained, iterations] = ICP(out_pointcloud, Target, sampleSize, sampleTech, normals);
            normals = [normals;(R_obtained * next_normals')'];
        else
            Target = removeNoise(Target, cutoffpoint_noise); %target
            [output, R_obtained, T_obtained, iterations] = ICP(out_pointcloud, Target, sampleSize, sampleTech);
        end
        
        iterations_stats = [iterations_stats; i iterations];

        out_pointcloud = [out_pointcloud;output];
        merges_done = merges_done+1;
        T(:, :, merges_done) = T_obtained;
        R(:, :, merges_done) = R_obtained;
    end
end
