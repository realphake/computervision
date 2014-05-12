function [ out_pointcloud, R, T ] = estimate_then_merge(useTheseNumbers, sampleSize, sampleTech)
    T = zeros(1, 3, length(useTheseNumbers));
    R = zeros(3, 3, length(useTheseNumbers));
    R(:, :, 1) = eye(3);
    last_T = [0, 0, 0];
    last_R = eye(3);
    cutoffpoint_noise = 2;
    Base = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    if ( strcmpi( sampleTech, 'normals' ) )
        normals = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '_normal.pcd']);
        [ Base, normals ] = removeNoise(Base, cutoffpoint_noise, normals); 
    else
        Base = removeNoise(Base, cutoffpoint_noise);
    end
    Base = removeNoise(Base, cutoffpoint_noise); 
    out_pointcloud = Base;
    merges_done = 1;
    for i = useTheseNumbers(2:length(useTheseNumbers))
        disp(num2str(i));
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        if ( strcmpi( sampleTech, 'normals' ) )
            next_normals = readPcd(['00000000', sprintf('%02d', i), '_normal.pcd']);
            [ Target, next_normals ] = removeNoise(Target, cutoffpoint_noise, next_normals); %target
            [output, R_obtained, T_obtained] = ICP(Base, Target, sampleSize, sampleTech, normals);
            normals = next_normals;
        else
            Target = removeNoise(Target, cutoffpoint_noise); %target
            [output, R_obtained, T_obtained] = ICP(Base, Target, sampleSize, sampleTech);
        end
        out_pointcloud = [out_pointcloud;(last_R * output')' + repmat(last_T,length(output),1)];
        last_T = last_T + T_obtained * last_R';
        last_R = last_R * R_obtained;
        Base = Target;
        merges_done = merges_done+1;
        T(:, :, merges_done) = last_T;
        R(:, :, merges_done) = last_R;
    end
end