function [ out_pointcloud, R, T ] = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech)
    T = zeros(1, 3, length(useTheseNumbers));
    R = zeros(3, 3, length(useTheseNumbers));
    R(:, :, 1) = eye(3);
    out_pointcloud = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    out_pointcloud = removeNoise(out_pointcloud); %base
    out_pointcloud = out_pointcloud(:,1:3);
    merges_done = 1;
    for i = useTheseNumbers(2:length(useTheseNumbers))
        disp(num2str(i));
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        Target = removeNoise(Target); %target
        Target = Target(:,1:3);
        [output, R_obtained, T_obtained] = ICP(out_pointcloud, Target, sampleSize, sampleTech, normals);
        out_pointcloud = [out_pointcloud;output];
        merges_done = merges_done+1;
        T(:, :, merges_done) = T_obtained;
        R(:, :, 1) = R_obtained;
    end
end
