function out_pointcloud = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech)
    last_T = [0, 0, 0];
    last_R = eye(3);
    out_pointcloud = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    out_pointcloud = removeNoise(out_pointcloud); %base
    out_pointcloud = out_pointcloud(:,1:3);
    for i = useTheseNumbers(2:length(useTheseNumbers))
        disp(num2str(i));
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        Target = removeNoise(Target); %target
        Target = Target(:,1:3);
        [output, R_obtained, T_obtained] = ICP(out_pointcloud, Target, sampleSize, sampleTech);
        out_pointcloud = [out_pointcloud;(last_R * output')' + repmat(last_T,length(output),1)];
        last_T = last_T + T_obtained * last_R';
        last_R = last_R * R_obtained;
    end
end
