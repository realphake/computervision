function out_pointcloud = estimate_then_merge(useTheseNumbers, sampleSize, sampleTech)
    last_T = [0, 0, 0];
    last_R = eye(3);
    Base = readPcd(['00000000', sprintf('%02d', useTheseNumbers(1)), '.pcd']);
    Base = removeNoise(Base); %target
    Base = Base(:,1:3);
    out_pointcloud = Base;
    pause on
    for i = useTheseNumbers(2:length(useTheseNumbers))
        disp(num2str(i));
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        Target = removeNoise(Target); %target
        Target = Target(:,1:3);
        [output, R_obtained, T_obtained] = ICP(Base, Target, sampleSize, sampleTech);
        out_pointcloud = [out_pointcloud;(last_R * output')' + repmat(last_T,length(output),1)];
        last_T = last_T + T_obtained * last_R';
        last_R = last_R * R_obtained;
        Base = Target;
        %displayPointCloud(out_pointcloud);
    end
    pause off
end