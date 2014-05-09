function out_pointcloud = merge_pointclouds(useTheseNumbers, sampleSize, sampleTech, type) 
    % useTheseNumbers = a list of image numbers to use example: [0:1:99]
    % sampleSize = total amount of points sampled, except for normal space,
    % which indicates the samples per bin
    % sampleTech = 'none', 'uniform', 'random', 'normals'
    % type = 'merge' which merges each consecutive frame into a new 
    % pointcloud and uses it as base or 'estimate' which finds the merge 
    % for consecutive frames and uses the individual transformation details
    % to create the final pointcloud.
    out_pointcloud = 0;
    if strcmpi( type, 'merge')
        out_pointcloud = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech);
    elseif strcmpi( type, 'estimate')
        out_pointcloud = estimate_then_merge(useTheseNumbers, sampleSize, sampleTech);
    else
        disp('Type unknown: exiting');
    end
    displayPointCloud(out_pointcloud);
end

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

function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end

function displayPointCloud(out_pointcloud)
    stepsize = ceil(size(out_pointcloud,1)/20000);
    uniform = out_pointcloud( 1:stepsize:size(out_pointcloud,1), : );
    % fscatter was created by Felix Morsdorf for which we take no credit
    fscatter3(uniform(:,1),uniform(:,2),uniform(:,3), [0 10]);
end
