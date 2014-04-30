function out_pointcloud = merge_pointclouds(stepSize, sampleSize, sampleTech, type) 
    out_pointcloud = 0;
    if strcmpi( type, 'merge')
        out_pointcloud = merge_then_estimate(stepSize, sampleSize, sampleTech);
    elseif strcmpi( type, 'estimate')
        out_pointcloud = estimate_then_merge(stepSize, sampleSize, sampleTech);
    else
        disp('Type unknown: exiting');
    end
    stepsize = ceil(size(out_pointcloud,1)/5000);
    uniform = out_pointcloud( 1:stepsize:size(out_pointcloud,1), : );
    scatter3(uniform(:,1),uniform(:,2),uniform(:,3), 3, [1, 0, 0]);
end

function out_pointcloud = merge_then_estimate(stepSize, sampleSize, sampleTech)
    T_total_new = [0,0,0];
    R_total_new = eye(3);
    out_pointcloud = readPcd(['0000000000.pcd']);
    out_pointcloud = removeNoise(out_pointcloud); %base
    out_pointcloud = out_pointcloud(:,1:3);
    for i = 1:stepSize:100-stepSize
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        Target = removeNoise(Target); %target
        Target = Target(:,1:3);
        [T_obtained, R_obtained] = ICP(out_pointcloud, Target, sampleSize, sampleTech);
        T_total_new = T_total_new + T_obtained;
        R_total_new = R_total_new * R_obtained;
        out_pointcloud = [out_pointcloud;(R_total_new * Target')' + repmat(T_total_new,length(Target),1)];
    end
end

function out_pointcloud = estimate_then_merge(stepSize, sampleSize, sampleTech)
    T_total_new = [0,0,0];
    R_total_new = eye(3);
    out_pointcloud = 0;
    Base = readPcd('0000000000.pcd');
    Base = removeNoise(Base); %target
    Base = Base(:,1:3);
    for i = stepSize:stepSize:100-stepSize
        disp(num2str(i));
        Target = readPcd(['00000000', sprintf('%02d', i), '.pcd']);
        Target = removeNoise(Target); %target
        Target = Target(:,1:3);
        [T_obtained, R_obtained] = ICP(Base, Target, sampleSize, sampleTech);
        %T_total_new = T_total_new + T_obtained * R_total_new';
        %R_total_new = R_total_new * R_obtained;
        if i == 1
            out_pointcloud = [Base;(R_obtained * Target')' + repmat(T_obtained,length(Target),1)];
        else
            out_pointcloud = [out_pointcloud;(R_obtained * Target')' + repmat(T_obtained,length(Target),1)];
        end
        Base = Target;
    end
end

function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end
