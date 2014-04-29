function [T_total R_total] = ICP( source, target, sampleSize, sampleTech )
    R = eye(3);
    T = [0, 0, 0];
    source = removeNoise(source); % will be sampled and transformed
    target = removeNoise(target); % will not be transformed or sample
    
    % selecting the first three columns that have the actual x y z
    % coordinates
    source = source(:,1:3);
    target = target(:,1:3);
    
    sourceCloud = subsampling(source, sampleSize, sampleTech);
    targetCloud = target; % target is not sampled: the search is done in all its points
    %targetCloud = subsampling(target, sampleSize, sampleTech);

    tree = kd_buildtree(targetCloud,0);
    vec_value_list = zeros(sampleSize,3);
    for i = 1:sampleSize
        point = sourceCloud(i,:);
        [index_vals,vec_vals,node_number] = kd_closestpointfast(tree,point);
        vec_value_list(i,:) = vec_vals;
    end

    sourceCloudCenter = mean(sourceCloud,1);
    targetCloudCenter = mean(vec_value_list,1);
    A = (sourceCloud - repmat(sourceCloudCenter, sampleSize,1))' * (vec_value_list - repmat(targetCloudCenter, sampleSize,1));
    [U,S,V] = svd(A);
    R_total = U*R*V';
    R = U*V'
    T_total = T + ( sourceCloudCenter - targetCloudCenter * R );
    T = sourceCloudCenter - targetCloudCenter * R
    targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
    target_new = zeros(size(target));
    target_center = mean(target,1);
    error = norm(mean(sourceCloud-targetCloud,1));
    pause on
    while error > 0.0012
        target_new = (R * target')' + repmat(T,length(target),1);
        figure(1);
        scatter3(sourceCloud(:,1),sourceCloud(:,2),sourceCloud(:,3), 3, [1, 0, 0]);
        scatter3(targetCloud(:,1),targetCloud(:,2),targetCloud(:,3), 3, [0, 0.5, 0]);
        drawnow();
        sourceCloud = subsampling(source, sampleSize, sampleTech);
        targetCloud = subsampling(target_new, sampleSize, sampleTech);

        tree = kd_buildtree(targetCloud,0);
        vec_value_list = zeros(sampleSize,3);
        for i = 1:sampleSize
            point = sourceCloud(i,:);
            [index_vals,vec_vals,node_number] = kd_closestpointfast(tree,point);
            vec_value_list(i,:) = vec_vals;
        end

        sourceCloudCenter = mean(sourceCloud,1);
        targetCloudCenter = mean(vec_value_list,1);
        A = (sourceCloud - repmat(sourceCloudCenter, sampleSize,1))' * (vec_value_list - repmat(targetCloudCenter, sampleSize,1));
        [U,S,V] = svd(A);
        R_total = U*R*V'
        R = U*V';
        T_total = T + ( sourceCloudCenter - targetCloudCenter * R )
        T = sourceCloudCenter - targetCloudCenter * R;
        targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
        error = norm(mean(sourceCloud-targetCloud,1))
    end
    pause off
end

function out = subsampling(in, sampleSize, technique)
    if ( strcmpi( technique, 'uniform' ) )
        out = uniformSubsampling(in, sampleSize);
    end
    if ( strcmpi( technique, 'random' ) )
        out = randsample(in, sampleSize);
    end
end

function uniform = uniformSubsampling(data, samples)
    stepsize = ceil(size(data,1)/samples);
    uniform = data( 1:stepsize:size(data,1), : );
end

function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end
