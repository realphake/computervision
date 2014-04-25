function [T_total R_total] = ICP( A1, A2, sampleSize, sampleTech )
    R = eye(3);
    T = [0, 0, 0];
    A1 = removeNoise(A1); %base
    A2 = removeNoise(A2); %target
    A1 = A1(:,1:3);
    A2 = A2(:,1:3);
    
    baseCloud = subsampling(A1, sampleSize, sampleTech);
    targetCloud = subsampling(A2, sampleSize, sampleTech);

    tree = kd_buildtree(targetCloud,0);
    vec_value_list = zeros(sampleSize,3);
    for i = 1:sampleSize
        point = baseCloud(i,:);
        [index_vals,vec_vals,node_number] = kd_closestpointfast(tree,point);
        vec_value_list(i,:) = vec_vals;
    end

    baseCloudCenter = mean(baseCloud,1);
    targetCloudCenter = mean(vec_value_list,1);
    A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (vec_value_list - repmat(targetCloudCenter, sampleSize,1));
    [U,S,V] = svd(A);
    R_total = U*R*V';
    R = U*V'
    T_total = T + ( baseCloudCenter - targetCloudCenter * R );
    T = baseCloudCenter - targetCloudCenter * R
    targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
    A2_new = zeros(size(A2));
    A2_center = mean(A2,1);
    error = norm(mean(baseCloud-targetCloud,1));
    pause on
    while error > 0.0012
        A2_new = (R * A2')' + repmat(T,length(A2),1);
        figure(1);
        scatter3(baseCloud(:,1),baseCloud(:,2),baseCloud(:,3), 3, [1, 0, 0]);
        scatter3(targetCloud(:,1),targetCloud(:,2),targetCloud(:,3), 3, [0, 0.5, 0]);
        drawnow();
        baseCloud = subsampling(A1, sampleSize, sampleTech);
        targetCloud = subsampling(A2_new, sampleSize, sampleTech);

        tree = kd_buildtree(targetCloud,0);
        vec_value_list = zeros(sampleSize,3);
        for i = 1:sampleSize
            point = baseCloud(i,:);
            [index_vals,vec_vals,node_number] = kd_closestpointfast(tree,point);
            vec_value_list(i,:) = vec_vals;
        end

        baseCloudCenter = mean(baseCloud,1);
        targetCloudCenter = mean(vec_value_list,1);
        A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (vec_value_list - repmat(targetCloudCenter, sampleSize,1));
        [U,S,V] = svd(A);
        R_total = U*R*V'
        R = U*V';
        T_total = T + ( baseCloudCenter - targetCloudCenter * R )
        T = baseCloudCenter - targetCloudCenter * R;
        targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
        error = norm(mean(baseCloud-targetCloud,1))
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
