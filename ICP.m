function [T_total, R_total] = ICP( base, target, sampleSize, sampleTech )
    R = eye(3);
    T = [0, 0, 0];
    base = removeNoise(base); %base
    target = removeNoise(target); %target
    base = base(:,1:3);
    target = target(:,1:3);
    
    % selecting the first three columns that have the actual x y z
    % coordinates
    baseCloud = subsampling(base, sampleSize, sampleTech);
    % target is not sampled: the search is done in all its points
    % targetCloud = subsampling(target, sampleSize, sampleTech);
    IDX = knnsearch(target,baseCloud, 'NSMethod','kdtree');
    targetCloud = target(IDX,:);
    size(targetCloud);
    
    baseCloudCenter = mean(baseCloud,1);
    targetCloudCenter = mean(targetCloud,1);
    A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (targetCloud - repmat(targetCloudCenter, sampleSize,1));
    [U,S,V] = svd(A);
    R_total = U*R*V';
    R = U*V'
    T_total = T + ( baseCloudCenter - targetCloudCenter * R );
    T = baseCloudCenter - targetCloudCenter * R
    targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
    error = norm(mean(baseCloud-targetCloud,1))
    pause on
    
    while error > 0.0012
        target_new = (R * target')' + repmat(T,length(target),1);
        figure(1);
        scatter3(baseCloud(:,1),baseCloud(:,2),baseCloud(:,3), 3, [1, 0, 0]);
        hold on
        scatter3(targetCloud(:,1),targetCloud(:,2),targetCloud(:,3), 3, [0, 0.5, 0]);
        hold off
        pause(3);
        
        baseCloud = subsampling(base, sampleSize, sampleTech);
        IDX = knnsearch(target_new,baseCloud, 'NSMethod','kdtree');
        targetCloud = target(IDX,:);
        size(targetCloud);
        
        baseCloudCenter = mean(baseCloud,1);
        targetCloudCenter = mean(targetCloud,1);
        A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (targetCloud - repmat(targetCloudCenter, sampleSize,1));
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
        out = in(randsample(length(in), sampleSize), :);
    end
end

function uniform = uniformSubsampling(data, samples)
    stepsize = ceil(size(data,1)/samples);
    uniform = data( 1:stepsize:size(data,1), : );
end

function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end
