function output = ICP( A1, A2, sampleSize, sampleTech )
    R = eye(3);
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
    monstrosity = [baseCloud';vec_value_list'];
    baseCloudCenter = mean(baseCloud,1)
    targetCloudCenter = mean(vec_value_list,1)
    [U,S,V] = svd(monstrosity);
    Vt = V';
    U = U(:,1:3)
    S = S(1:3,1:3)
    Vt = Vt(1:3,:)
    R
    R = Vt*R'*U;
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
