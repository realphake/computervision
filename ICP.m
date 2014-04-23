function output = ICP( A1, A2, sampleSize, sampleTech )
    R = eye(3);
    A1 = removeNoise(A1); %base
    A2 = removeNoise(A2); %target
    disp('A1');
    A1 = subsampling(A1, sampleSize, sampleTech);
    disp('A2');
    A2 = subsampling(A2, sampleSize, sampleTech);

    tree = kd_buildtree(A2,0);
    vec_value_list = zeros(sampleSize,4);
    for i = 1:sampleSize
        point = A1(i,:);
        [index_vals,vec_vals,node_number] = kd_closestpointfast(tree,point);
        vec_value_list(i,:) = vec_vals;
    end
    monstrosity = [A1';vec_value_list'];
    SVD(monstrosity)
    
    % compute distances
    %disp('calculating distances...');
    %dd = distances(A1, A2);

    %disp('calculating the point in the target set that is nearest to each of the base set...');
    %[mindist target_index] = min(dd,[],2);

    %disp('all done.');
    %output = [mindist target_index];
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
