function output = ICP( A1, A2, sampleSize, sampleTech )
    R = eye(3);
    A1 = removeNoise(A1); %base
    A2 = removeNoise(A2); %target
    disp('A1');
    size(A1)
    A1 = subsampling(A1, sampleSize, sampleTech);
    disp('A2');
    size(A2)
    A2 = subsampling(A2, sampleSize, sampleTech);

    % compute distances
    disp('calculating distances...');
    dd = distances(A1, A2);

    disp('calculating the point in the target set that is nearest to each of the base set...');
    [mindist target_index] = min(dd,[],2);

    disp('all done.');
    output = [mindist target_index];
end

function ret = flatten(in)
    n_rows = size(in,1);
    n_cols = size(in,2);
    n_cells = n_rows * n_cols;
    ret = zeros(n_cells, 3);
    i = 0;
    for r = 1:n_rows
        for c = 1:n_cols
            i = i + 1;
            ret(i,1) = r;
            ret(i,2) = c;
            ret(i,3) = in(r,c);
        end
    end
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
