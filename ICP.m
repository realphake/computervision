function output = ICP( A1, A2, sampleTech )
    R = eye(3);
    A1 = removeNoise(A1);
    A2 = removeNoise(A2);
    A1 = subsampling(A1, sampleTech);
    A2 = subsampling(A2, sampleTech);
    size(A1)
end

function out = subsampling(in, technique)
    if ( strcmpi( technique, 'uniform' ) )
        out = uniformSubsampling(in, 100); end
    if ( strcmpi( technique, 'random' ) )
        out = randsample(in, 100); end
end

function uniform = uniformSubsampling(data, samples)
    stepsize = ceil(size(data,1)/samples);
    uniform = data( 1:stepsize:size(data,1), : );
end

function noNoise = removeNoise(noise)
    noNoise = noise(noise(:,3) < 2,:);
end