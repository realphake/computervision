function [output, R_total, T_total] = ICP( base, target, sampleSize, sampleTech )
    R = eye(3);
    R_total = eye(3);
    T = [0, 0, 0];
    T_total = [0, 0, 0];

    
    % selecting the first three columns that have the actual x y z
    % coordinates
    baseCloud = subsampling(base, sampleSize, sampleTech);
    % target is not sampled: the search is done in all its points
    % targetCloud = subsampling(target, sampleSize, sampleTech);
    IDX = knnsearch(target,baseCloud, 'NSMethod','kdtree');
    targetCloud = target(IDX,:);
    
    baseCloudCenter = mean(baseCloud,1);
    targetCloudCenter = mean(targetCloud,1);

    A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (targetCloud - repmat(targetCloudCenter, sampleSize,1));
    [U,S,V] = svd(A);
    R_total = R_total*(U*V');
    R = U*V';
    T = baseCloudCenter - (targetCloudCenter * R');
    T_total = T_total + T;
    targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
    old_error = -1;
    error = norm(mean(baseCloud-targetCloud,1));
    pause on
    target_new = target;
    iterations = 0;
    while error ~= old_error && iterations ~= 20 
        target_new = (R * target_new')' + repmat(T,length(target),1);
        
        if ( strcmpi( sampleTech, 'random' ) )
            baseCloud = subsampling(base, sampleSize, sampleTech);
        end
        IDX = knnsearch(target_new,baseCloud, 'NSMethod','kdtree');
        targetCloud = target_new(IDX,:);
        
        baseCloudCenter = mean(baseCloud,1);
        targetCloudCenter = mean(targetCloud,1);
        A = (baseCloud - repmat(baseCloudCenter, sampleSize,1))' * (targetCloud - repmat(targetCloudCenter, sampleSize,1));
        [U,S,V] = svd(A);
        R_total = R_total*(U*V');
        R = U*V';
        T = baseCloudCenter - (targetCloudCenter * R');
        T_total = T_total + T;
        targetCloud = (R * targetCloud')' + repmat(T,sampleSize,1);
        old_error = error;
        error = rms(baseCloud-targetCloud,1);
%       [x, y, z] = decompose_rotation(R_total);
        iterations = iterations + 1;
        %displayPointClouds(baseCloud, targetCloud);
        %pause(2);
    end
    output = (R * target_new')' + repmat(T,length(target),1);
    pause off
end

function out = subsampling(in, sampleSize, technique)
    if ( strcmpi( technique, 'uniform' ) )
        out = uniformSubsampling(in, sampleSize);
    end
    if ( strcmpi( technique, 'random' ) )
        out = in(randsample(length(in), sampleSize), :);
    end
    if ( strcmpi( technique, 'normals' ) )
        out = sampleNormalSpace(in, sampleSize);
    end
    if ( strcmpi( technique, 'none' ) )
        out = in;
    end
end

function uniform = uniformSubsampling(data, samples)
    stepsize = ceil(size(data,1)/samples);
    uniform = data( 1:stepsize:size(data,1), : );
end

% decompose copied over from http://nghiaho.com/uploads/code/rotation_matrix_demo.m
function [x,y,z] = decompose_rotation(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end

% function noNoise = removeNoise(noise)
%     noNoise = noise(noise(:,3) < 2,:);
% end

function displayPointClouds(PointCloud1, PointCloud2)
    figure(1);
    clf();
    hold on
    stepsize = ceil(size(PointCloud1,1)/5000);
    uniform = PointCloud1( 1:stepsize:size(PointCloud1,1), : );
    scatter3(uniform(:,1),uniform(:,2),uniform(:,3), 3, [1, 0, 0]);
    stepsize = ceil(size(PointCloud2,1)/5000);
    uniform = PointCloud2( 1:stepsize:size(PointCloud2,1), : );
    scatter3(uniform(:,1),uniform(:,2),uniform(:,3), 3, [0, 0.5, 0]);
    hold off
end

function out = sampleNormalSpace(in, sampleSize)
    mesh = pointCloud2mesh(in);
    mesh.vertexNormals
    icosahedron = createIcosahedron(); % 20 sides as bins for the normals
    TRI = TriRep(icosahedron.faces, icosahedron.vertices(:,1),icosahedron.vertices(:,2),icosahedron.vertices(:,3));
    faceN
    
end

