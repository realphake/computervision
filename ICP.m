function [output, R_total, T_total] = ICP( base, target, sampleSize, sampleTech )

    % Expects base to be of size N by 3, target to be of size M by 3
    % Sample techniques:
    % None - all points are used
    % Uniform - Samplesize of points is taken from the base with a certain
    %           stepsize
    % Random - Samplesize of points randomly taken from the base
    % Normals - Samplesize taken uniformly per bin, normals are divided 
    %           into 20 bins

    R = eye(3);
    R_total = eye(3);
    T = [0, 0, 0];
    T_total = [0, 0, 0];
    
    % selecting the first three columns that have the actual x y z
    % coordinates
    baseCloud = subsampling(base, sampleSize, sampleTech);
    base_length = length(baseCloud);
    % target is not sampled: the search is done in all its points
    IDX = knnsearch(target,baseCloud, 'NSMethod','kdtree');
    targetCloud = target(IDX,:);
    
    baseCloudCenter = mean(baseCloud,1);
    targetCloudCenter = mean(targetCloud,1);

    A = (baseCloud - repmat(baseCloudCenter, base_length,1))' * (targetCloud - repmat(targetCloudCenter, base_length,1));
    [U,S,V] = svd(A);
    R_total = R_total*(U*V');
    R = U*V';
    T = baseCloudCenter - (targetCloudCenter * R');
    T_total = T_total + T;
    targetCloud = (R * targetCloud')' + repmat(T,base_length,1);
    old_error = -1;
    error = pdist([rms(baseCloud-targetCloud,1); 0 0 0], 'euclidean');
    pause on
    target_new = target;
    iterations = 0;
    while error ~= old_error && iterations ~= 50 
        target_new = (R * target_new')' + repmat(T,length(target),1);
        
        if strcmpi( sampleTech, 'random' ) 
            baseCloud = subsampling(base, sampleSize, sampleTech);
        end
        IDX = knnsearch(target_new,baseCloud, 'NSMethod','kdtree');
        targetCloud = target_new(IDX,:);
        
        baseCloudCenter = mean(baseCloud,1);
        targetCloudCenter = mean(targetCloud,1);
        A = (baseCloud - repmat(baseCloudCenter, base_length,1))' * (targetCloud - repmat(targetCloudCenter, base_length,1));
        [U,S,V] = svd(A);
        R_total = R_total*(U*V');
        R = U*V';
        T = baseCloudCenter - (targetCloudCenter * R');
        T_total = T_total + T;
        targetCloud = (R * targetCloud')' + repmat(T,base_length,1);
        old_error = error;
        error = pdist([rms(baseCloud-targetCloud,1); 0 0 0], 'euclidean');
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
    % pointCloud2mesh is created by Ajmal Saeed Mian
    % see specific .m file for license and credits
    mesh = pointCloud2mesh(in);
    % createIcosahedron is created by David Legland
    % see specific .m file for license and credits
    icosahedron = createIcosahedron(); % 20 sides as bins for the normals
    normals_icosahedron = surfaceNormals( icosahedron );
    similarity = mesh.vertexNormals * normals_icosahedron';
    [values, binning] = max(similarity, [], 2);
    out = [];
    for i = 1:20
        out = [ out ; uniformSubsampling(in(binning == i, :), sampleSize)];
    end
end

% copied from pointCloud2mesh.m and adjusted to only get the surface
% normals
% see specific .m file for license and credits
function normals_out = surfaceNormals( mesh )
    normals_out = zeros(length(mesh.faces), 3);
    for ii = 1:length(mesh.faces) %find normals of all triangles
        %indices of the points from which the triangle is made
        pointIndex1 = mesh.faces(ii,1);
        pointIndex2 = mesh.faces(ii,2);
        pointIndex3 = mesh.faces(ii,3);

        %coordinates of the points
        point1 = mesh.vertices(pointIndex1,:);
        point2 = mesh.vertices(pointIndex2,:);
        point3 = mesh.vertices(pointIndex3,:);

        vector1 = point2 - point1;
        vector2 = point3 - point2;

        normal = cross(vector1,vector2);
        normal = normal / norm(normal);
        
        % create_icosahedron is already set up to have the normals pointing
        % outward

        normals_out(ii,:)=normal;      
    end
end

