function [output, R_total, T_total, iterations] = ICP( base, target, sampleSize, sampleTech, targetNormals)

    % Expects base to be of size N by 3, target to be of size M by 3
    % Sample techniques:
    % None - all points are used
    % Uniform - Samplesize of points is taken from the base with a certain
    %           stepsize
    % Random - Samplesize of points randomly taken from the base
    % Normals - Samplesize taken uniformly per bin, normals are divided 
    %           into 20 bins, total number of samples is at most 20 * sampleSize
    %           as some bins might have datapoints < sampleSize

    R = eye(3);
    R_total = eye(3);
    T = [0, 0, 0];
    T_total = [0, 0, 0];
    if ( strcmpi( sampleTech, 'normals' ) ) && nargin == 5
        baseCloud = subsampling(base, sampleSize, sampleTech, targetNormals);
    else
        baseCloud = subsampling(base, sampleSize, sampleTech);
    end
    base_length = length(baseCloud);
    % target is not sampled: the search is done in all its points
    
    error = 0;
    pause on
    target_new = target;
    iterations = 0;
    epsilon_err = 0.000000000000001;
    good_enough = false;
    while ~ good_enough
        %tmp = [base;target_new];
        %displayPointCloud(tmp);
        %input('continue?');
        target_new = (R * target_new')' + repmat(T,length(target),1);
        
        % only needs resampling if random, the other methods are deterministic
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
        %[x, y, z] = decompose_rotation(R_total);
        iterations = iterations + 1
        %displayPointClouds(baseCloud, targetCloud);
        %pause(2);
        diff_err_normalized = abs(( error - old_error ) / old_error)
        good_enough = diff_err_normalized < epsilon_err; %|| iterations == 50;
    end
    output = (R * target_new')' + repmat(T,length(target),1);
    pause off
end

function out = subsampling(in, sampleSize, technique, normals)
    if ( strcmpi( technique, 'uniform' ) )
        out = uniformSubsampling(in, sampleSize);
    end
    if ( strcmpi( technique, 'random' ) )
        out = in(randsample(length(in), sampleSize), :);
    end
    if ( strcmpi( technique, 'normals' ) ) && nargin == 4
        out = sampleNormalSpace(in, sampleSize, normals);
    end
    if ( strcmpi( technique, 'normals' ) ) && nargin == 3
        out = sampleNormalSpace(in, sampleSize);
    end
    if ( strcmpi( technique, 'none' ) )
        out = in;
    end
end

function uniform = uniformSubsampling(data, samples)
    if length(data) < samples
        uniform = data;
    else
	stepsize = ceil(size(data,1)/samples);
	uniform = data( 1:stepsize:size(data,1), : );
    end
end

% decompose copied over from http://nghiaho.com/uploads/code/rotation_matrix_demo.m
function [x,y,z] = decompose_rotation(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end

function displayPointClouds(PointCloud1, PointCloud2)
    figure(1);
    clf();
    hold on
    stepsize = ceil(size(PointCloud1,1)/5000);
    uniform = PointCloud1( 1:stepsize:size(PointCloud1,1), : );
    fscatter3(uniform(:,1),uniform(:,2),uniform(:,3), [0 5]);
    stepsize = ceil(size(PointCloud2,1)/5000);
    uniform = PointCloud2( 1:stepsize:size(PointCloud2,1), : );
    fscatter3(uniform(:,1),uniform(:,2),uniform(:,3), [5 10]);
    hold off
end

function out = sampleNormalSpace(in, sampleSize, normals)
    % were apparently not needed since the normals were supplied
    % pointCloud2mesh is created by Ajmal Saeed Mian
    % see specific .m file for license and credits
    if nargin == 2
        mesh = pointCloud2mesh(in);
        normals = mesh.vertexNormals;
    end
    % createIcosahedron is created by David Legland
    % see specific .m file for license and credits
    icosahedron = createIcosahedron(); % 20 sides as bins for the normals
    normals_icosahedron = surfaceNormals( icosahedron );
    similarity = normals * normals_icosahedron';
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

