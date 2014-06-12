function [output, R_total, T_total, iterations] = ICP( base, target, sampleSize, sampleTech, maxIterations, targetNormals )

    % In this implementation we fit the target cloud onto the base cloud,
    % as implemented in the supplied paper: 
    % Implementation of a 3D ICP-based Scan Matcher
    
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
    if ( strcmpi( sampleTech, 'normals' ) ) && nargin == 6
        baseCloud = subsampling(base, sampleSize, sampleTech, targetNormals);
    else
        baseCloud = subsampling(base, sampleSize, sampleTech);
    end
    searchAssist = baseCloud;
    base_length = length(baseCloud);
    % target is not sampled: the search is done in all its points
    KD_treesearcher = KDTreeSearcher(target);
    error = 0;
    iterations = 0;
    epsilon_err = 0.000000000000001;
    good_enough = false;
    while ~ good_enough
        % update target cloud
        target = (target * R') + repmat(T,length(target),1);
        % only needs resampling if random, the other methods are deterministic
        if strcmpi( sampleTech, 'random' ) 
            baseCloud = subsampling(base, sampleSize, sampleTech);
        end
        % search for the nearest neighbor for every point in the baseCloud
        searchAssist = (searchAssist - repmat(T,base_length,1)) * R ;
        IDX = knnsearch(KD_treesearcher,searchAssist);
        %IDX = knnsearch(target,baseCloud);
        targetCloud = target(IDX,:);
        % compute mean for centering
        baseCloudCenter = mean(baseCloud,1);
        targetCloudCenter = mean(targetCloud,1);
        % perform SVD on the matrix as explained in the paper
        A = (baseCloud - repmat(baseCloudCenter, base_length,1))' * (targetCloud - repmat(targetCloudCenter, base_length,1));
        [U,S,V] = svd(A);
        % get rotation from the SVD
        R = U*V';
        R_total = R_total*R';
        % after rotating the target cloud we calculate the distance between
        % matched clouds
        T = baseCloudCenter - (targetCloudCenter * R');
        T_total = T_total + T;
        % we update the target cloud to compute the error
        targetCloud = (targetCloud * R') + repmat(T,base_length,1);
        old_error = error;
        % error is approximated as the root mean squared distance of the
        % two clouds
        error = pdist([rms(baseCloud-targetCloud,1); 0 0 0], 'euclidean');
        iterations = iterations + 1;
		disp(['iteration: ', num2str(iterations)]);
        % the to be evaluated distance is normalized such that we can use
        % our epsilon
        diff_err_normalized = abs(( error - old_error ) / old_error);
        if maxIterations == 0
            % if no max iterations is given then continue untill error is
            % small enough
            good_enough = diff_err_normalized < epsilon_err;
        else
            good_enough = diff_err_normalized < epsilon_err || iterations == maxIterations;
        end
    end
    output = (target * R') + repmat(T,length(target),1);
end

function out = subsampling(in, sampleSize, technique, normals)
    % main function for our subsampling methods
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
    % take the specified number of samples from the data with a equal
    % stepsize
    if length(data) < samples
        uniform = data;
    else
	stepsize = ceil(size(data,1)/samples);
	uniform = data( 1:stepsize:size(data,1), : );
    end
end

% decompose copied over from http://nghiaho.com/uploads/code/rotation_matrix_demo.m
% used to make the rotation matrix interpretable
function [x,y,z] = decompose_rotation(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end

function out = sampleNormalSpace(in, sampleSize, normals)
    % createIcosahedron is created by David Legland
    % see specific .m file for license and credits
    icosahedron = createIcosahedron(); % 20 sides as bins for the normals
    normals_icosahedron = surfaceNormals( icosahedron );
    % by taking the dot product we can measure the similarity between
    % normals
    similarity = normals * normals_icosahedron';
    [values, binning] = max(similarity, [], 2);
    out = [];
    for i = 1:20
        out = [ out ; uniformSubsampling(in(binning == i, :), sampleSize)];
    end
end

% copied from pointCloud2mesh.m and adjusted to only get the surface
% normals
% Author : Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
%           Computer Science. Univ of Western Australia
% http://www.csse.uwa.edu.au/~ajmal/code/pointCloud2mesh.m
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

