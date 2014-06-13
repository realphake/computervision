% tested using matlab 2012b
% parameters
% pointcloud used
PC = readPcd('Norbert Heijne.pcd'); % or arjen swellengrebel.pcd
% ballradius used for BPA, 0.004 is default, if 0 or not given then the
% algorithm will use a heuristic to estimate the radius
ballradius = 0.004;
% the amount of neighbors used for normal smoothing
neighbors_for_smoothing = 16;
% sample density from 0 to 1
density = 0.5;
% outlier tolerance, points that have a distance from the mean of the data
% larger than the mean_distance + tolerance * std_distance are cut off
tolerance = 2;
%first 1000 points for illustration
%raw_datapoints = PC(1:1000, 1:3);
%raw_normals = PC(1:1000, 5:7);
% the head also makes for a good illustration
head = PC(PC(:, 2) < -0.5,:);
raw_datapoints = head(:, 1:3); % change to PC for entire mesh: PC(1:3, 5:7);
raw_normals = head(:, 5:7); % change to PC for entire mesh: PC(1:3, 5:7);
% sample uniformly with density as parameter
stepsize = ceil(size(raw_datapoints,1)/ceil(size(raw_datapoints,1)*density));
datapoints = raw_datapoints( 1:stepsize:size(raw_datapoints,1), : );
normals = raw_normals(1:stepsize:size(raw_normals,1), :);
% noisy data comes with outliers, we filter them out
outliers = findPcdOutliers(datapoints, tolerance);
datapoints(outliers,:) = [];
normals(outliers,:) = [];
% smoothing the normals, kinect normals aren't all that reliable
smoothedNormals = smoothNormals( datapoints, normals, neighbors_for_smoothing );
% performing the BPA
[tri, boundaries,b_empty, b_alreadyused, b_invalid] = ballpivot(datapoints, smoothedNormals, ballradius);
% stitching together the boundaries to make it watertight
boundaryfaces = stitchBoundaries(boundaries, datapoints);
boundaryfaces = unique(boundaryfaces, 'rows');

% the resulting mesh
end_result = [tri;boundaryfaces];

% visualize
hold on
% the mesh from BPA
trimesh(tri, datapoints(:, 1), datapoints(:, 2), datapoints(:, 3), 'FaceColor',[0.8, 0.8, 0.8]);
% the mesh from the stitching
trimesh(boundaryfaces, datapoints(:, 1), datapoints(:, 2), datapoints(:, 3), 'FaceColor',[0.4, 0.4, 0.4]);
% visualizing the normals
quiver3(datapoints(:,1),datapoints(:,2),datapoints(:,3),normals(:,1),normals(:,2),normals(:,3), 0.5);
% visualizing the boundaries
A = datapoints(boundaries(:, 1), :);
B = datapoints(boundaries(:, 2), :);
x = [A(:,1), B(:,1)]';
y = [A(:,2), B(:,2)]';
z = [A(:,3), B(:,3)]';
plot3(x, y, z, 'k', 'LineWidth', 2);
% Types of boundaries for debugging purposes
% if ~isempty(b_empty)
%     A = datapoints(b_empty(:, 1), :);
%     B = datapoints(b_empty(:, 2), :);
%     x = [A(:,1), B(:,1)]';
%     y = [A(:,2), B(:,2)]';
%     z = [A(:,3), B(:,3)]';
%     plot3(x, y, z, 'r', 'LineWidth', 2);
% end
% if ~isempty(b_alreadyused)
%     A = datapoints(b_alreadyused(:, 1), :);
%     B = datapoints(b_alreadyused(:, 2), :);
%     x = [A(:,1), B(:,1)]';
%     y = [A(:,2), B(:,2)]';
%     z = [A(:,3), B(:,3)]';
%     plot3(x, y, z, 'g', 'LineWidth', 2);
% end
% if ~isempty(b_invalid)
%     A = datapoints(b_invalid(:, 1), :);
%     B = datapoints(b_invalid(:, 2), :);
%     x = [A(:,1), B(:,1)]';
%     y = [A(:,2), B(:,2)]';
%     z = [A(:,3), B(:,3)]';
%     plot3(x, y, z, 'b', 'LineWidth', 2);
% end
% show outliers
%scatter3(PC(outliers, 1),PC(outliers, 3),PC(outliers, 3),4, [0,0,0]);

hold off