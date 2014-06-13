% tested using matlab 2012b
PC = readPcd('Norbert Heijne.pcd');

% sample uniformly

% the head also makes for a good illustration
head = PC(PC(:, 2) < -0.5,:);
datapoints = head(:, 1:3);
normals = head(:, 5:7);
%first 1000 points for illustration
%datapoints = PC(1:1000, 1:3);
%normals = PC(1:1000, 5:7);
ballradius = 0.004;
outliers = findPcdOutliers(datapoints);
datapoints(outliers,:) = [];
normals(outliers,:) = [];
smoothedNormals = smoothNormals( datapoints, normals, 16 );
[tri, boundaries,b_empty, b_alreadyused, b_invalid] = ballpivot(datapoints, smoothedNormals, ballradius);
boundaryfaces = stitchBoundaries(boundaries, datapoints);
boundaryfaces = unique(boundaryfaces, 'rows');
% also want to show normals and points which have been missed

hold on
trimesh(tri, datapoints(:, 1), datapoints(:, 2), datapoints(:, 3), 'FaceColor',[0.8, 0.8, 0.8]);
trimesh(boundaryfaces, datapoints(:, 1), datapoints(:, 2), datapoints(:, 3), 'FaceColor',[0.4, 0.4, 0.4]);
quiver3(datapoints(:,1),datapoints(:,2),datapoints(:,3),normals(:,1),normals(:,2),normals(:,3), 0.5);
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