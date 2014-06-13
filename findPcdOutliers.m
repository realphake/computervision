function [ indices ] = findPcdOutliers( datapoints, tolerance )
% expects datapoints to be of size N by 3
% first centers the datapoints
centered_datapoints = datapoints - repmat(mean(datapoints, 1), length(datapoints), 1);
% measure the distance to each point
distances = sqrt(sum(centered_datapoints.^2,2));
mean_distance = mean(distances);
std_distance = std(distances);
indices = (1:length(datapoints))';
% designate any point that is beyond the mean+std in distance as outlier 
indices = indices(distances > mean_distance + tolerance*std_distance);
end

