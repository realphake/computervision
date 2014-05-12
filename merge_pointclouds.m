function [ out_pointcloud, R, T ] = merge_pointclouds(useTheseNumbers, sampleSize, sampleTech, type) 
    % useTheseNumbers = a list of image numbers to use example: [0:1:99]
    % sampleSize = total amount of points sampled, except for normal space,
    % which indicates the samples per bin
    % sampleTech = 'none', 'uniform', 'random', 'normals'
    % type = 'merge' which merges each consecutive frame into a new 
    % pointcloud and uses it as base or 'estimate' which finds the merge 
    % for consecutive frames and uses the individual transformation details
    % to create the final pointcloud.
    if strcmpi( type, 'merge')
        [out_pointcloud, R, T] = merge_then_estimate(useTheseNumbers, sampleSize, sampleTech);
    elseif strcmpi( type, 'estimate')
        [out_pointcloud, R, T] = estimate_then_merge(useTheseNumbers, sampleSize, sampleTech);
    else
        disp('Type unknown: exiting');
    end
    displayPointCloud(out_pointcloud);
end
