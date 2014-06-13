function [noNoise, normals] = removeNoise(noise, cutoffPoint, normals )
    % noise: Noisy input data
    % cutoffPoint: Determines the value at which Z value we remove points
    % normals: if normals are supplied then the points which do not have
    % valid normals will also be removed
    noNoiseIndices = noise(:,3) < cutoffPoint;
    noNoise = noise(noNoiseIndices,:);
    noNoise = noNoise(:,1:3);
    if nargin == 3
        normals = normals(noNoiseIndices,:);
        usable_normals = ~any(isnan(normals), 2);
        noNoise = noNoise(usable_normals, :);
        normals = normals(usable_normals, :);
        normals = normals(:,1:3);
    else
        normals = 'not_used';
    end
end
