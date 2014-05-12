function [noNoise, normals] = removeNoise(noise, cutoffPoint, normals )
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
