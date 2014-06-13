function smoothedNormals  = smoothNormals( datapoints, normals, neighbors )
    %SMOOTHNORMALS Summary of this function goes here
    %   Detailed explanation goes here
    KD_treesearcher = KDTreeSearcher(datapoints);
    NN = knnsearch(KD_treesearcher, datapoints, 'K', neighbors);
    NN_normals = zeros(length(datapoints),3,neighbors);
    parfor i=1:neighbors
        NN_normals(:,:,i) = normals(NN(:, i), :);
    end
    smoothedNormals = mean(NN_normals, 3);
    smoothedNormals = smoothedNormals ./ repmat(sqrt(sum(smoothedNormals.^2,2)),1,3);
end

