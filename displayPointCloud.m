function displayPointCloud(out_pointcloud)
    stepsize = ceil(size(out_pointcloud,1)/20000);
    uniform = out_pointcloud( 1:stepsize:size(out_pointcloud,1), : );
    % fscatter was created by Felix Morsdorf for which we take no credit
    fscatter3(uniform(:,1),uniform(:,2),uniform(:,3), [0 10]);
end
