function ret = distances(cloud1, cloud2)
    ret = zeros(size(cloud1,1),size(cloud1,1));
    for i = 1:size(cloud1,1)
        for j = 1:size(cloud2,1)
            p1 = cloud1(i,1:3);
            p2 = cloud2(j,1:3);
            d = pdist([p1;p2],'euclidean');
            ret(i,j) = d;
        end;
    end;
end
