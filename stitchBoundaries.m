function triangles = stitchBoundaries(boundaryList, data)
    % if the ballpivot works correctly it should find a surface all the mesh
    % and the boundaries should be the holes in the mesh and the outside
    % perimeter
    triangles = zeros(length(boundaryList), 3);
    counter = 0;
    while ~isempty(boundaryList)
        % find boundaries that are stringed together
        
        stitchThis = boundaryList(1,:);
        boundaryList(1,:) = [];
        startNewStitch = 1;
        while startNewStitch || found
            startNewStitch = 0;
            found = 0;
            [row1, col1] = find(boundaryList == stitchThis(1));
            [row2, col2] = find(boundaryList == stitchThis(2));
            min_dist = inf;
            min_row = 0;
            if ~isempty(row1)
                for i=length(row1)
                    dist = sqrt(sum((data(boundaryList(row1(i), mod(col1(i), 2)+1), :) - data(stitchThis(2))).^2, 2));
                    if dist < min_dist && boundaryList(row1(i), mod(col1(i), 2)+1) ~= stitchThis(2)
                        min_dist = dist;
                        stitchFrom = stitchThis(2);
                        stitchTo = boundaryList(row1(i), mod(col1(i), 2)+1);
                        found = 1;
                        min_row = row1(i);
                    end
                end
            end
            if ~isempty(row2)
                for i=length(row2)
                    dist = sqrt(sum((data(boundaryList(row2(i), mod(col2(i), 2)+1), :) - data(stitchThis(1))).^2, 2));
                    if dist < min_dist && boundaryList(row2(i), mod(col2(i), 2)+1) ~= stitchThis(1)
                        min_dist = dist;
                        stitchFrom = stitchThis(1);
                        stitchTo = boundaryList(row2(i), mod(col2(i), 2)+1);
                        found = 1;
                        min_row = row2(i);
                    end
                end
            end
            if found
                triangles(counter+1, :) = [stitchThis, stitchTo];
                counter = counter +1;
                stitchThis = [stitchFrom, stitchTo];
                boundaryList(min_row, :)= [];
            end
        end
    end
    triangles = triangles(1:counter, :);
    triangles = unique(triangles, 'rows');
end


