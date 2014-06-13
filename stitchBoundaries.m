function triangles = stitchBoundaries(boundaryList, data)
    % if the ballpivot works correctly it should find a surface all the mesh
    % and the boundaries should be the holes in the mesh and the outside
    % perimeter
    triangles = zeros(length(boundaryList), 3);
    counter = 0;
    while ~isempty(boundaryList)
        % find boundaries that are stringed together
        % we take the first boundary
        stitchThis = boundaryList(1,:);
        boundaryList(1,:) = [];
        startNewStitch = 1;
        % as long as we are still stitching up a single hole
        while startNewStitch || found
            found = 0;
            % we find the boundaries that link to the stitchThis boundary
            [row1, col1] = find(boundaryList == stitchThis(1));
            [row2, col2] = find(boundaryList == stitchThis(2));
            min_dist = inf;
            min_row = 0;
            % determine the would be created edge with the smallest length
            if ~isempty(row1)
                for i=1:length(row1)
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
            % by only using row1 we alternate our stitches
            if (startNewStitch || ~found) && ~isempty(row2)
                for i=1:length(row2)
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
                % create a triangle with the found boundary edge
                triangles(counter+1, :) = [stitchThis, stitchTo];
                counter = counter +1;
                stitchThis = [stitchFrom, stitchTo];
                boundaryList(min_row, :)= [];
            end
            startNewStitch = 0;
        end
    end
    triangles = triangles(1:counter, :);
    triangles = unique(triangles, 'rows');
end


