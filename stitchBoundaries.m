function triangles = stitchBoundaries(boundaryList, data)
    % if the ballpivot works correctly it should find a surface all the mesh
    % and the boundaries should be the holes in the mesh
    triangles = zeros(length(boundaryList), 3);
    counter = 0;
    while ~isempty(boundaryList)
        % find boundaries that are stringed together
        found = 0;
        stitchThis = boundaryList(1,:);
        boundaryList(1,:) = [];
        [row, col] = find(boundaryList == stitchThis(1));
        if ~isempty(row)
            stitchThisConnection = 1;
            stitchFrom = 2;
            found = 1;
        else
            [row, col] = find(boundaryList == stitchThis(2));
            if ~isempty(row)
                stitchThisConnection = 2;
                stitchFrom = 1;
                found = 1;
            end
        end
        diff = data(boundaryList(row, col)) - data(stitchThis(stitchThisConnection));
        [~, min_row] = min(sqrt(sum(diff .^2, 2)));
        if found
            if col(1) == 1; stitchTo = 2; end;
            if col(1) == 2; stitchTo = 1; end;
            toThis = boundaryList(row(1), :);
            if stitchThis(stitchFrom) ~= toThis(stitchTo)
                triangles(counter+1, :) = [stitchThis, toThis(stitchTo)];
                counter = counter +1;
                boundaryList(row(min_row), :)= [stitchThis(stitchFrom), toThis(stitchTo)];
            end
        end
        
    end
    triangles = triangles(1:counter, :);
end


