function triangles = stichBoundaries(boundaryList)
    % if the ballpivot works correctly it should find a surface all the mesh
    % and the boundaries should be the holes in the mesh
    triangles = zeros(length(boundaryList), 3);
    counter = 0;
    while ~isempty(boundaryList)
        % find boundaries that are stringed together
        found = 0;
        stichThis = boundaryList(1,:);
        boundaryList(1,:) = [];
        [row, col] = find(boundaryList == stichThis(1));
        if ~isempty(row)
            stichFrom = 2;
            found = 1;
        else
            [row, col] = find(boundaryList == stichThis(2));
            if ~isempty(row)
                stichFrom = 1;
                found = 1;
            end
        end
        
        if found
            if col == 1; stichTo = 2; end;
            if col == 2; stichTo = 1; end;
            toThis = boundaryList(row(1), :);
            if stichThis(stichFrom) ~= toThis(stichTo)
                triangles(counter+1, :) = [stichThis, toThis(stichTo)];
                counter = counter +1;
                boundaryList(row(1), :)= [stichThis(stichFrom), toThis(stichTo)];
            end
        end
        
    end
    triangles = triangles(1:counter, :);
end


