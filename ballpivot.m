function [tri, boundary_edges, b_empty, b_alreadyused, b_invalid] = ballpivot( data, normals, ro)
    % POINTCLOUD2TRIANGLES 
    % delaunay doesn't work properly, trying ball pivoting algorithm
    % patented by IBM, does this matter?
    % http://www.research.ibm.com/vistechnology/pdf/bpa_tvcg.pdf
    % ro = radius of the ball that pivots
    KD_treesearcher = KDTreeSearcher(data);
    if nargin < 3 || ro == 0;
        % heuristic, take the average point distance between consecutive
        % points
        ro_vec = mean(sqrt((data(1:length(data)-1, 1:3) - data(2:length(data), 1:3)).^2));
        ro = sqrt(sum(ro_vec.^2));
    end
    % variables to be initialized
    toBeEvaluated = 1:length(data);
    boundary_edges = zeros(length(data), 2);
    boundary_counter = 0;
    front = [];
    used_indices = zeros(length(data), 1);
    used_indices_counter = 0;
    tri = zeros(length(data)*3, 3);
    tri_counter = 0;
    stillLeft = toBeEvaluated(toBeEvaluated~=0);
    done = 0;
    b_empty = [];
    b_alreadyused = [];
    b_invalid = [];
    pause on
    while ~done
        % find seed triangle
        % pick a point somewhere in the middle
        found = 0;
        counter = 0;
        while ~found && ~done
            if counter == length(stillLeft)
                % when we've checked everything and haven't found a
                % triangle then we're done
                done = 1;
            else
                % as long as there are points not yet tried we keep finding
                % seed triangles
                idx = stillLeft(mod(counter+ceil(length(stillLeft)/2),length(stillLeft))+1);
                toBeEvaluated(idx) = 0;
                % find all points within a 2*ro distance
                I = rangesearch(KD_treesearcher, data(idx,:), 2*ro);
                I = cell2mat(I);
                neighbor_data_all = data(I,:);
                % we don't want to make triangles with points that we've
                % already used
                dontUseThese = ismember(I, [idx;used_indices(1:used_indices_counter); front]);
                neighbor_data = neighbor_data_all(dontUseThese == 0,:);
                I(dontUseThese) = [];
                if ~isempty(I)
                    distance = sqrt(sum((neighbor_data - repmat(data(idx,:), size(neighbor_data, 1), 1)).^2, 2));
                    % sort in order of distance from small to large
                    [~, ind_sort]=sort(distance);
                    I = I(ind_sort);
                    % find a good triangle
                    for i=I(1:length(I)-1)
                        for j=I(2:length(I))
                            if i ~= j && i ~= idx && j ~= idx
                                A = data(idx,:);
                                B = data(i,:);
                                C = data(j,:);
                                % has to intersect with 3 spheres with radius
                                % ro
                                [solution0, solution1, err] = sphereIntersection( A, B, C, ro, ro, ro );
                                if err == 0
                                    A_normal = normals(idx,:);
                                    B_normal = normals(i,:);
                                    C_normal = normals(j,:);
                                    % normal has to coincide with the face normal
                                    % of the 3 points
                                    selection = checkNormals(solution0, solution1, [A;B;C], [A_normal; B_normal; C_normal]);
                                    if selection ~= -1
                                        % consider all neighbors when
                                        % checking if the ball contains no
                                        % other points
                                        if selection == 0
                                            distance = sqrt(sum((neighbor_data_all - repmat(solution0, size(neighbor_data_all, 1), 1)).^2, 2));
                                        else
                                            distance = sqrt(sum((neighbor_data_all - repmat(solution1, size(neighbor_data_all, 1), 1)).^2, 2));
                                        end
                                        if sum(distance <= ro) == 3
                                            % when we find a good triangle
                                            % we're going to continue
                                            good_triangle = [idx, i , j];
                                            found = 1;
                                            break;
                                        end
                                    end
                                end
                            end
                        end
                        if found; break; end;
                    end
                end
                if ~found
                    counter = counter+1;
                end
            end
            
        end
        if ~done
        % keep queue of active edges (i.e. edges that have not been pivoted around yet 
        front = [front;good_triangle'];
        % active edges are built like this (edge from (1) to (2) that are 
        % both connected to (3) which is used to calculate the center of the sphere)
        active_edges = [good_triangle([1, 2, 3]); good_triangle([2, 3, 1]); good_triangle([3, 1, 2])];
        % create triangles with the random point where distance between two
        % neighbors is < 2ro
        tri(tri_counter+1,:) = good_triangle;
        tri_counter = tri_counter+1;
        % ball pivoting
        active_index = 1;
        % while we have active edges, select next edge
        while active_index <= length(active_edges)
            % find closest pivoted point from edge
            A = data(active_edges(active_index,1),:);
            B = data(active_edges(active_index,2),:);
            C = data(active_edges(active_index,3),:);
            A_normal = normals(active_edges(active_index,1),:);
            B_normal = normals(active_edges(active_index,2),:);
            C_normal = normals(active_edges(active_index,3),:);
            midpoint = (A+B)./2;
            % search for neighbors from the midpoint in a 2*ro distance
            neighbors = rangesearch(KD_treesearcher, midpoint, 2* ro);
            neighbors = cell2mat(neighbors);
            disp(['Has ', num2str(length(neighbors)), ' neighbors.']);
            % we evaluate all the neighbors except the 3 initial points
            evaluateThese = neighbors;
            evaluateThese(ismember(evaluateThese, active_edges(active_index,:))) = [];
            if isempty(evaluateThese)
                % if the initial triangle is are the only points within
                % reach then add to boundary
                b_empty = [b_empty;active_edges(active_index,1), active_edges(active_index,2)];
                disp('No valid points within reach, adding to boundary');
                boundary_edges(boundary_counter+1,:) = [active_edges(active_index,1), active_edges(active_index,2)];
                boundary_counter = boundary_counter+1;
            else
                disp(['Has ', num2str(length(evaluateThese)), ' left to be evaluated.']);
                % compute the intersections of three spheres of radius ro
                % centered in A B and C
                [solution0, solution1, ~] = sphereIntersection( A, B, C, ro, ro, ro );
                % for both solutions we check which of the normals is pointing in
                % the same direction as the vertex normals
                selection = checkNormals(solution0, solution1, [A;B;C], [A_normal; B_normal; C_normal]);
                if selection == 0
                    active_edge_center = solution0 - midpoint;
                elseif selection == 1
                    active_edge_center = solution1 - midpoint;
                else
                    disp('This shouldnt happen, selection = -1');
                end
                % To properly compute the directed angle between the found
                % centers and our initial ball center:
                % first find the projection of C onto the midpoint
                % edge
                norm_aec = active_edge_center ./ sqrt(sum(active_edge_center.^2, 2));
                length_AB = sqrt(sum((A-B).^2, 2));
                point_on_AB = project_point_to_line_segment((A-midpoint)*ro/(length_AB/2),(B-midpoint)*ro/(length_AB/2),C-midpoint);
                % get normalized vector from that and compute cross product
                % for a vector that can be rotated to point in the Z
                % direction
                centered_C = C-midpoint-point_on_AB;
                norm_C = centered_C ./ sqrt(sum(centered_C.^2, 2));
                cross_aec_C = cross(norm_aec, norm_C);
                % rotate cross product of norm_C and norm_aec to [0, 0, 1]
                R_vec = vrrotvec(cross_aec_C ./ sqrt(sum(cross_aec_C.^2, 2)),[0, 0, 1]);
                R = vrrotvec2mat(R_vec);
                % keep rotation for every center found
                % then determine rotation via atan2 using x and y
                % of the found centers
                rotated_C = norm_C * R';
                rotated_aec = norm_aec * R';
                % from vector 1 to vector 2
                %angle = atan2(vector2.y, vector2.x) - atan2(vector1.y, vector1.x);
                ref_angle = atan2(rotated_C(2), rotated_C(1)) - atan2(rotated_aec(2), rotated_aec(1));
                % we now have a reference angle that we can use to
                % determine which way we should rotate
                % we SHOULD rotate away from this reference point
                % reference angle should also be the smallest directed
                % angle from C (the back point of the triangle) to the
                % active edge center (center of the sphere that lies on top of the triangle)
                if ref_angle > pi
                    ref_angle = ref_angle-(2*pi);
                elseif ref_angle < -pi
                    ref_angle = ref_angle+(2*pi);
                end
                % intialize these variables
                centers = zeros(length(evaluateThese), 3)*NaN;
                angles = zeros(length(evaluateThese), 1)*NaN;
                evaluatePoints = data(evaluateThese, :);
                evaluateNormals = normals(evaluateThese, :);
                valid = ones(length(evaluateThese), 1);
                neighbor_data = data(neighbors,:);
                % start a parallelized for loop to find each center and
                % compute the angle between the initial ball center and the
                % found center
                parfor i=1:length(evaluateThese)
                    [solution0, solution1, err] = sphereIntersection( A, B, evaluatePoints(i,:), ro, ro, ro );
                    % errorflag 0 means OK, 1 is not OK
                    if err == 1
                        valid(i) = 0;
                    elseif err == 0
                        % check the normals again
                        %
                        selection = checkNormals(solution0, solution1, [A;B;evaluatePoints(i,:)], [A_normal; B_normal; evaluateNormals(i, :)]);
                        % if there are other points inside the sphere it is
                        % invalid
                        if selection == 0
                            distance = sqrt(sum((neighbor_data - repmat(solution0, size(neighbor_data, 1), 1)).^2, 2));
                            if sum(distance <= ro) > 3
                                valid(i) = 0;
                            end
                        elseif selection == 1
                            distance = sqrt(sum((neighbor_data - repmat(solution1, size(neighbor_data, 1), 1)).^2, 2));
                            if sum(distance <= ro) > 3
                                valid(i) = 0;
                            end
                        else
                            valid(i) = 0;
                        end
                        if selection == 0
                            centers(i, :) = solution0 - midpoint;
                        elseif selection == 1
                            centers(i, :) = solution1 - midpoint;
                        else
                            % we still need an angle
                            centers(i, :) = solution0 - midpoint;
                        end
                        % we rotate the found center with R to cancel
                        % out the Z direction
                        normalized_center = centers(i, :) ./ sqrt(sum(centers(i, :).^2, 2));
                        normalized_center = normalized_center * R';
                        % compute the directed angle between the
                        % initial center and the found center
                        angle = atan2(normalized_center(2), normalized_center(1)) - atan2(rotated_aec(2), rotated_aec(1));
                        % if our reference angle is <0 then the angle
                        % we want is >0 and vice versa
                        if ref_angle < 0 && angle < 0
                            angle = angle+(2*pi);
                        elseif ref_angle > 0 && angle > 0
                            angle = angle-(2*pi);
                        end
                        angles(i) = angle;
                        % for the case that the other angle might be
                        % smaller
                        if selection == -1
                            other_center = solution1 - midpoint;
                            normalized_center = other_center ./ sqrt(sum(other_center.^2, 2));
                            normalized_center = normalized_center * R';
                            other_angle = atan2(normalized_center(2), normalized_center(1)) - atan2(rotated_aec(2), rotated_aec(1));
                            if ref_angle < 0 && other_angle < 0
                                other_angle = other_angle+(2*pi);
                            elseif ref_angle > 0 && other_angle > 0
                                other_angle = other_angle-(2*pi);
                            end
                            if abs(other_angle) < abs(angle)
                                angles(i) = other_angle;
                            end
                        end
                    end
                end
                
                min_point = 0;
                invalid_point = 0; 
                % minimum angle neighbor selection
                if ~isempty(centers(~isnan(centers)))
                    [~, ind] = min(abs(angles));
                    if valid(ind)
                        min_point = evaluateThese(ind);
                    else
                        invalid_point = 1; 
                    end
                end
                
                % if the first encountered point is not usable then we will
                % add this edge to the boundary instead
                foundAlreadyUsed = ismember(min_point, used_indices(1:used_indices_counter));
                % foundAlreadyUsed can be empty in the case that
                % used_indices has no entries yet
                if isempty(foundAlreadyUsed)
                    foundAlreadyUsed = 0;
                end
                
                 if min_point && ~foundAlreadyUsed
                     % we now create a triangle with the new point;
                     disp(['Found minimal point: ', num2str(min_point)]);
                     tri(tri_counter+1,:) = [active_edges(active_index,1), min_point, active_edges(active_index,2)];
                     tri_counter = tri_counter+1;
                     %trimesh(tri(1:tri_counter,:), data(:, 1), data(:, 2), data(:, 3));
                     %pause(3);
                     if ~ismember(min_point, front)
                         % if the point is new we add both new edges to the
                         % list
                         disp('Adding point to front and active edges');
                         % is a new point so we don't have to glue
                         active_edges = [active_edges; active_edges(active_index,1), min_point, active_edges(active_index,2); min_point, active_edges(active_index,2), active_edges(active_index,1)];
                         front = [front; min_point];
                     else
                         disp('Glueing');
                         % we'll have to glue
                         % we remove any edge still in the queue that now
                         % overlaps with the newly formed triangle
                         % remove the point that is behind the front and
                         % add it to the used_indices
                         if ismember(active_edges(active_index,1), front)
                             used_indices(used_indices_counter+1) = active_edges(active_index,1);
                             used_indices_counter = used_indices_counter+1;
                             toBeEvaluated(idx) = 0;
                             front(front == active_edges(active_index,1)) = [];
                         end
                         if ismember(active_edges(active_index,2), front)
                             used_indices(used_indices_counter+1) = active_edges(active_index,2);
                             used_indices_counter = used_indices_counter+1;
                             toBeEvaluated(idx) = 0;
                             front(front == active_edges(active_index,2)) = [];
                         end
                         output1 = quickfind2d(active_edges(:,[1,2]), min_point, active_edges(active_index,1));
                         if isempty(output1)
                             active_edges = [active_edges; min_point, active_edges(active_index,1), active_edges(active_index,2) ];
                         elseif output1 > active_index
                             active_edges(output1, :) = [];
                         end
                         output2 = quickfind2d(active_edges(:,[1,2]), min_point, active_edges(active_index,2));
                         if isempty(output2)
                             active_edges = [active_edges; active_edges(active_index,2), min_point, active_edges(active_index,1)];
                         elseif output2 > active_index
                             active_edges(output2, :) = [];
                         end
                         output3 = quickfind2d(boundary_edges(:,[1,2]), min_point, active_edges(active_index,1));
                         if ~isempty(output3)
                             boundary_edges(output3, :) = [];
                         end
                         output4 = quickfind2d(boundary_edges(:,[1,2]), min_point, active_edges(active_index,2));
                         if ~isempty(output4)
                             boundary_edges(output4, :) = [];
                         end
                     end
                 else
                     if foundAlreadyUsed
                         b_alreadyused = [b_alreadyused;active_edges(active_index,1), active_edges(active_index,2)];
                     elseif invalid_point
                         b_invalid = [b_invalid;active_edges(active_index,1), active_edges(active_index,2)];
                     end
                     disp('Could not find point within reach, adding to boundary');
                     boundary_edges(boundary_counter+1,:) = [active_edges(active_index,1), active_edges(active_index,2)];
                     boundary_counter = boundary_counter+1;
                 end
            end
            active_index = active_index +1;
            disp('active index now:');
            disp(active_index);
        end
        % we've finished with the initial seed and constructed the boundary,
        % all points on the front should now correspond with the boundary, all
        % glued points have already been added to used_indices
        toBeEvaluated(front) = 0;
        stillLeft = toBeEvaluated(toBeEvaluated~=0);
        %used_indices(used_indices_counter+1:used_indices_counter+length(front)) = front;
        %used_indices_counter = used_indices_counter+length(front);
        end
    end
    pause off
    tri = tri(1:tri_counter,:);
    boundary_edges = boundary_edges(1:boundary_counter, :);
    [row, col] = find(boundary_edges==0);
    row = unique(row);
    boundary_edges(row,:) = [];
    [row, col] = find(b_empty==0);
    row = unique(row);
    b_empty(row,:) = [];
    [row, col] = find(b_invalid==0);
    row = unique(row);
    b_invalid(row,:) = [];
    [row, col] = find(b_alreadyused==0);
    row = unique(row);
    b_alreadyused(row,:) = [];
    %tri = [tri;stitchBoundaries(boundary_edges)];
    tri = unique(tri, 'rows');
end

function selection = checkNormals(solution0, solution1, points, normals)
    % normal has to coincide with the face normal
    % of the 3 points
    centered_sol0 = solution0 - mean(points);
    surfaceNormal0 = centered_sol0 ./ sqrt(sum(centered_sol0.^2, 2));
    centered_sol1 = solution1 - mean(points);
    surfaceNormal1 = centered_sol1 ./ sqrt(sum(centered_sol1.^2, 2));
    if size(normals, 1) == 1
        average_face_dir = normals;
    else
        average_face_dir = mean(normals, 1);
    end
    average_face_normal = average_face_dir ./ sqrt(sum(average_face_dir.^2, 2));
    selection = -1;
    if dot(surfaceNormal0, average_face_normal) >= 0
        selection = 0;
    elseif dot(surfaceNormal1, average_face_normal) >= 0
        selection = 1;
    end
end
