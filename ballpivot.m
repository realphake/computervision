function [tri, data] = ballpivot( input )
    % POINTCLOUD2TRIANGLES 
    % delaunay doesn't work properly, trying ball pivoting algorithm
    % patented by IBM, does this matter?
    % http://www.research.ibm.com/vistechnology/pdf/bpa_tvcg.pdf
    % ro = radius of the ball that pivots
    data = input(1:1000, 1:3);
    normals = input(1:1000, 5:7);
    KD_treesearcher = KDTreeSearcher(data);
    %ro_vec = mean(sqrt((data(1:length(data)-1, 1:3) - data(2:length(data), 1:3)).^2));
    %ro = sqrt(sum(ro_vec.^2));
    ro = 0.004;
    toBeEvaluated = 1:length(data);
    % keep list of boundary edges
    boundary_edges = zeros(length(data), 2);
    boundary_counter = 0;
    % keep list of frozen edges
    % frozen_edges = [];
    used_indices = zeros(length(data));
    used_indices_counter = 0;
    tri = zeros(length(data), 3);
    tri_counter = 0;
    stillLeft = toBeEvaluated(toBeEvaluated~=0);
    done = 0;
    pause on
    while ~done
        % find seed triangle
        % pick random point
        found = 0;
        counter = 0;
        while ~found && ~done
            if counter == length(stillLeft)
                done = 1;
            else
                idx = stillLeft(mod(counter+ceil(length(stillLeft)/2),length(stillLeft))+1);
                I = rangesearch(KD_treesearcher, data(idx,:), ro);
                I = cell2mat(I);
                I(ismember(I, used_indices(1:used_indices_counter))) = [];
                % find the smallest triangle
                min_dist = ro*3;
                smallest_triangle = [idx, 0, 0];
                for i=I(1:length(I)-1)
                    for j=I(2:length(I))
                        if i ~= j && i ~= idx && j ~= idx
                            dist = sqrt(sum(sum([data(idx,:)-data(i,:);data(idx,:)-data(j,:);data(j,:)-data(i,:);].^2, 1), 2));
                            dist12 = sqrt(sum((data(idx,:)-data(i,:)).^2, 2));
                            dist23 = sqrt(sum((data(i,:)-data(j,:)).^2, 2));
                            dist31 = sqrt(sum((data(j,:)-data(idx,:)).^2, 2));
                            if dist < min_dist && dist12 < ro && dist23 < ro && dist31 < ro
                                min_dist = dist;
                                smallest_triangle = [idx, i , j];
                            end
                        end
                    end
                end
                if ~ismember(0, smallest_triangle)
                    found = 1;
                else
                    counter = counter+1;
                end
            end
            
        end
        if ~done
        % keep queue of active edges (i.e. edges that have not been pivoted around yet 

        front = smallest_triangle;
        active_edges = [front([1, 2, 3]); front([2, 3, 1]); front([3, 1, 2])];
        % create triangles with the random point where distance between two
        % neighbors is < 2ro
        tri(tri_counter+1,:) = smallest_triangle;
        tri_counter = tri_counter+1;
        % ball pivoting
        active_index = 1;
        % while we have active edges, select next edge
        while active_index <= length(active_edges)
            % find closest pivoted point from edge
            A = data(active_edges(active_index,1),:);
            B = data(active_edges(active_index,2),:);
            C = data(active_edges(active_index,3),:);
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
                disp('No valid points within reach, adding to boundary');
                boundary_edges(boundary_counter+1,:) = [active_edges(active_index,1), active_edges(active_index,2)];
                boundary_counter = boundary_counter+1;
            else
                disp(['Has ', num2str(length(evaluateThese)), ' left to be evaluated.']);
                min_point = 0;
                % compute the intersections of three spheres of radius ro
                % centered in A B and C
                % trilateration was implemented by Bruno Luong
                % see the function for license notes
                [solution0, solution1, ~] = trilateration(A', B', C', [ro, ro, ro]');
                solution0 = solution0';
                solution1 = solution1';
                % for both solutions we check if the normal is pointing in
                % the same direction as the vertex normals
                A_normal = normals(active_edges(active_index,1),:);
                B_normal = normals(active_edges(active_index,2),:);
                C_normal = normals(active_edges(active_index,3),:);
                centered_sol0 = solution0 - mean([A; B; evaluatePoints(i,:)]);
                surfaceNormal0 = centered_sol0 ./ sqrt(sum(centered_sol0.^2, 2));
                centered_sol1 = solution1 - mean([A; B; evaluatePoints(i,:)]);
                surfaceNormal1 = centered_sol1 ./ sqrt(sum(centered_sol1.^2, 2));
                average_face_normal = mean([A_normal;B_normal;C_normal]);
                if dot(surfaceNormal0, average_face_normal) >= 0
                    active_edge_center = solution0 - midpoint;
                elseif dot(surfaceNormal1, average_face_normal) >= 0
                    active_edge_center = solution1 - midpoint;
                end
                % To properly compute the directed angle between the found
                % centers and our initial ball center:
                % first find the projection of C onto the midpoint
                % edge
                norm_aec = active_edge_center ./ sqrt(sum(active_edge_center.^2, 2));
                length_AB = sqrt(sum((A-B).^2, 2));
                point_on_AB = project_point_to_line_segment((A-midpoint)*ro/(length_AB*2),(B-midpoint)*ro/(length_AB*2),C-midpoint);
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
                if ref_angle > pi
                    ref_angle = ref_angle-(2*pi);
                elseif ref_angle < -pi
                    ref_angle = ref_angle+(2*pi);
                end
                % intialize these variables
                centers = zeros(length(evaluateThese), 3);
                angles = zeros(length(evaluateThese), 1)+2*pi;
                evaluatePoints = data(evaluateThese, :);
                evaluateNormals = normals(evaluateThese, :);
                % start a parallelized for loop to find each center and
                % compute the angle between the initial ball center and the
                % found center
                parfor i=1:length(evaluateThese)
                    % trilateration was implemented by Bruno Luong
                    % see the function for license notes
                    [solution0, solution1, errorflag] = trilateration(A', B', evaluatePoints(i,:)', [ro, ro, ro]');
                    solution0 = solution0';
                    solution1 = solution1';
                    % errorflag 0 means OK, < 0 is not OK
                    if errorflag == 0
                        % check the normals again
                        centered_sol0 = solution0 - mean([A; B; evaluatePoints(i,:)])
                        surfaceNormal0 = centered_sol0 ./ sqrt(sum(centered_sol0.^2, 2));
                        centered_sol1 = solution1 - mean([A; B; evaluatePoints(i,:)])
                        surfaceNormal1 = centered_sol1 ./ sqrt(sum(centered_sol1.^2, 2));
                        average_face_normal = mean([A_normal;B_normal;evaluateNormals(i, :)]);
                        invalid = 0;
                        if dot(surfaceNormal0, average_face_normal) >= 0
                            centers(i, :) = solution0 - midpoint;
                        elseif dot(surfaceNormal1, average_face_normal) >= 0
                            centers(i, :) = solution1 - midpoint;
                        else
                            centers(i, :) = [NaN, NaN, NaN];
                            invalid = 1;
                        end
                        % do not accept the point if the normals do not
                        % coincide
                        if ~invalid
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
                            angles(i) = angle
                        end
                    else
                        centers(i, :) = [NaN, NaN, NaN];
                        angles(i) = NaN;
                    end
                end
                
                % minimum angle neighbor selection
                if ~isempty(centers(~isnan(centers)))
                    [~, ind] = min(angles);
                    min_point = evaluateThese(ind);
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
                         front = [front, min_point];
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
                         output2 = quickfind2d(active_edges(:,[1,2]), min_point, active_edges(active_index,2));
                         if isempty(output1)
                             active_edges = [active_edges; min_point, active_edges(active_index,2), active_edges(active_index,1) ];
                         elseif output1 > active_index
                             active_edges(output1, :) = [];
                         end
                         if isempty(output2)
                             active_edges = [active_edges; active_edges(active_index,1), min_point, active_edges(active_index,2)];
                         elseif output2 > active_index
                             active_edges(output2, :) = [];
                         end
                     end
                 else
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
        used_indices(used_indices_counter+1:used_indices_counter+length(front)) = front;
        used_indices_counter = used_indices_counter+length(front);
        end
    end
    pause off
    tri = tri(1:tri_counter,:);
    trimesh(tri(1:tri_counter,:), data(:, 1), data(:, 2), data(:, 3));
end
