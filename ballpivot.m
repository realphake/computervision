function [tri, data] = ballpivot( data )
    % POINTCLOUD2TRIANGLES 
    % delaunay doesn't work properly, trying ball pivoting algorithm
    % patented by IBM, does this matter?
    % http://www.research.ibm.com/vistechnology/pdf/bpa_tvcg.pdf
    % ro = radius of the ball that pivots
    data = data(1:1000, 1:3);
    data = data - repmat(mean(data), length(data), 1); % centered
    KD_treesearcher = KDTreeSearcher(data);
    ro_vec = mean(sqrt((data(1:length(data)-1, 1:3) - data(2:length(data), 1:3)).^2));
    ro = sqrt(sum(ro_vec.^2));
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
    while ~done
        % find seed triangle
        % pick random point
        found = 0;
        counter = 1;
        while ~found
            if counter == length(stillLeft)+1
                done = 1;
            else
                idx = stillLeft(counter);
                % find nearest neighbors in 2ro radius
                %I = nearestneighbour(idx, data_t, 'Radius', 2*ro); % in order of distance
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
                            if dist < min_dist
                                min_dist = dist;
                                smallest_triangle = [idx, i , j];
                            end
                        end
                    end
                end
                if ~ismember(smallest_triangle,0)
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
        % TODO check if triangle has normal facing outward

    % ball pivoting
        active_index = 1;
        % while we have active edges, select next edge
        while active_index <= length(active_edges)
            % find closest pivoted point from edge
            A = data(active_edges(active_index,1),:);
            B = data(active_edges(active_index,2),:);
            C = data(active_edges(active_index,3),:);
            midpoint = (A+B)./2;
            %neighbors = nearestneighbour(midpoint', data_t, 'Radius', 2*ro);
            neighbors = rangesearch(KD_treesearcher, midpoint, 2* ro);
            neighbors = cell2mat(neighbors);
            disp(['Has ', num2str(length(neighbors)), ' neighbors.']);
            % exclude those that have been used (glued)
            evaluateThese = neighbors;
            evaluateThese(ismember(evaluateThese, [used_indices(1:used_indices_counter), active_edges(active_index,:)])) = [];
            if isempty(evaluateThese)
                disp('No valid points within reach, adding to boundary');
                boundary_edges(boundary_counter+1,:) = [active_edges(active_index,1), active_edges(active_index,2)];
                boundary_counter = boundary_counter+1;
            else
                disp(['Has ', num2str(length(evaluateThese)), ' left to be evaluated.']);
                min_point = 0;
                % interx was implemented by Hrishi Shah
                % see the function for license notes
                solution0 = interx(A, B, C, ro, ro, ro, 0);
                solution1 = interx(A, B, C, ro, ro, ro, 1);
                if pdist([solution0(1:3)' ; 0, 0, 0]) > pdist([solution1(1:3)'; 0, 0, 0])
                    active_edge_center = solution0(1:3)' - midpoint;
                else
                    active_edge_center = solution1(1:3)' - midpoint;
                end
                centers = zeros(length(evaluateThese), 3);
                evaluatePoints = data(evaluateThese, :);
                parfor i=1:length(evaluateThese)
                    solution0 = interx(A, B, evaluatePoints(i,:), ro, ro, ro, 0);
                    solution1 = interx(A, B, evaluatePoints(i,:), ro, ro, ro, 1);
                    if ~isnan(solution0)
                        if pdist([solution0(1:3)' ; 0, 0, 0]) > pdist([solution1(1:3)'; 0, 0, 0])
                            centers(i, :) = solution0(1:3)' - midpoint;
                        else
                            centers(i, :) = solution1(1:3)' - midpoint;
                        end
                        if  ~isempty(cell2mat(rangesearch(data(neighbors,:), centers(i, :), ro)))
                            centers(i, :) = [NaN, NaN, NaN];
                        end
                    else
                        centers(i, :) = [NaN, NaN, NaN];
                    end
                end
                % angles and minimum angle neighbor selection
                if ~isempty(centers(~isnan(centers)))
                    active_edge_matrix = repmat(active_edge_center, length(evaluateThese), 1);
                    ae_norm = sqrt(sum((active_edge_center).^2,2));
                    c_norm_vec = sqrt(sum((centers).^2,2));
                    cos_diff = diag(active_edge_matrix * centers') ./ (ae_norm * c_norm_vec);
                    degrees_diff = acos(cos_diff)*180/pi;
                    [~, ind] = min(degrees_diff);
                    min_point = evaluateThese(ind);
                end

                 if min_point ~= 0
                     disp(['Found minimal point: ', num2str(min_point)]);
                     tri(tri_counter+1,:) = [active_edges(active_index,1), min_point, active_edges(active_index,2)];
                     tri_counter = tri_counter+1;
                     if ~ismember(front, min_point)
                         disp('Adding point to front and active edges');
                         % is a new point so we don't have to glue
                         active_edges = [active_edges; active_edges(active_index,1), min_point, active_edges(active_index,2); min_point, active_edges(active_index,2), active_edges(active_index,1)];
                         front = [front, min_point];
                     else
                         disp('Glueing');
                         % we'll have to glue
                         output = quickfind2d(active_edges(:,[1,2]), min_point, active_edges(active_index,1));
                         if isempty(output)
                             active_edges = [active_edges; active_edges(active_index,1), min_point, active_edges(active_index,2)];
                         elseif ismember(front, active_edges(active_index,1))
                             idx = front(front == active_edges(active_index,1));
                             used_indices(used_indices_counter+1) = idx;
                             used_indices_counter = used_indices_counter+1;
                             toBeEvaluated(idx) = 0;
                             front(front == active_edges(active_index,1)) = [];
                         end
                         output = quickfind2d(active_edges(:,[1,2]), min_point, active_edges(active_index,2));
                         if isempty(output)
                             active_edges = [active_edges; min_point, active_edges(active_index,2), active_edges(active_index,1) ];
                         elseif ismember(front, active_edges(active_index,2))
                             idx = front(front == active_edges(active_index,2));
                             used_indices(used_indices_counter+1) = idx;
                             used_indices_counter = used_indices_counter+1;
                             toBeEvaluated(idx) = 0;
                             front(front == active_edges(active_index,2)) = [];
                         end
                     end
                 else
                     disp('Could not find point within reach, adding to boundary');
                     boundary_edges(boundary_counter+1,:) = [active_edges(active_index,1), active_edges(active_index,2)];
                     boundary_counter = boundary_counter+1;
                 end
                % find ball pivots for each neighboring point and determine
                % if point is unused or in the front
                % then 
                 % create triangle, try to get the normal to point outward
                 % add new edges to active edge list
                 % if edge already exists in active edge list, delete it
                 % from active edge list instead

                % cannot find, designate edge as boundary
                active_index = active_index +1;
                disp('active index now:');
                disp(active_index);
            end
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
    
    trimesh(tri, data(:, 1), data(:, 2), data(:, 3));
end

function output = quickfind2d(x,u,v)
    % quick intersect taken from http://www.mathworks.com/matlabcentral/answers/24414-find-row-with-certain-values
    [a, ~]=find(x==u);
    [c, ~]=find(x==v);
    index = c.*NaN;
    for kk=1:length(c)
      check =  a(a==c(kk))';
      if ~isempty ( check )
        index(kk) = check;
      end
    end
    output = index(~isnan(index));
    % quick intersect taken from http://www.mathworks.com/matlabcentral/answers/24414-find-row-with-certain-values
end
