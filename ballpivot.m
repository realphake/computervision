function  tri  = ballpivot( data )
% POINTCLOUD2TRIANGLES 
% delaunay doesn't work properly, trying ball pivoting algorithm
% patented by IBM, does this matter?
% http://www.research.ibm.com/vistechnology/pdf/bpa_tvcg.pdf
% ro = radius of the ball that pivots
ro = 0.002; % somewhere between 0.001 and 0.002
data = data(:, 1:3);
data_t = data';
% find seed triangle
    % pick random point
    idx = ceil(rand * size(data,1));
    % find nearest neighbors in 2ro distance
    I = nearestneighbour(idx, data_t, 'Radius', 2*ro); % in order of distance
    % keep queue of active edges (i.e. edges that have not been pivoted around yet 
    front = [idx, I(1:2)];
    active_edges = [front([1, 2]); front([2, 3]); front([3, 1])];
    used_indices = [];
    % create triangles with the random point where distance between two
    % neighbors is < 2ro
    tri = [idx, I(1:2)];
    % TODO check if triangle has normal facing outward

% ball pivoting
    % keep list of boundary edges
    boundary_edges = [];
    % keep list of frozen edges
    frozen_edges = [];
    active_index = 1;
    % while we have active edges, select next edge
    while active_index <= length(active_edges)
        % find closest pivoted point from edge
        midpoint = (data(active_edges(active_index,1),:)+data(active_edges(active_index,2),:))./2;
        neighbors = nearestneighbour(midpoint', data_t, 'Radius', 2*ro);
        % exclude those that have been used (glued)
        [~, ia] = intersect(neighbors, used_indices);
        neighbors(ia) = [];
        % find ball pivots for each neighboring point and determine
            % if point is unused or in the front
            % then 
                % create triangle, try to get the normal to point outward
                % add new edges to active edge list
                % if edge already exists in active edge list, delete it
                % from active edge list instead
            
            % cannot find, designate edge as boundary
        active_index = active_index +1;
    end

end

