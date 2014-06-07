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
