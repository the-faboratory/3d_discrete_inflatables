function [L] = cumulative_arc_len(X)
% Radius of curvature and curvature vector for 2D or 3D curve
%   X:   2 or 3 column array of x, y (and possibly z) coordiates
%   L:   Cumulative arc length


N = size(X,1);
dims = size(X,2);
if dims == 2
    X = [X,zeros(N,1)];  % Do all calculations in 3D
end
L = zeros(N,1);

for i = 2:N-1
    L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));
end
i = N;
L(i) = L(i-1)+norm(X(i,:)-X(i-1,:));

end

