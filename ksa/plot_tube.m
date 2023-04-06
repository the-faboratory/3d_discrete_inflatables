function [x, y, z] = plot_tube(curve, radius, n,  join_radius, face_alpha, varargin) 
% 
% plot_tube constructs a tube, or warped cylinder, along any 3D curve, much like the build in cylinder function. If no output are requested, the tube is plotted. Otherwise, you can plot by using surf(x, y, z);

% Janus H. Wesenberg, july 2004
% Edited by Robert Baines, 2022

% Function parser described here https://www.mathworks.com/help/matlab/matlab_prog/parse-function-inputs.html
% In brief: add[type](inputParser,name,check function)
p = inputParser;
addRequired(p, 'curve', @isnumeric)
addRequired(p, 'radius', @isnumeric) % Positional, required
addRequired(p, 'n', @isnumeric) % Positional, required
addParameter(p, 'join_radius', 0.1, @isnumeric) % Not positional, marked by a flag
addParameter(p, 'face_alpha', 1, @isnumeric)

parse(p, curve, radius, n, varargin{:})
join_radius = p.Results.join_radius;
face_alpha = p.Results.face_alpha;

% Join nearby points
n_points = 1;
for k = 2:(size(curve, 2) - 1)
  if norm(curve(:, k) - curve(:, n_points))>join_radius;
    n_points = n_points + 1;
    curve(:, n_points) = curve(:, k);
  end
end

% Always include endpoint
if norm(curve(:, end) - curve(:, n_points)) > 0
  n_points = n_points + 1;
  curve(:, n_points) = curve(:, end);
end

% Deltavecs: average for internal points. First stretch for endpoitns.
tangent = curve(:, [2:end, end]) - curve(:, [1, 1:end - 1]);

% Make psuedo_normal not parallel to tangent(:, 1)
psuedo_normal = zeros(3, 1);
[buf, idx] = min(abs(tangent(:, 1)));
psuedo_normal(idx) = 1;

xyz = repmat([0], [3, n + 1, n_points + 2]);

% Precalculate cos and sine factors:
cfact = repmat(cos(linspace(0, 2*pi, n + 1)), [3, 1]);
sfact = repmat(sin(linspace(0, 2*pi, n + 1)), [3, 1]);

% Main loop: propagate the normal (psuedo_normal) along the tube
for k = 1:n_points
  binormal = cross(psuedo_normal, tangent(:, k));
  binormal = binormal./norm(binormal);
  psuedo_normal = cross(tangent(:, k), binormal);
  psuedo_normal = psuedo_normal./norm(psuedo_normal);
  % Update xyz:
  xyz(:, :, k + 1) = repmat(curve(:, k), [1, n + 1]) + cfact.*repmat(radius*psuedo_normal, [1, n + 1]) + sfact.*repmat(radius*binormal, [1, n + 1]);
end;

% Cap the ends
xyz(:, :, 1) = repmat(curve(:, 1), [1, n + 1]);
xyz(:, :, end) = repmat(curve(:, end),[1, n + 1]);

% Extract results
x = squeeze(xyz(1, :, :));
y = squeeze(xyz(2, :, :));
z = squeeze(xyz(3, :, :));

% Plot
if nargout<3
    hSurface=surf(x, y, z,'FaceColor',[1 0 0],'FaceLighting','gouraud','EdgeColor','none');
  hold on
  set(hSurface,'FaceColor',[1 0 0], ...
      'FaceAlpha',1,'FaceLighting','gouraud','EdgeColor','none')
    garbage = 0;
end;
