function int = FireInitializer(locations, n, IC)
% FIREINITIALIZER  Initial temperature field with circular spot fires.
%   int = FireInitializer(locations, n, IC) returns an n-by-n grid on the
%   unit square (rows = y, columns = x, matching meshgrid) that is IC inside
%   each spot fire and 0 elsewhere.
%   locations: k-by-3 matrix, one row per fire: [x_center, y_center, radius].
x = linspace(0, 1, n);
[X, Y] = meshgrid(x, x);
int = zeros(n);
for i = 1:size(locations, 1)
    r2 = (X - locations(i,1)).^2 + (Y - locations(i,2)).^2;
    int = max(int, IC * (r2 <= locations(i,3)^2));
end
end
