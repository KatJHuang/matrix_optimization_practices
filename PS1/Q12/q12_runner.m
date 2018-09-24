n_points = 201;
lo_bound = -3;
hi_bound = 3;

dx = (hi_bound - lo_bound) / (n_points - 1);

t = linspace(lo_bound, hi_bound, n_points)
x = [3, 2]';

v0 = [-2 -4]'; % offset
v = [-2 5]'; % slope

A = v .* t + v0
figure; hold on
axis equal

plot(x(1, :), x(2, :), 'ro', 'MarkerSize', 5);
plot(A(1, :), A(2, :));

% l2-norms
l2_norm = sqrt(dot(A - x, A - x));
[~, pos] = min(l2_norm);

plot(A(1, pos), A(2, pos), 'bo', 'MarkerSize', 5);
circle(x(1, :), x(2, :), l2_norm(pos))

% L1-norm  
l1_norm = sum(abs(A - x), 1);
[~, pos] = min(l1_norm);

plot(A(1, pos), A(2, pos), 'bo', 'MarkerSize', 5);
diamond(x(1, :), x(2, :), l1_norm(pos))

% L-inf norm
l_inf_norm = max(abs(A - x));
[min_dist, pos] = min(l_inf_norm);
plot(A(1, pos), A(2, pos), 'bo', 'MarkerSize', 5);
square(x(1, :), x(2, :), l_inf_norm(pos))
%% Linear Programming
%%
norms = [1 2 inf];
for ii = 1:size(norms, 2)
    fprintf('L-%d norm\n', norms(ii));
    [y, r] = proj_cvx (x, v0 , v, norms(ii));
    y
end
%%
function circle(x, y, r)
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)
    theta = 0:0.01:2 * pi; 
    xp = r * cos(theta);
    yp = r * sin(theta);
    plot(x + xp, y + yp);
end

function diamond(x, y, r)
% by diamond I mean a rhombus of equal sides, which is essentially
% the shape of L1 norm ball.
% x, y - centre of the diamond
% r - distance between the centre to any of the vertices
    vertices = [x-r y; x y-r; x+r y; x y+r];
    poly = polyshape(vertices);
    plot(poly);
end

function square(x, y, r)
% x, y - centre of the square
% r - distance from centre to any side. Equal to half of the side length
    side_length = 2 * r;
    rectangle('Position', ...
        [x - side_length/2, y - side_length/2, side_length, side_length]);
end

function [y, r] = proj_cvx (x, v0 , v, nrm)%% x, v0 and v must be column vectors
    objtv = @(y) norm (x-y, nrm); %% objective is L_2 norm
    cvx_begin
        variable y(2) %% 2-d variable we are optimizing over
        variable t(1) %% real valued parameter that defines
        minimize ( objtv (y)) %% defining the objective
        subject to
        v0 + t*v == y %% the projection y must be in set A
    cvx_end
    r = objtv (y); %% minimum value of the objective
end