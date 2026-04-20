function plot_circle(x, y, r, varargin)
    theta = linspace(0, 2*pi, 100);
    plot(x + r*cos(theta), y + r*sin(theta), varargin{:});
end