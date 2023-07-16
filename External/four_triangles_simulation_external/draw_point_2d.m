function draw_point_2d(pt,varargin);

if size(pt,1)==0
    return;
end

h = plot(pt(:,1),pt(:,2),'.','MarkerSize',20);
if length(varargin)>0
    set(h, varargin{:});
end