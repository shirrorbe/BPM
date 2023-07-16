function h=show_mesh(mesh, vals, varargin)

p = inputParser;
    p.KeepUnmatched=true;
    addOptional(p,'FaceColor','w');
    addOptional(p,'EdgeColor','k');
    addOptional(p,'CW',[]);
    addOptional(p,'UVt',[]);
    addOptional(p,'FUVt',[]);
    addOptional(p,'Texture','cross.png');
    addOptional(p,'SingularPts',0);
    addOptional(p,'SingularPtsR',20);
    addOptional(p,'FaceAlpha',.8);
    addOptional(p,'EdgeAlpha',1);

%     addOptional(p,'locs',{[1:M.nf]'});
    addOptional(p,'LineWidth',.5);
    addOptional(p,'docked',1);

parse(p,varargin{:});


V = mesh.V; F = mesh.F;

colors = [1 1 1];
if nargin >= 2
    cmin = min(vals); cmax = max(vals);
    if cmax == cmin 
        colors = ones(size(F,1),3);
    else
        vals_inds = floor(255*(vals - cmin)/(cmax-cmin)) + 1;
        cmap = colormap(summer);
        colors = cmap(vals_inds,:)/255;
    end
       
end
h=patch('faces',F,'vertices',V, ...
                        'FaceColor',p.Results.FaceColor, ...
                        'EdgeColor',p.Results.EdgeColor, ...
                        'FaceAlpha',p.Results.FaceAlpha, ...
                        'linewidth',p.Results.LineWidth, ...
                        'EdgeAlpha',p.Results.EdgeAlpha);

% patch('vertices',V,'faces',F,'facecolor',colors,'edgecolor','r','facealpha',0.5); axis equal; cameratoolbar; 
cameratoolbar('SetCoordSys','none'); hold on