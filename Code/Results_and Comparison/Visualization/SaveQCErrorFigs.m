function f = SaveQCErrorFigs(mesh, error,climit,cm_data ,filename, save_figs_flag, show_figs_flag)
    f = figure('visible','off');
    h = function_visualization_on_mesh(mesh.V,mesh.F,error);
    set(h,'edgecolor','none')
    clim(climit)
    set(colorbar,'visible','off')
    colormap(cm_data)
    axis off
    set(gcf,'color','w')
    
    if save_figs_flag
        export_fig(filename,'-nocrop','-m3');
%         savefig(gcf,[filename,'_fig.fig'])
    end
   
    if show_figs_flag
        figure(f)
        set(gcf, 'WindowStyle', 'docked')
        cameratoolbar
    end

 


end