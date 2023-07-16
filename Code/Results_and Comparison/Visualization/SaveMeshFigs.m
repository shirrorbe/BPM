function f = SaveMeshFigs(mesh, filename, save_figs_flag, show_figs_flag)
    
    f = figure('visible','off');
    show_mesh(mesh); set(gcf,'color','w')
    axis equal off; 


    if save_figs_flag
        export_fig(filename,'-nocrop','-m3');
    end
   
    if show_figs_flag
        figure(f)
        set(gcf, 'WindowStyle', 'docked')
        cameratoolbar
    end

end