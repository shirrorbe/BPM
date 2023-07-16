function Bijectivity_interactive(Mz,Mw,n_subd,show_QCerr,show_Wt)
Vz = Mz.V; Vw = Mw.V; F=Mz.F;
% [Mw_fine, Mz_fine,t_weights_ijkt] = BPMContinuousParam(Vz, Vw, F ,n_subd);
[Mw_fine, Mz_fine,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMContinuousParam(Vz, Vw, F ,n_subd);

Vz_fine = Mz_fine.V; Vw_fine = Mw_fine.V; F_fine = Mz_fine.F;
nf_fine = size(F_fine,1); nff1 = floor(nf_fine/4);
% original mesh figure
figure
title('Mz');
patch('Faces',F_fine,'Vertices',Vz_fine(:,1:2),'FaceColor','w'); axis equal; hold on; 
draw_point_2d([Vz(:,1),Vz(:,2)],'MarkerSize',30);
aa = axis; aa = aa + 0.5*[-1,1,-1,1]; axis(aa);
set(gcf,'WindowStyle','docked')
if show_Wt
    % show weights 
    w_ind = 1; % weights of Mij
    hp_wt=patch('Faces',Mz_fine.F,'Vertices',Mw_fine.V,'FaceVertexCData',full(t_weights_ijkt(:,w_ind)),'FaceColor','interp','EdgeColor','none');
    colorbar
end


% transformed mesh figure
figure('WindowButtonDownFcn',@wbdcb,'Interruptible','off','BusyAction','cancel')
% figure('WindowButtonDownFcn',@wbdcb2,'Interruptible','off','BusyAction','cancel')
title('Mw');
hp1 = patch('Faces',F_fine,'Vertices',Vw_fine(:,1:2),'FaceColor','none'); axis equal; hold on; 
hf1 = plot(Vw(:,1),Vw(:,2),'.','MarkerSize',30);


prec='%.2f';
disp(['(' num2str(Vw(1,1),prec) ',' num2str(Vw(1,2),prec) ')',...
        '(' num2str(Vw(2,1),prec) ',' num2str(Vw(2,2),prec) ')',...
        '(' num2str(Vw(3,1),prec) ',' num2str(Vw(3,2),prec) ')',...
        '(' num2str(Vw(4,1),prec) ',' num2str(Vw(4,2),prec) ')',...
        '(' num2str(Vw(5,1),prec) ',' num2str(Vw(5,2),prec) ')',...
        '(' num2str(Vw(6,1),prec) ',' num2str(Vw(6,2),prec) ')']);

% 
% % selection point
% spnt = [0.5, 0.5]; f=1;
% pp_select = plot(spnt(1),spnt(2),'MarkerSize',10,'Marker','.','MarkerFaceColor','r');


if show_QCerr
    % Quasi conformal error
    limits_qc =[1 3];% [1 1.1];
    [qc_error, ~, ~, ~] = ComputeQuasiConformalError(Mz_fine.V,Mw_fine.V,Mz_fine.F);
    hp_QCerr=patch('Faces',Mz_fine.F,'Vertices',Mw_fine.V,'FaceVertexCData',full(qc_error),'FaceColor','flat','EdgeColor','none');
    clim(limits_qc)
    colorbar
end



aa = axis; aa = aa + 0.5*[-1,1,-1,1]; axis(aa);

% draw original ijk edges
line([Vz(1,1); Vz(2,1)],[Vz(1,2); Vz(2,2)], 'Color',[0 0 1],'LineWidth',1)
line([Vz(1,1); Vz(3,1)],[Vz(1,2); Vz(3,2)], 'Color',[0 0 1],'LineWidth',1)
line([Vz(2,1); Vz(3,1)],[Vz(2,2); Vz(3,2)], 'Color',[0 0 1],'LineWidth',1)


set(gcf,'WindowStyle','docked')

x = Vw(:,1:2);
w=Vw(:,1)+1i*Vw(:,2);
n = length(w);

l1 = 0;


function update_location(w)  
%     disp('update_location')
    x = [real(w), imag(w)];
    Vw = [real(w), imag(w),zeros(size(w,1),1)];
            
    
    % plot weights near the point
    [Mw_fine, Mz_fine,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMContinuousParam(Vz, Vw, F ,n_subd);
    Vz_fine = Mz_fine.V; Vw_fine = Mw_fine.V; F_fine = Mz_fine.F;

    
    % check bijectivity
    Vw_fine_i = Vw_fine(F_fine(:,1),:);
    Vw_fine_j = Vw_fine(F_fine(:,2),:);
    Vw_fine_k = Vw_fine(F_fine(:,3),:);
    
    tri_normals = cross(Vw_fine_j-Vw_fine_i,Vw_fine_k-Vw_fine_i);
    flipped = tri_normals(:,3)<0;
    if sum(flipped>0)
        disp('Bijectivity violated!')
    end
    set(hp1,'Vertices',Vw_fine(:,1:2));

     % color flipped trianges in red
    fc = ones(size(F_fine,1),3); fc(flipped,2:3) = 0*fc(flipped,2:3);


    set(hf1,'XData',Vw(:,1),'YData',Vw(:,2));drawnow
    set(hp1,'FaceVertexCData',fc,'FaceColor','flat');drawnow

    set(hp1,'FaceAlpha',0.5);drawnow

%         % update chosen points
%         for ii=1:size(chosen_pnts,2)
%             ind = chosen_pnts(ii);
%             set(hPntDat(ii),'XData',Vw_fine(ind,1),'YData',Vw_fine(ind,2));drawnow
%             set(hPntTxt(ii),'Position',[Vw_fine(ind,1)+0.01,Vw_fine(ind,2)]);        
%         end
    
%     % update location of the 4 triangles points
%     set(ht1,'Position',[Vw(1,1)+0.01,Vw(1,2)]);
%     set(ht1,'string',['(' num2str(Vw(1,1),prec) ',' num2str(Vw(1,2),prec) ')']);
% 
%     set(ht2,'Position',[Vw(2,1)+0.01,Vw(2,2)]);
%     set(ht2,'string',['(' num2str(Vw(2,1),prec) ',' num2str(Vw(2,2),prec) ')']);
%     
%     set(ht3,'Position',[Vw(3,1)+0.01,Vw(3,2)]);
%     set(ht3,'string',['(' num2str(Vw(3,1),prec) ',' num2str(Vw(3,2),prec) ')']);
%      
%     set(ht4,'Position',[Vw(4,1)+0.01,Vw(4,2)]);
%     set(ht4,'string',['(' num2str(Vw(4,1),prec) ',' num2str(Vw(4,2),prec) ')']);
%     
%     set(ht5,'Position',[Vw(5,1)+0.01,Vw(5,2)]);
%     set(ht5,'string',['(' num2str(Vw(5,1),prec) ',' num2str(Vw(5,2),prec) ')']);
%     
%     set(ht6,'Position',[Vw(6,1)+0.01,Vw(6,2)]);
%     set(ht6,'string',['(' num2str(Vw(6,1),prec) ',' num2str(Vw(6,2),prec) ')']);
     
    if show_QCerr
        %QC error
        [qc_error, ~, ~, ~] = ComputeQuasiConformalError(Mz_fine.V,Mw_fine.V,Mz_fine.F);
        set(hp_QCerr,'Vertices',Vw_fine(:,1:2));
        set(hp_QCerr,'FaceVertexCData',full(qc_error),'FaceColor','flat');
        clim(limits_qc)
        colorbar
    end
% 
%     % update the location of the seclected point
% %     [Mw_fine, Mz_fine,t_weights_ijkt] = BPMContinuousParam(Vz, Vw, F ,n_subd);
%     [w,t_weights_ijkt,Mijk,Mij,Mjk,Mki] = BPMMatrixInterpolator_wrap([spnt(1),spnt(2)],Vz, Vw,f);
% %         [w,t_weights_ijkt] = BPMMatrixInterpolator_wrap([spnt(1),spnt(2)],Vz, Vw,f);
% 
%     set(pp_select,'XData',real(w),'YData',imag(w));drawnow
%     Mijk
%     disp(['weight Mt=Mijk: ',num2str(t_weights_ijkt(4))])
%     Mij
%     disp(['weight Mij: ',num2str(t_weights_ijkt(1))])
%     Mjk
%     disp(['weight Mjk: ',num2str(t_weights_ijkt(2))])
%     Mki
%     disp(['weight Mki: ',num2str(t_weights_ijkt(3))])
% 
%     disp(['chosen point loc: ',num2str(w),' w: ',num2str(t_weights_ijkt)])
%    

end

function wbdcb2(src,evnt)
%     disp('wbdcb2')
%     if strcmp(get(src,'SelectionType'),'normal') %left click
%         set(src,'pointer','hand')
%         set(src,'WindowButtonMotionFcn',@wbmcb2)
%         set(src,'WindowButtonDownFcn',@wbdcb3)
%     end
    if strcmp(get(src,'SelectionType'),'normal') %left click
        set(src,'pointer','hand')
        set(src,'WindowButtonMotionFcn',@wbdcb)
        set(src,'WindowButtonDownFcn',@wbdcb)
    end
end

function wbdcb3(src,evnt)
%         disp('wbdcb3')
        if strcmp(get(src,'SelectionType'),'normal') %left click
            set(src,'pointer','hand')
            cp = get(gca,'CurrentPoint');
            spnt = [cp(1,1), cp(1,2)];
            set(pp_select,'XData',spnt(1),'YData',spnt(2));drawnow
            set(src,'WindowButtonDownFcn',@wbdcb)
            set(src,'WindowButtonUpFcn',@wbucb2) % mouse bottom released

        end
    end

function wbdcb(src,evnt)

        set(src,'WindowButtonMotionFcn',@wbdcb2)

        if strcmp(get(src,'SelectionType'),'normal') %left click
            set(src,'pointer','circle')
            cp = get(gca,'CurrentPoint');
            dd = normv(x - repmat(cp(1,1:2),n,1));
            [m,l1] = min(dd);
            w(l1) = cp(1,1) + 1i*cp(1,2);
            update_location(w);
            set(src,'WindowButtonMotionFcn',@wbmcb)% mouse moves
            set(src,'WindowButtonUpFcn',@wbucb) % mouse bottom released
        end 

              
end

function wbmcb(src,evnt)
%         disp('wbmcb')
        cp = get(gca,'CurrentPoint');
        w(l1) = cp(1,1) + 1i*cp(1,2);
        update_location(w);
end

function wbucb(src,evnt)
%         disp('wbucb')
        disp(['(' num2str(Vw(1,1),prec) ',' num2str(Vw(1,2),prec) ')',...
        '(' num2str(Vw(2,1),prec) ',' num2str(Vw(2,2),prec) ')',...
        '(' num2str(Vw(3,1),prec) ',' num2str(Vw(3,2),prec) ')',...
        '(' num2str(Vw(4,1),prec) ',' num2str(Vw(4,2),prec) ')',...
        '(' num2str(Vw(5,1),prec) ',' num2str(Vw(5,2),prec) ')',...
        '(' num2str(Vw(6,1),prec) ',' num2str(Vw(6,2),prec) ')']);

        set(src,'Pointer','arrow')
        set(src,'WindowButtonMotionFcn','')
        set(src,'WindowButtonUpFcn','')
end

function wbmcb2(src,evnt)
%                 disp('wbmcb2')
                set(src,'Pointer','hand')
%                 set(src,'WindowButtonDownFcn',@wbdcb2)

%                 cp = get(gca,'CurrentPoint');
%                 w(l1) = cp(1,1) + 1i*cp(1,2);
%                 update_location(w);
      end

function wbucb2(src,evnt)
%     disp('wbucb2')
    in1 = inpolygon(spnt(1),spnt(2),Vw(F(1,[1 2 3 1]),1)',Vw(F(1,[1 2 3 1]),2)');
    in2 = inpolygon(spnt(1),spnt(2),Vw(F(2,[1 2 3 1]),1)',Vw(F(2,[1 2 3 1]),2)');
    in3 = inpolygon(spnt(1),spnt(2),Vw(F(3,[1 2 3 1]),1)',Vw(F(3,[1 2 3 1]),2)');
    in4 = inpolygon(spnt(1),spnt(2),Vw(F(4,[1 2 3 1]),1)',Vw(F(4,[1 2 3 1]),2)');
    f=in1+2*in2+3*in3+4*in4;

    set(src,'Pointer','arrow')
    set(src,'WindowButtonMotionFcn','')
    set(src,'WindowButtonUpFcn','')
    set(src,'WindowButtonDownFcn',@wbdcb)
end

function nn = normv(V)
        nn = sqrt(sum(V.^2,2));
end



%     function wbdcb(src,evnt)
%         disp('wbdcb')
%         if strcmp(get(src,'SelectionType'),'normal') %left click
%             set(src,'pointer','circle')
%             cp = get(gca,'CurrentPoint');
%             dd = normv(x - repmat(cp(1,1:2),n,1));
%             [m,l1] = min(dd);
%             w(l1) = cp(1,1) + 1i*cp(1,2);
%             update_location(w);
%             set(src,'WindowButtonMotionFcn',@wbmcb)
%             set(src,'WindowButtonUpFcn',@wbucb)
%         end
% 
%         function wbmcb(src,evnt)
%             disp('wbmcb')
%             cp = get(gca,'CurrentPoint');
%             w(l1) = cp(1,1) + 1i*cp(1,2);
%             update_location(w);
%         end
% 
%         function wbucb(src,evnt)
%             disp('wbucb')
%             set(src,'Pointer','arrow')
%             set(src,'WindowButtonMotionFcn','')
%             set(src,'WindowButtonUpFcn','')
%         end
%         
%         function nn = normv(V)
%             nn = sqrt(sum(V.^2,2));
%         end
%     end
end
