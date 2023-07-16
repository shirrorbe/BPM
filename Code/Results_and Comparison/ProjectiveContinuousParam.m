function [Mw_fine, Mz_fine] = ProjectiveContinuousParam(Vz, Vw, F ,n_subd)
  % split the mesh into seperate triangles to render each triangle seperately
        nf = size(F,1); nv = size(Vz,1);
        Vz_new = Vz(reshape(F',[],1),:); 
        Vw_new = Vw(reshape(F',[],1),:); 
        F_new = (reshape((1:(3*nf)),[3,nf]))';        
        u = ComputeApproximatedScaleFactors_faces(Vw,Vz,F);

            
        % sampling the domain by subdividing the mesh
        if n_subd > 0
            [Vz1, F1] = meshSubdivision(Vz_new, F_new);
%             [Vz1, F1] = meshSubdivision([Vz_new,zeros(size(Vz_new,1),1)], F_new);
            for i=1:(n_subd-1)
                [Vz1, F1] = meshSubdivision(Vz1, F1);
            end
            orig_F = repmat((1:size(F_new,1)),4^n_subd,1);
            orig_F = orig_F(:);
        else
            Vz1 = Vz_new; F1=F_new;
            orig_F = (1:size(F_new,1))';
        end
        orig_F = repmat(orig_F,3,1); % the original face of each new face

        Mz_fine.V = Vz1; Mz_fine.F = F1;
        u_new = u(orig_F,:);
        u1_i = u_new(:,1); u1_j = u_new(:,2); u1_k = u_new(:,3);
    
%         u1_i = ui(orig_F);u1_j = uj(orig_F);u1_k = uk(orig_F);

        zi = complex(Vz_new(F_new(:,1),1),Vz_new(F_new(:,1),2)); 
        zj = complex(Vz_new(F_new(:,2),1),Vz_new(F_new(:,2),2));
        zk = complex(Vz_new(F_new(:,3),1),Vz_new(F_new(:,3),2));
        
        wi = complex(Vw_new(F_new(:,1),1),Vw_new(F_new(:,1),2)); 
        wj = complex(Vw_new(F_new(:,2),1),Vw_new(F_new(:,2),2));
        wk = complex(Vw_new(F_new(:,3),1),Vw_new(F_new(:,3),2));
        

        Vz1_faces = Vz1([F1(:,1);F1(:,2);F1(:,3)],:);

        z1 = complex(Vz1_faces(:,1),Vz1_faces(:,2));
        z1_i = zi(orig_F); z1_j = zj(orig_F); z1_k = zk(orig_F);
        w1_i = wi(orig_F); w1_j = wj(orig_F); w1_k = wk(orig_F);

        [w] = arrayfun(@ProjectiveInterp, z1 ,z1_i, z1_j, z1_k, w1_i,w1_j,w1_k,u1_i,u1_j,u1_k);
        
        Vw1=zeros(size(Vz1));      
        Vw1(F1(:),:) = [real(w) imag(w) zeros(size(w,1),1)];
        Mw_fine.V = Vw1; Mw_fine.F = F1;


   
end

