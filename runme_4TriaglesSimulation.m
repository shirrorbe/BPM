startup

% create 4 triangles setup - Mz
vi = [0 0 0]; vj = [1 0 0]; vk = [0.5 cos(pi/6) 0]; 
vl = [1.5 cos(pi/6) 0]; vm = [-0.5 cos(pi/6) 0]; vn = [0.5 -cos(pi/6) 0];
Vz = [vi;vj;vk;vl;vm;vn];
F = [1 2 3; 3 2 4; 5 1 3; 1 6 2];
Mz = []; Mz.V = Vz; Mz.F = F;

% set the initial transform - Mw
vi_tr = vi; vj_tr = vj; vk_tr = vk; 
vl_tr = vl; vm_tr = vm; vn_tr = vn;
Vw = [vi_tr;vj_tr;vk_tr;vl_tr;vm_tr;vn_tr];
Mw = []; Mw.V = Vw; Mw.F = F;

% simulation parameters
n_subd = 4;
show_QCerr = false;
show_Wt = false;
Bijectivity_interactive(Mz,Mw,n_subd,show_QCerr,show_Wt);
cameratoolbar

% Interactive simultaion instructions:
% To move a point, left click on it and move the cursor. when getting to 
% the required location, press the right bottom.



