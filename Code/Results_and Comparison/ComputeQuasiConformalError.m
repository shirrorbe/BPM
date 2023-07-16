function [qc_error, min_qc, max_qc, avg_qc] = ComputeQuasiConformalError(Vz,Vw,F)

% F = mesh.T;
% Vz = mesh.X;


p1 = Vz(F(:,1),:); p2 = Vz(F(:,2),:); p3 = Vz(F(:,3),:);

q1 = Vw(F(:,1),:); q2 = Vw(F(:,2),:); q3 = Vw(F(:,3),:);

% compute edge vectors
u1 = p2 - p1;
u2 = p3 - p1;

v1 = q2 - q1;
v2 = q3 - q1;

% compute orthonormal bases
e1 = normalize(u1,2,'norm');
e2 = normalize(u2 - dot(u2,e1,2) .* e1,2,'norm');

f1 = normalize(v1,2,'norm');
f2 = normalize(v2 - dot(v2,f1,2) .* f1,2,'norm');

% project onto bases
p1 = zeros(size(p1));
p2 = [dot(u1,e1,2), dot(u1,e2,2), zeros(size(u1,1),1)];
p3 = [dot(u2,e1,2), dot(u2,e2,2), zeros(size(u2,1),1)];

q1 = zeros(size(q1,1),3);
q2 = [dot(v1,f1,2), dot(v1,f2,2), zeros(size(f1,1),1)];
q3 = [dot(v2,f1,2), dot(v2,f2,2), zeros(size(f2,1),1)];


A = 2 * normrow(cross(u1, u2, 2));

Ss = (q1.*(p2(:,2) - p3(:,2)) + q2.*(p3(:,2) - p1(:,2)) + q3.*(p1(:,2) - p2(:,2)))./A;
St = (q1.*(p3(:,1) - p2(:,1)) + q2.*(p1(:,1) - p3(:,1)) + q3.*(p2(:,1) - p1(:,1)))./A;
a = dot(Ss,Ss,2);
b = dot(Ss,St,2);
c = dot(St,St,2);
det = sqrt((a-c).^2 + 4.*(b.^2));
Gamma = sqrt(0.5*(a + c + det));
gamma = sqrt(0.5*(a + c - det));

Gamma_ordered = Gamma;
Gamma_ordered(Gamma < gamma) = gamma(Gamma < gamma);

gamma_ordered = gamma;
gamma_ordered(Gamma < gamma) = Gamma(Gamma < gamma);

qc_error = Gamma_ordered./gamma_ordered;
min_qc = min(qc_error);
max_qc = max(qc_error);

[f_areas,~] = calculate_face_vertex_areas  (Vz,F);
avg_qc = sum(qc_error.*f_areas) / sum(f_areas);
end

