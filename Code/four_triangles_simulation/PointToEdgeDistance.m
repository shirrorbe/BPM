function d = PointToEdgeDistance(pnt, vi, vj)
    Vz_fine = pnt;
    a =  repmat(vi - vj,size(Vz_fine,1),1); b = Vz_fine - vj;
    d = normrow(cross(a,b)) ./ normrow(a);
end