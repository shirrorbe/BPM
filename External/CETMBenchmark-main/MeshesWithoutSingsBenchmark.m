%% CETM generating benchmark from meshes without singularities
clear all
close all
clc

%CutMeshFromSingularity executable needs to be in the main drectory of this
%script


%inputs
%%%%%change to directory of input OBJ files
%datafolder = ['.' filesep 'OdedBenchmark' filesep 'No_UVs' filesep 'Uncut'];
datafolder = ['.' filesep 'ToyBenchmark'];
%%%%%change to directory of output files
resultsFolder = ['.' filesep 'ToyBenchmark' filesep 'results'];

files = dir([datafolder filesep '*.obj']);
benchMeshes = {files.name}';
benchMeshes = cellfun(@(f) [datafolder filesep f], benchMeshes, 'UniformOutput', false);
success=zeros(length(benchMeshes),1);
%success variable shows who succeeded
for meshNum=1:length(benchMeshes)
    meshName=benchMeshes{meshNum};
    [V, T, ~,~] = readObj(meshName);
    mesh = create_mesh(V,T);
    mesh.name = meshName;
    meshName
    meshNum
    phidiff = 1; % Threshold for iterative singularity placement
    
    % Place singularities and compute PHI
    [locS, KS, PHI] = find_singularities(mesh, phidiff);
    locS=locS';
    if ~isempty(intersect(mesh.BV, locS))
        disp([meshName, ': singularity on boundary']);
        continue;
    end
    rawMeshFullDir = split(meshName,filesep);
    onlyName = rawMeshFullDir{end};
    onlyName = onlyName(1:end-4);
    singFileName = [resultsFolder filesep onlyName '.sing'];
    write_singularities(singFileName, locS);
    %figure; show_func(mesh,PHI);
    ufixed = [find(mesh.BV) zeros( size(find(mesh.BV),1), 1)];
    %theta = ones(size(V,1),1)*2*pi;
    %theta(locS) = 2*pi-KS;
    KPres = zeros(size(V,1),1);
    KPres(locS) = KS;
    %[y, nbrokentris, lcrerr, flaterr]= dcflatten_wrap(T, size(V,1), sqrt( meshFaceEdgeLen2s(V, T) ), theta, ufixed);
    [UV,u,ELnew,Knew, success(meshNum)] = cetm(meshName, singFileName, resultsFolder, V,T, KPres);
    sumSuccess=sum(success)
end