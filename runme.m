% BPM: Blended Piecewise Mobius

startup

% Global parameters
n_subd = 4;
texture_filename = 'texture.png';

in_dir = 'Data\Inputs\';

in_data_files = dir([in_dir,'*.obj']);
out_dir = 'Data\data_out\';

name_ext = '-output';


for i = 1:size(in_data_files,1)

    RunBPM(out_dir,in_data_files(i),n_subd,name_ext)
    
end



