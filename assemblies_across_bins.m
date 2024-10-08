function [As_across_bins, As_across_bins_index] = assemblies_across_bins(assembly,BinSizes)
% Assemblies of all bins are rearranged in a different structure (no assembly 
% modification, only formal). 
%
% ARGUMENTS:
%   assembly  - structure containing assembly information
%   BinSizes  - list of investigated bin sizes
%
% RETURNS:
%   As_across_bins         - structure containing assembly information 
%   As_across_bins_index   - information to link assemblies in "As_across_bins"
%                          back to the structure "assembly" (output of Main_assemblies_detection.m):
%                          assembly As_across_bins{i} is contained in 
%                          assembly.bin{As_across_bins_index{i}(1)}.n{As_across_bins_index{i}(2)}.
%
%
%  � 2020 Russo
%  for information please contact eleonora.russo@zi-mannheim.de


empty=1;
e=1;
while empty
% for the first bin size tested
     if isfield(assembly.bin{e},'n')
        nAss=length(assembly.bin{e}.n);    
        As_across_bins=assembly.bin{e}.n;
        for i=1:nAss
            As_across_bins{i}.bin=BinSizes(e);
            As_across_bins_index{i}=[e,i];
        end
        j=nAss+1;    
        identity=1:nAss;
        id=nAss;
        empty=0;
     else
         e=e+1;
     end
end  


for gg=e+1:size(assembly.bin,2)
    if isfield(assembly.bin{gg},'n')
        nAss=length(assembly.bin{gg}.n);    
        for i=1:nAss             
            As_across_bins{j}.elements=assembly.bin{gg}.n{i}.elements;
            As_across_bins{j}.pr=assembly.bin{gg}.n{i}.pr;
            As_across_bins{j}.lag=assembly.bin{gg}.n{i}.lag;
            As_across_bins{j}.Time=assembly.bin{gg}.n{i}.Time;
            As_across_bins{j}.Noccurrences=assembly.bin{gg}.n{i}.Noccurrences;
            As_across_bins{j}.bin=BinSizes(gg);
            As_across_bins_index{j}=[gg,i];
            id=id+1;  % assembly identity                    
            identity(j)=id;           
            j=j+1;                   
        end
    end
end




% REORDER LAGS AND SHIFT ASSEMBLY'S OCCURRENCE
[As_across_bins,As_across_bins_index]=restyle_assembly_lags_time(As_across_bins,As_across_bins_index);



end



