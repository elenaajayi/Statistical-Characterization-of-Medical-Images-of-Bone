function [patch_string,pair_vecs]=nibs_patch2pairstring(patch_num)
% [patch_string,pair_vecs]=nibs_patch2pairstring(patch_num) is a utility that converts one
% or more patch numbers to patch pair strings
%
%   patch_num: patch number, or vector of patch numbers
%
%   patch_string: string for display of pair number and figure/ground tag
%   pair_vecs: 2-column array, first column is pair number, second is 1 or 2 (indicating fig or ground)
%
fg_tag={'fig','gnd'};
patch_string=cell(1,length(patch_num));
pair_vecs=zeros(length(patch_num),2);
for ip=1:length(patch_num)
    ifg=1+mod(patch_num(ip)-1,2);
    pair_vecs(ip,:)=[1+(patch_num(ip)-ifg)/2,ifg];
    patch_string{ip}=sprintf('%s%4.0f',fg_tag{ifg},pair_vecs(ip,1));
end
return


