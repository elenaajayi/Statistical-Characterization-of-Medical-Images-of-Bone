function [patches,nimg_dir_used,read_errs,imrange]=nimg_read_patches(m,patch_list,nimg_dir)
% [patches,nimg_dir_used,read_errs,imrange]=nimg_read_patches(m,patch_list,nimg_dir) reads one or more natural-image files
%   and parses them into patches
%
%  Also checks that files exist, and if not, gets a new directory
%
%  see C:\Users\jdvicto\Dropbox\Textures\RawNIstats\ImagePatchStats\code\extractingPatches\README.rtfd
%  for cutting patches from within images
%
% m: metadata structure returned by nimg_readmetadata
% patch_list: list of patch numbers, from 1 to m.npatches, typically supplied by nimg_select_patches
% nimg_dir: patch directory, if empty or not given, efaults to defaults to 'C:\Users\jdvicto\Desktop\PennNaturalImageDatabase\';
%
% patches: a stack of patch data, size is prod(m.NR_orig) x prod(m.NR_orig) x length(patch_list), NaN if file found
% nimg_dir_used: directory used
%     PID_num (listed in PID_nums)
%     view_num (listed in view_nums, always 1)
%     examp_num (listed in examp_nums)
%     gal_num (listed in gal_nums)
% read_errs: number of read errors (missing files)
% imrange: min and max of raw image values
%
% %See also:  FFDM_BTC_CALC_GEN, NIMG_READMETADATA, NIPATCHES_GETGALS, NIMG_SELECT_PATCHES.
%
if (nargin<=2)
    nimg_dir=[];
end
%
NR=prod(m.NR_orig);
patches=NaN(NR,NR,length(patch_list));
if isempty(nimg_dir)
    nimg_dir='C:\Users\jdvicto\Desktop\PennNaturalImageDatabase\';
end
%check to see if first and last files are present
ifok=0;
while ifok==0
    ifok=1;
    for ifl=1:2
        if (ifl==1)
            patch_no=patch_list(1);
        else
            patch_no=patch_list(2);
        end
        full_name=cat(2,nimg_dir,m.image_file{patch_no});
        if_file=exist(full_name,'file');
        if (if_file==0)
            ifok=0;
            disp(sprintf(' %s not found.',full_name));
        end
    end
    if (ifok==0)
        nimg_dir=getinp('path for natural image database, followed by a slash','s',[],strrep(nimg_dir,'\','/'));
        nimg_dir=strrep(strrep(nimg_dir,'/',filesep),'\',filesep);
    end
end
nimg_dir_used=nimg_dir;
%read requested files, avoiding re-reading a file if previous patch was from same file
img_read=0;
read_errs=0;
maxmax=-Inf;
minmin=Inf;
for ipatch=1:length(patch_list)
    patch_no=patch_list(ipatch);
    img_no=m.PID_num(patch_no);
    if (img_no~=img_read)
        full_name=cat(2,nimg_dir,m.image_file{patch_no});
        if_file=exist(full_name,'file');
        if (if_file==0)
            read_errs=read_errs+1;
            disp(sprintf(' %s not found.',full_name));
        else
            image_data=getfield(load(full_name),'LUM_Image');
        end
        img_read=img_no;
    end
    if (if_file~=0)
        patches(:,:,ipatch)=image_data([0:NR-1]+m.XYPos(patch_no,1),[0:NR-1]+m.XYPos(patch_no,2));
    end
end
if (read_errs>0)
    disp(sprintf(' %4.0f read errors.',read_errs));
end
imrange=[min(patches(:)) max(patches(:))];
return
