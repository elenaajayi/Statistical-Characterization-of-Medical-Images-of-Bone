%nibs_bpdata_demo:  create a border patch database
% invokes Suniyya Waraich's software, typically on Berkeley Segmented Database,
%   starting with segmentImagesBSDB.
%
% 29Mar22: added flexibility for PIDlists
% 25May22: fixed bug so that pathOutputDir can be non-default, changed output dir default
% 06Jun22: add all-PID list, add if_objlist
% 08Jun22: minor reorg of dialogs
% 15Jun22: start support for extraction of surrounding patches for whitening
% 20Jun22: add option for save in -v7.3 and nibs_bpdata_prune for large datasets and 
%
%   See also:  SEGMENTIMAGESBSDB, NIBS_BPDATA_SHOW, NIBS_DEFOPTS, NIBS_GET_META.
%
if ~exist('bpdata_filename') bpdata_filename='bpdata_BSDS_03Jun22b'; end %bpdata_BSDB_[18|23|29]Mar22 are earlier versions
if ~exist('pathImagesDir') pathImagesDir='C:\Users\jdvicto\Dropbox\ObjectBorderStats\Images'; end
if ~exist('pathOutputDir') pathOutputDir='nibs-output'; end
if ~exist('quantile_list') quantile_list=[.1 .25 .5 .75 .9]; end
%
if ~exist('PIDlists')
    PIDlists=struct();
    PIDlists.ver16Dec21=...
        [2 9 10 11 13 26 27 28 47 48 54 59 68 69 79 102 112 121 143 147 166 167 186 192 195];
    PIDlists.ver29Mar22=...
        [2 4 9 10 11 13 20 26 27 28 47 48 54 59 67 68 69 79 102 112 121 143 146 147 166 167 186 192 195];
    PIDlists.all=...
        [1:200];
end
if ~exist('PIDmax') PIDmax=200; end
%
bpdata_filename=getinp('bpdata (border-patch database) file name to read or create','s',[],bpdata_filename);
pathOutputDir=getinp('output directory','s',[],pathOutputDir);
bpdata_fullname=cat(2,pathOutputDir,filesep,bpdata_filename,'.mat');
bpdata_fullname=strrep(bpdata_fullname,'.mat.mat','.mat');
if_bpdata_make=~exist(bpdata_fullname,'file');
if_bpdata_make=getinp('1 to create a border-patch database from scratch','d',[0 1],if_bpdata_make);
if ~exist('if_small_poly_warn') if_small_poly_warn=1; end
if ~exist('if_prune') if_prune=1; end
%
%in case this was not put in nibs_bpdata_demo
if ~exist('opts_nibs') opts_nibs=struct; end
opts_nibs=nibs_defopts(opts_nibs);
%
if ~exist('opts_show') opts_show=struct; end
if ~exist('opts_patchpair') opts_patchpair=struct; end
%
opts_show.if_objlist=getinp('1 to show object list for each image','d',[0 1],0);
if_warnings=getinp('1 to show warnings','d',[0 1],1);
%
if (if_bpdata_make)
    if ~exist('patchProperties')
        patchProperties=struct;
    end
    patchProperties=filldefault(patchProperties,'heightWidthRatio',[1 1.5]); %<1: short and wide, >1: tall and thin; can be 1 x n array
    patchProperties.heightWidthRatio=getinp('aspect ratio (height/width), can be vector','f',[.1 10],patchProperties.heightWidthRatio);
    %
    patchProperties=filldefault(patchProperties,'pixelsPerPatch',256);
    patchProperties.pixelsPerPatch=getinp('pixels per patch','d',[64 8192],patchProperties.pixelsPerPatch);
    if ~exist('algoParams')
        algoParams=struct;
    end
    algoParams=filldefault(algoParams,'smoothness',[1 2 3]); %higher b leads to smoother patches
    algoParams.smoothness=getinp('smoothness, can be vector','d',[1 64],algoParams.smoothness);
    if ~exist('imageSet')
        PIDvers=fieldnames(PIDlists);
        for ilist=1:length(PIDvers)
            PIDver=PIDvers{ilist};
            disp(sprintf('base image list option %1.0f (%10s): %4.0f PIDs, %4.0f to %4.0f',...
                ilist,PIDver,length(PIDlists.(PIDver)),min(PIDlists.(PIDver)),max(PIDlists.(PIDver))));
        end
        ilist=getinp('choice','d',[1 length(PIDvers)]);
        imageSet=getinp('final image list (should be a subset of above)','d',[1 PIDmax],PIDlists.(PIDvers{ilist}));
    end
    if_small_poly_warn=getinp('1 to turn off warnings for small polygons','d',[0 1],if_small_poly_warn);
    if if_small_poly_warn
        warning('off','MATLAB:inpolygon:ModelingWorldLower')
    else
        warning('on','MATLAB:inpolygon:ModelingWorldLower')
    end
    if_prune=getinp('1 to prune the database','d',[0 1],if_prune);
    %[bpdata,bpseg,patchPropertiesUsed,algoParamsUsed]=...
    %    segmentImagesBSDB(patchProperties, algoParams, imageSet, pathImagesDir, pathOutputDir);
    %can use simpler call:
    bpdata=segmentImagesBSDB(patchProperties, algoParams, imageSet, pathImagesDir, pathOutputDir);
    if (if_prune)
         bpdata=nibs_bpdata_prune(bpdata);
    end
    %bpseg=bpdata.segmentation
    %patchPropertiesUsed=bpdata.PatchProperties;
    %algoParamsUsed=bpdata.algoParams;
    ifsave=getinp('1 to save file (-1 to use matlab -v7.3 for large files)','d',[-1 1],1);
    if (ifsave>0)
        save(bpdata_fullname,'bpdata');
        disp(sprintf('bpdata saved in %s',bpdata_fullname));
    elseif (ifsave<0)
        save(bpdata_fullname, '-v7.3','bpdata');
        disp(sprintf('bpdata saved in %s',bpdata_fullname));
    end
else
    bpdata=getfield(load(bpdata_fullname),'bpdata');
    disp(sprintf('bpdata loaded from %s',bpdata_fullname));
end
bpdata=filldefault(bpdata,'mask_dict',opts_nibs.mask_dict_bpdata);
%
%survey for each value of patchProperties and algoParams
%
disp(sprintf('summary of %s',bpdata_fullname));
disp(bpdata);
%
n_aspects=length(bpdata.patchProperties.heightWidthRatio);
n_smooths=length(bpdata.algoParams.smoothness);
%
opts_show=filldefault(opts_show,'quantile_list',quantile_list);
[opts_show_used,warnings]=nibs_bpdata_show(bpdata,opts_show);
disp(sprintf('warnings encountered: %4.0f',size(warnings,1)));
if (if_warnings)
    disp(warnings);
end
if_done=0;
opts_patchpair=filldefault(opts_patchpair,'if_log',1);
opts_patchpair=filldefault(opts_patchpair,'pathImagesDir',pathImagesDir);
opts_patchpair=filldefault(opts_patchpair,'if_showimages',1);
opts_patchpair=filldefault(opts_patchpair,'rect_size',128);
while (if_done==0)
    disp('aspect ratios available:')
    disp(bpdata.patchProperties.heightWidthRatio);
    i_aspect=getinp('index of aspect ratio to show (0 to end)','d',[0 length(bpdata.patchProperties.heightWidthRatio)]);
    if (i_aspect==0)
        if_done=1;
    else
        disp('smoothness ratios available:')
        disp(bpdata.algoParams.smoothness)
        i_smoothness=getinp('index of smoothness parameter to show','d',[1 length(bpdata.algoParams.smoothness)]);
        %
        meta_use=nibs_get_meta(bpdata,i_aspect,i_smoothness);
        patch_pair_ID=Inf;
        while (patch_pair_ID>0)
            patch_pair_ID=getinp('patch pair ID to show (0 to end)','d',[0 meta_use.npatches]);
            if (patch_pair_ID>0)
                rect_size=getinp('size of square (or rectangular) region to extract','d',[0 4096],opts_patchpair.rect_size);
                if max(rect_size)>0
                    opts_patchpair.rect_size=rect_size;
                else
                    opts_patchpair.rect_size=[];
                end
                [patchpair,h,opts_patchpair_used]=nibs_get_patchpair(bpdata,i_aspect,i_smoothness,patch_pair_ID,opts_patchpair);
            end
        end %patch_pair_ID
    end %not done
end
