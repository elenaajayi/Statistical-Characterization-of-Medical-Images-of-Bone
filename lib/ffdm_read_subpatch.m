function [subpatch,meta_patch,preproc_used]=ffdm_read_subpatch(patchID,meta,patch_path,preproc);
% [subpatch,mpatch,preproc_used]=ffdm_read_subpatch(patchID,meta,patch_path,preproc); reads and prcesses a patch
% from a tiff file
%
% patchID: an index into the patch metadata structure meta
% meta: patch metadata structure, returned by ffdm_readmetadata
% patch_path: path to patch files
% preproc: preprocessing, [] if not supplied
%    preproc.XYLen: length (x,y) of patch to sample, if not supplied, defaults meta.XYLen(patchID,:)
%    preproc.XYOff: offset into patch, defaults to [0 0] if not supplied
%    prepoc.downsample: downsampling factor, defaults to 1, must divide XYLen
%    preproc.scale: full grayscale range, typically 0 to 4095 for 12-bit TIFF
% Here, X is first coord, Y is second, considered as viewing a matrix, which is imagesc conventions:
%    X is 1 at top and increases down
%    Y is 1 at left and increases right
%
% subpatch: patch following preprocessing, scaled to [0 1]
% meta_patch: metadata of patch, drawn from meta, same fields (but cell arrays beome strings)
% preproc_used: preprocessed
%
%  See also: FFDM_SELECT_PATCHES, FILLDEFAULT, FFDM_SUBPATCH_DEMO, FFDM_READMETADATA
%
meta_exclude={'headers','view_names','npatches','patch_count_table'}; %metadata fields to exclude from meta_patch
%
if (nargin<=3)
    preproc=[];
end
preproc=filldefault(preproc,'scale',[0 4095]);
preproc=filldefault(preproc,'XYOff',[0 0]);
preproc=filldefault(preproc,'XYLen',meta.XYLen(patchID,:));
preproc=filldefault(preproc,'downsample',1);
preproc_used=preproc;
%
XY=preproc.XYLen;
d=preproc.downsample;
%
if any(mod(XY,d)>0)
    error(sprintf('requested downsampling (%3.0f %3.0f region size, downsampling by %3.0f) requires non-integer number of pixels per patch',XY,d));
end
if any(preproc.XYOff+XY>meta.XYLen(patchID,:)) | any(preproc.XYOff<0)
    disp(preproc)
    error(sprintf('region requested from patch exceeds patch boundaries, patch is of size %3.0f %3.0f',meta.XYLen(patchID,:)));
end
patch_fullname=cat(2,patch_path,filesep,meta.FName{patchID});
if ~exist(patch_fullname,'file')
    error(sprintf('cannot find patch data file %s',patch_fullname));
end
patch=imread(patch_fullname);
region=patch(preproc.XYOff(1)+[1:XY(1)],preproc.XYOff(2)+[1:XY(2)]); %cut out requested region
region_scaled=(double(region)-preproc.scale(1))/diff(preproc.scale); %rescale
%downsample by block averaging
region_scaled=reshape(region_scaled,[d,XY(1)/d,d,XY(2)/d]);
subpatch=reshape(sum(sum(region_scaled,3),1)/d/d,XY/d);
%
%get metadata
meta_patch=struct();
fns=fieldnames(meta);
for ifn=1:length(fns)
    fn=fns{ifn};
    if isempty(strmatch(fn,meta_exclude,'exact'))
        v=meta.(fn);
        if isnumeric(v)
            meta_patch.(fn)=v(patchID,:);
        else
            if iscell(v)
                meta_patch.(fn)=v{patchID,1};
            end
        end
    end
end
meta_patch.FName_short=strrep(strrep(strrep(meta_patch.FName,'Patch','P'),'.tiff',''),'.tif','');
return
