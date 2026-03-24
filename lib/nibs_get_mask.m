function [mask_composite,mask_fullsize,mask_sampscaled]=nibs_get_mask(bpdata,i_aspect,i_smoothness,patch_pair_ID,img_size,rect_size,regtype_masks,opts_ds)
% [mask_composite,mask_fullsize,mask_sampscaled]=nibs_get_mask(bpdata,i_aspect,i_smoothness,patch_pair_ID,img_size,rect_size,regtype_masks,opts_ds)
% retrieves a mask from the boundary-parsed metadata, combines according to
% specified region types (figure near border, ground, etc), and optionally  selects and downsamples
%
% input:
% bpdata: database of border patches, creted with segmentImagesBSDB.
% i_aspect: aspect ratio index, index into bpdata.patchProperties.heightWidthRatio
% i_smoothness: smoothing, index into bpdata.algoParams.smoothness
% patch_pair_ID: index into patch pair ID
% img_size: original image size (height, width), as returned by nibs_get_patchpair
% rect_size: rectangle size (height, width) for region surrounding the patch
% regtype_masks: a list of integer tags for the components of the composite map to be used
%    if empty, then mask_fullsize is not computed
% opts_ds: downsampling options -- used to crete mask_sampled from mask_fullsize
%    opts_ds.xs: values of first  coordinate to select (in ffdm_btc_calc_gen, Rsample{ixs}), defaults to [1:length of dim1 of mask]=[1:rect_size(1)]
%    opts_ds.ys: values of second coordinate to select (in ffdm_btc_calc_gen, Rsample{iys}), defaults to [2:length of dim1 of mask]=[1:rect_size(2)]
%    opts_ds: blocking (in ffdm_btc_calc_gen, S), defaults to 1
%
% output:
% mask_composite: a mask of integers, corresponding to the different region types, of size rect_size
% mask_fullsize: (only computed if regtype_masks is non-empty) -- 0 where mask selects image, 1 where it de-selects,
%    same size as mask_composite
% mask_sampscaled: sampled and scaled mask_fullsize (only computed if opts_ds is non-empty)
%    size is [length(opts_ds.xs)/opts_ds.S),[length(opts_ds.ys)/opts_ds.S)
%
%
%   See also:  NIBP_BPDATA_DEMO, SEGMENTIMAGESBSB, NIBS_GET_PATCHPAIR, FFDM_BTC_CALC_GEN, FILLDEFAULT.
%
if (nargin<=6)
    regtype_masks=[];
end
if (nargin<=7)
    opts_ds=[];
end
mask_fullsize=[];
mask_sampscaled=[];
if length(rect_size)==1
    rect_size=[rect_size,rect_size];
end
%
opts_patchpair=struct;
opts_patchpair.if_noread=1;
opts_patchpair.rect_size=rect_size;
opts_patchpair.img_size=img_size;
%retrieve the composite match from the patch pair
mask_composite=getfield(getfield(nibs_get_patchpair(bpdata,i_aspect,i_smoothness,patch_pair_ID,opts_patchpair),....
    'rect_mask'),'composite');
%
if ~isempty(regtype_masks) %create a mask with 0 in selected pixels, 1 in de-selected pixels
    mask_fullsize=ones(size(mask_composite)); %assume nothing is to be used
    for im=1:length(regtype_masks) %set to zero any pixel to use 
        mask_fullsize(mask_composite==regtype_masks(im))=0;
    end
end
%
if ~isempty(opts_ds)
    opts_ds=filldefault(opts_ds,'xs',[1:rect_size(1)]);
    opts_ds=filldefault(opts_ds,'ys',[1:rect_size(2)]);
    opts_ds=filldefault(opts_ds,'S',1);
    S=opts_ds.S;
    mask_sampled=mask_fullsize(opts_ds.xs,opts_ds.ys);
    mask_sampscaled=squeeze(mean(mean(reshape(mask_sampled,[S,size(mask_sampled,1)/S,S,size(mask_sampled,2)/S]),3),1));
    mask_sampscaled=(mask_sampscaled~=0); %any nonzero value goes to 1
end
return


