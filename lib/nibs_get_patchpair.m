function [patchpair,h,opts_patchpair_used]=nibs_get_patchpair(bpdata,i_aspect,i_smoothness,patch_pair_ID,opts_patchpair);
% [patchpair,h,opts_patchpair_used]=nibs_get_patchpair(bpdata,i_aspect,i_smoothness,patch_pair_ID,opts_patchpair) retrieves
% a figure-ground patch pair from the border patch database and the corresponding image from the image library
%
%  also extracts a rectangular (typically square) region of the image that is centered as much as possible
%    on the fig_near_border and ground_near_border regions, but adjusted so
%    as not to overflow the image, to be used for whitening
%
% bpdata: database of border patches, creted with segmentImagesBSDB.
% i_aspect: aspect ratio index, index into bpdata.patchProperties.heightWidthRatio
% i_smoothness: smoothing, index into bpdata.algoParams.smoothness
% patch_pair_ID: index into patch pair ID
% opts_patchpair: options for display
%  if_log: 1 to log key values to console (defaults to 0)
%  if_showimages: 1 to show image and patches (defaults to 0)
%  if_noread: 1 to prevent reading of image file, so mask_image and rect_image are returned as zeros (defaults to zero)
%     if set to 1, then opts_patchpair.img_size must have image dimensions
%  img_size: size of image, in pixels [height, width]  (not size of patch), must be supplied if if_noread=1
%  pathImagesDir: path to image directory (defaults to [])
%  hfig: handle to existing figure
%  image_file_suffix: image file suffix, defaults to _LUM
%  image_var_name: name of image variable, defaults to 'LUM_Image'
%  image_range: range of image, defaults to [0 256], autoscaled if empty
%  roi_color: color for border of the ROI of the requested patch, defaults to 'r', not drawn if empty
%  roi_color_others: color for border of the other ROIs in the image, defaults to 'y', not drawn if empty
%  roi_color_otherobj: color for border of the other ROIs of other objects in the image, defaults to 'g', not drawn if empty
%  roi_color_rect: color for rect_image, defaults to 'm', not drawn if empty or rect_size absent or empty
%  roi_linewidth:line with for roi borders
%  composite_range: range for composite image, typically [0 6], autoscaled if empty
%  rect_size: size of a rectangular region centered on the combined
%     centroid of fig_near_border and ground_near_border.  If empty,
%     rect_image is not computed.  If scalar, same value used for both sides
%  composite_marker_size: marker size for composite map, defaults to 2
%
% patchpair: structure with the retrieved patches, with added fields:
%   patchpair.mask_image: the image cropped to the regions of the masks
%   patchpair.rect_image: the image cropped to the rectangular region specified by rect_size
%   patchpair.rect_mask: masks, one field for each field of patchpair.mask, resized for rect_image
% h: handle to plot(s) made
% opts_patchpair_used: opts_patchpair, after defaults, and including a field img_size
%
% 02Jun22: plot and tally patch pairs from other objects; nice colorbar with key, change plot order; adjustable line width
% 15Jun22: add mask_image, rect_image, rect_mask; rearrange plot
% 30Jun22: add opts_patchpair.if_noread, opts_patchpair.img_size
%
%   See also:  NIBP_BPDATA_DEMO, SEGMENTIMAGESBSB, NIBS_BPDATA_SHOW, NIBS_DEFOPTS, NIBS_GET_META,
%    NIBS_PATCH2PAIRSTRING, NIBS_DEFOPTS, FFDM_BTC_CALC_GEN, NIBS_GET_MASK.
%
if (nargin<=4)
    opts_patchpair=struct;
end
fg_string={'figure','ground'};
%
nibs_opts=nibs_defopts;
mask_labels=fieldnames(nibs_opts.mask_dict_composite);
mask_tokens=zeros(1,length(mask_labels));
mask_labels_nice=cell(1,length(mask_labels));
for ilabel=1:length(mask_labels)
    mask_tokens(ilabel)=nibs_opts.mask_dict_composite.(mask_labels{ilabel}); %known map tokens
    mask_labels_nice{ilabel}=strrep(mask_labels{ilabel},'_',' ');
end
[mask_tokens,mask_idx]=sort(mask_tokens);
mask_labels=mask_labels(mask_idx);
mask_labels_nice=mask_labels_nice(mask_idx);
%
opts_patchpair=filldefault(opts_patchpair,'if_log',0);
opts_patchpair=filldefault(opts_patchpair,'if_showimages',0);
opts_patchpair=filldefault(opts_patchpair,'if_noread',0);
opts_patchpair=filldefault(opts_patchpair,'img_size',[0 0]);
opts_patchpair=filldefault(opts_patchpair,'pathImagesDir',[]);
opts_patchpair=filldefault(opts_patchpair,'hfig',[]);
opts_patchpair=filldefault(opts_patchpair,'image_file_suffix','_LUM');
opts_patchpair=filldefault(opts_patchpair,'image_var_name','LUM_Image');
opts_patchpair=filldefault(opts_patchpair,'roi_color','r');
opts_patchpair=filldefault(opts_patchpair,'roi_color_others','y');
opts_patchpair=filldefault(opts_patchpair,'roi_color_otherobj','g');
opts_patchpair=filldefault(opts_patchpair,'roi_color_rect','m');
opts_patchpair=filldefault(opts_patchpair,'roi_linewidth',2);
opts_patchpair=filldefault(opts_patchpair,'image_range',[0 256]);
opts_patchpair=filldefault(opts_patchpair,'composite_range',[min(mask_tokens) max(mask_tokens)]);
opts_patchpair=filldefault(opts_patchpair,'composite_marker_size',2);
%
if ~isfield(opts_patchpair,'rect_size')
    opts_patchpair.rect_size=[];
end
if length(opts_patchpair.rect_size)==1
    opts_patchpair.rect_size=repmat(opts_patchpair.rect_size,1,2);
end
%
h=struct;
patchpair=struct;
%
meta_use=nibs_get_meta(bpdata,i_aspect,i_smoothness);
%
PID_num=meta_use.PID_num(patch_pair_ID);
Examp_num=meta_use.Examp_num(patch_pair_ID);
image_file=meta_use.image_file{patch_pair_ID};
%
patchpair.patch_pair_ID=patch_pair_ID;
patchpair.PID_num=PID_num;
patchpair.Examp_num=Examp_num;
patchpair.image_file=image_file;
%
seg=bpdata.segmentation{PID_num}.patches{i_aspect,i_smoothness};
patchpair.XYPos=seg.pairedPatches{Examp_num}.XYPos;
patchpair.XYLen=seg.pairedPatches{Examp_num}.XYLen;
patchpair.mask=seg.pairedPatches{Examp_num}.mask;
%
for ifg=1:2
    ifg_no=2*(Examp_num-1)+ifg;
    patchpair.order{1,ifg}=fg_string{ifg};
    patchpair.centroids(ifg,:)=seg.patch{ifg_no}.centroid;
    patchpair.areas(ifg,:)=seg.patch{ifg_no}.area;
    patchpair.neighbor_patches{ifg}=seg.patch{ifg_no}.neighbors;
end
PID_patchpairs=find(meta_use.PID_num==PID_num);
tstring_count=sprintf('pair %3.0f of %3.0f pairs in image %s',...
    patch_pair_ID-min(PID_patchpairs)+1,length(PID_patchpairs),image_file);
tstring_obj=[];
%if object type fields are available, then make a descriptor and determine
%what patch pairs are in other objects of this image
obj_id=[];
Examps_otherobj=[];
if isfield(meta_use,'obj_ids') & isfield(meta_use,'obj_type')
    obj_id=meta_use.obj_ids{patch_pair_ID};
    PID_patchpairs_thisobj=strmatch(obj_id,meta_use.obj_ids,'exact');
    PID_patchpairs_otherobj=setdiff(PID_patchpairs,PID_patchpairs_thisobj);
    Examps_otherobj=meta_use.Examp_num(PID_patchpairs_otherobj);
    tstring_obj=sprintf(', obj %7s (%3.0f of %3.0f pairs)',...
        obj_id,min(find(PID_patchpairs_thisobj==patch_pair_ID)),length(PID_patchpairs_thisobj));
end
patchpair.obj_id=obj_id;
tstring=sprintf('patch pair ID %4.0f:  PID %4.0f example %4.0f',patch_pair_ID,PID_num,Examp_num);
tstring_params=sprintf('aspect ratio %1.0f->%5.3f smoothness %1.0f->%4.2f',...
    i_aspect,bpdata.patchProperties.heightWidthRatio(i_aspect),...
    i_smoothness,bpdata.algoParams.smoothness(i_smoothness));
if (opts_patchpair.if_log)
    disp(sprintf(' %s XYPos %5.0f %5.0f XYLen %5.0f %5.0f file name %s',tstring,patchpair.XYPos,patchpair.XYLen,image_file));
    disp(cat(2,tstring_count,' ',tstring_obj));
    for ifg=1:2
        neighbor_strings=nibs_patch2pairstring(patchpair.neighbor_patches{ifg});
        disp(sprintf('   %8s: centroid %7.2f %7.2f   area %6.0f  neighbor patches %s  from pairs %s',...
            patchpair.order{ifg},patchpair.centroids(ifg,:),patchpair.areas(ifg,:),...
            sprintf('%5.0f',patchpair.neighbor_patches{ifg}(:)),...
            sprintf('%s    ',neighbor_strings{:})));
    end
end
%
%read the image
%
%uniformize file name
fullname=cat(2,opts_patchpair.pathImagesDir,filesep,image_file,opts_patchpair.image_file_suffix);
fullname=strrep(fullname,'/',filesep);
fullname=strrep(fullname,'\',filesep);
fullname=strrep(fullname,cat(2,filesep,filesep),filesep);
fullname=cat(2,fullname,'.mat');
fullname=strrep(fullname,'.mat.mat','.mat');
opts_patchpair.fullname=fullname;
%
if (opts_patchpair.if_noread==0)
    img=getfield(load(fullname),opts_patchpair.image_var_name);
    opts_patchpair.img_size=size(img);
else
    img=zeros(opts_patchpair.img_size); %do not attempt to read
end
patchpair.mask_image=img(patchpair.XYPos(1)+[0:patchpair.XYLen(2)],patchpair.XYPos(2)+[0:patchpair.XYLen(1)]);
patchpair.rect_image=[];
patchpair.rect_XYPos=[];
patchpair.rect_mask=struct;
if ~isempty(opts_patchpair.rect_size)
    if any(size(img)<opts_patchpair.rect_size)
        warning(sprintf('requested rectangle size (%4.0f %4.0f) exceeds image size (%4.0f x %4.0f), so no rectangle extracted',...
            opts_patchpair.rect_size,size(img)));
    else
        %disp(opts_patchpair.rect_size)
        %establish centroids relative to patchpair
        fnb_centroid=nibs_get_patchpair_centroid(patchpair.mask.fig_near_border);
        gnb_centroid=nibs_get_patchpair_centroid(patchpair.mask.ground_near_border);
        rect_centroid=(fnb_centroid+gnb_centroid)/2;
        rect_XYPos=round(patchpair.XYPos+rect_centroid-opts_patchpair.rect_size/2); %without regard to image boundary
        rect_XYPos=rect_XYPos+max(0,1-rect_XYPos); %adjust for image boundary on left and top
        rect_XYPos=rect_XYPos-max(0,rect_XYPos+opts_patchpair.rect_size-size(img)-1); %adjust for image boundary on right and bottom
        patchpair.rect_image=img(rect_XYPos(1)+[0:opts_patchpair.rect_size(1)-1],rect_XYPos(2)+[0:opts_patchpair.rect_size(2)-1]);
        patchpair.rect_XYPos=rect_XYPos;
        %make masks for rectangular region adjusting source and detsination based on 
        % positions (rect_XYPos, patchpair.XYPos) and sizes opts_patchpair.rect_size, patchpair.XYLen
        %XYsrc_range=[ones(2,1) flipud(patchpair.XYLen(:))]; %flip because XYLen has coords reversed
        %XYdst_range=[ones(2,1) opts_patchpair.rect_size(:)];
        XYsrc_range=zeros(2,2); %first row: X, second row: Y; columns: beg and end
        XYdst_range=zeros(2,2); %first row: X, second row: Y; columns: beg and end
        %
        offset_beg=rect_XYPos(:)-patchpair.XYPos(:);
        XYsrc_range(:,1)=1+max(+offset_beg,0);
        XYdst_range(:,1)=1+max(-offset_beg,0);
        %
        offset_end=rect_XYPos(:)+opts_patchpair.rect_size(:)-(patchpair.XYPos(:)+flipud(patchpair.XYLen(:))); %flip because XYLen has coords reversed
        XYsrc_range(:,2)=flipud(patchpair.XYLen(:))-max(-offset_end,0);
        XYdst_range(:,2)=opts_patchpair.rect_size(:)-max(+offset_end,0);
        %
        mask_list=fieldnames(patchpair.mask);
        for im=1:length(mask_list)
            fn=mask_list{im};
            mask_orig=patchpair.mask.(fn);
            rect_mask=zeros(opts_patchpair.rect_size,'uint8');
            rect_mask([XYdst_range(1,1):XYdst_range(1,2)],[XYdst_range(2,1):XYdst_range(2,2)])=...
                mask_orig([XYsrc_range(1,1):XYsrc_range(1,2)],[XYsrc_range(2,1):XYsrc_range(2,2)]);
            patchpair.rect_mask.(fn)=rect_mask;
        end
        patchpair.rect_mask.XYsrc_range=XYsrc_range;
        patchpair.rect_mask.XYdst_range=XYdst_range;
    end
end
%
if (opts_patchpair.if_showimages) %5.0f file name %s',...
    if (~isempty(opts_patchpair.hfig))
        hfig=opts_patchpair.hfig;
    else
        hfig=figure;
        opts_patchpair.hfig=hfig;
        set(gcf,'Position',[100 100 1200 800]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
    end
    %
    if isempty(opts_patchpair.composite_range)
        comp_lohi=[0 1]; %should never happen
    else
        comp_lohi=opts_patchpair.composite_range;
    end
    cmap_gray=colormap(gray);
    cmap_jet=colormap(jet); 
    cmap_step=cmap_jet; %modify cmap_jet to be stepwise, to only have colors actually used for tokens
    njet=size(cmap_jet,1);
    nsteps=max(comp_lohi)-min(comp_lohi);
    for k=1:njet
        nearest=floor((k-1)/njet*(nsteps+1));
        cmap_step(k,:)=cmap_jet(floor((nearest/nsteps)*(njet-1))+1,:);
    end
    ncmap=2*floor(size(cmap_gray,1)/2);
    cmap_use=zeros(ncmap,3);
    cmap_use(1:ncmap/2,:)=cmap_gray([2:2:ncmap],:); %use bottom half of the color map for gray
    cmap_use(ncmap/2+[1:ncmap/2],:)=cmap_step([2:2:ncmap],:); %use upper half of colormap for other colors
    set(gcf,'Colormap',cmap_use);
    %
    if ~isempty(opts_patchpair.image_range)
        im_lohi=opts_patchpair.image_range;
    else
        im_lohi=[min(img(:)) max(img(:))];
    end
    if diff(im_lohi)~=0
        im_scaled=min(1-2/ncmap,max(0,(img-im_lohi(1))/diff(im_lohi)));
        im_rect_scaled=min(1-2/ncmap,max(0,(patchpair.rect_image-im_lohi(1))/diff(im_lohi)));
    else
        im_scaled=zeros(size(img));
        im_rect_scaled=zeros(size(patchpair.rect_image));
    end
    %
    %show entire image
    %
    subplot(2,3,1);
    imagesc(im_scaled,[0 2]); %plot in first half of colormap
    hold on;
    nExamps=length(seg.pairedPatches);
    %draw a rectangle on all the other patchpairs in other objects
    if ~isempty(opts_patchpair.roi_color_otherobj)
        nibs_get_patchpair_showroi(patchpair,seg,Examps_otherobj,opts_patchpair.roi_color_otherobj,opts_patchpair.roi_linewidth);
    end %other rectangles
    %draw a rectangle on all the other patchpairs in same object
    if ~isempty(opts_patchpair.roi_color_others)
        nibs_get_patchpair_showroi(patchpair,seg,setdiff([1:nExamps],union(Examp_num,Examps_otherobj)),opts_patchpair.roi_color_others,opts_patchpair.roi_linewidth);
    end %other ROIs in this object
    %draw a rectangle around the surrounding rectangle
    if ~isempty(opts_patchpair.roi_color_rect) & ~isempty(patchpair.rect_image) & ~isempty(patchpair.rect_XYPos);
        hp=plot(...
            patchpair.rect_XYPos(2)+[0 1 1 0 0]*size(patchpair.rect_image,2),...
            patchpair.rect_XYPos(1)+[0 0 1 1 0]*size(patchpair.rect_image,1),'k');
        set(hp,'Color',opts_patchpair.roi_color_rect);
        set(hp,'LineWidth',opts_patchpair.roi_linewidth);
    end
    %draw a rectangle corresponding to the requested patchpair
    if ~isempty(opts_patchpair.roi_color)
        nibs_get_patchpair_showroi(patchpair,seg,Examp_num,opts_patchpair.roi_color,opts_patchpair.roi_linewidth);
    end
    axis equal;
    axis tight;
    title(image_file,'Interpreter','none');
    %
    %show the ROI that contains the patch
    %
    subplot(2,3,2);
    imagesc(im_scaled(patchpair.XYPos(1)+[0:patchpair.XYLen(2)],patchpair.XYPos(2)+[0:patchpair.XYLen(1)]),[0 2]);
    axis equal;
    axis tight;
    title(image_file,'Interpreter','none');
    %
    %show the rectangular region that contains the patch, with dots
    %corresponding to figure and ground
    %
    if ~isempty(patchpair.rect_image) & ~isempty(patchpair.rect_XYPos);
        rect_composite=double(patchpair.rect_mask.composite);
        rect_composite_scaled=min(1,max(0,(rect_composite-comp_lohi(1))/diff(comp_lohi)));
        subplot(2,3,4);
        imagesc(im_rect_scaled,[0 2]);
        axis equal;
        axis tight;
        hold on;
        nibs_get_patchpair_compmap(rect_composite_scaled,ncmap,cmap_use,opts_patchpair);
        title(cat(2,image_file,sprintf(' %3.0f x %3.0f rect',size(im_rect_scaled))),'Interpreter','none');
    end  
    %
    %show composite map
    %
    subplot(2,3,3);
    composite=double(patchpair.mask.composite);
    composite_scaled=min(1,max(0,(composite-comp_lohi(1))/diff(comp_lohi)));
    imagesc(1+composite_scaled,[0 2]); %map composite range into [1 2]
    axis equal;
    axis tight;
    %
    %show the rectangular region that contains the patch, with dots
    %corresponding to figure and ground
    %
    subplot(2,3,5);
    imagesc(im_scaled(patchpair.XYPos(1)+[0:patchpair.XYLen(2)],patchpair.XYPos(2)+[0:patchpair.XYLen(1)]),[0 2]);
    axis equal;
    axis tight;
    hold on;
    nibs_get_patchpair_compmap(composite_scaled,ncmap,cmap_use,opts_patchpair);
    title(cat(2,image_file,' ROI'),'Interpreter','none');
    %
    %make a key
    %
    subplot(2,3,6);
    axis off;
    hcolor=colorbar;
    set(hcolor,'Limits',[.5 1]);
    set(hcolor,'Ticks',0.5+([0:diff(comp_lohi)]+0.5)/(2*diff(comp_lohi)+2));
    set(hcolor,'TickLabels',mask_labels_nice);
    set(hcolor,'Position',[.78 .12 .05 .3]);
    hold on;
    %
    %add global legends
    %
    axes('Position',[0.02,0.08,0.01,0.01]);
    text(0,0,tstring,'Interpreter','none');
    axis off;
    axes('Position',[0.02,0.05,0.01,0.01]);
    text(0,0,cat(2,tstring_count,' ',tstring_obj),'Interpreter','none');
    axis off;
    axes('Position',[0.02,0.02,0.01,0.01]);
    text(0,0,tstring_params,'Interpreter','none');
    axis off;
end
%
opts_patchpair_used=opts_patchpair;
return

function nibs_get_patchpair_showroi(patchpair,seg,examp_list,roi_color,roi_linewidth)
%utility to draw a box around an ROI
for ix=1:length(examp_list)
    iexamp=examp_list(ix);
    XYPos=seg.pairedPatches{iexamp}.XYPos;
    XYLen=seg.pairedPatches{iexamp}.XYLen;
%note that image is plotted with first coordinate vertical -- applies to XYPos but not XYLen
    hp=plot(XYPos(2)+[0 1 1 0 0]*XYLen(1),XYPos(1)+[0 0 1 1 0]*XYLen(2),'k');
    hold on;
    set(hp,'Color',roi_color);
    set(hp,'LineWidth',roi_linewidth);
end
return

function c=nibs_get_patchpair_centroid(umap)
%utility to find a centroid of a uint8 map
c=nan(1,2);
if all(umap(:)==0)
    return
else
    dmap=double(umap);
    %
    colsum=sum(dmap,1);
    colidx=[1:size(dmap,2)];
    %
    rowsum=sum(dmap,2)';
    rowidx=[1:size(dmap,1)];
    c=[sum(rowsum.*rowidx),sum(colsum.*colidx)]/sum(dmap(:));
end
return

function nibs_get_patchpair_compmap(comp_scaled,ncmap,cmap_use,opts_patchpair)
%utility to plot a composite map
ic_uniques=unique(comp_scaled);
for ic_ptr=1:length(ic_uniques)
    ic=ic_uniques(ic_ptr);
    colormap_val=cmap_use(max(ncmap/2+1,ceil(ncmap*(1+ic)/2)),:);
    [ix,iy]=find(comp_scaled'==ic);
    hp=plot(ix(:),iy(:),'k.');
    set(hp,'color',colormap_val);
    set(hp,'MarkerSize',opts_patchpair.composite_marker_size);
end
return

