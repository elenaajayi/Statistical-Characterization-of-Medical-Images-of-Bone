function [patch_list,sel_used,avail,patch_list_sel]=nimg_select_patches(m,sel,if_auto)
% [patch_list,sel_used,avail,patch_list_sel]=nimg_select_patches(m,sel,if_auto) selects a list of
% patches according to patch number, PID_num (image number, takes place of patient ID), example number, and gallery number
%
% m: metadata structure returned by nimg_readmetadata
% sel: structure with fields corresponding to default choices for
%     patch_num (listed in patch_nums)
%     PID_num (listed in PID_nums)
%     view_num (listed in view_nums, always 1)
%     examp_num (listed in examp_nums)
%     gal_num (listed in gal_nums)
%  In the above, empty means ask with no default, 0 means all
% if_auto: 1 to not ask, just uses defaults (or all), defaults to 0
%
% patch_list: list of selected patches
% sel_used: selection choices
% avail: structure of what is available
% patch_list_sel: structure of patch lists that meet each criterion
%
%  There is a lot of awkwardness here because of capitalization and naming conventions for header rows of original csv sheets
%  Names in input and output arguments, sel, sel_used, and avail are systematic (differ by "_nums"); fields in m are irregular.
%
%  26Oct21: minor typos fixed
%
% See also:  FFDM_SELECT_PATCHES, NIMG_READMETADATA, NIPATCHES_GETGALS, NIMG_READ_PATCHES, FFDM_BTC_CALC_GEN.
%
if (nargin<=1) sel=[]; end
if (nargin<=2) if_auto=0; end
sel_fields={'patch_nums','PID_nums','view_nums','examp_nums','gal_nums'}; %allowed field names for sel
m_fields={'patch_num','PID_num','View_num','Examp_num','Gal_num'}; %corresponding field names in m
avail_fields={'patch','PID','view','examp','gallery'}; %corresponding field names in avail
sel_defaults={0,0,1,0,0}; %corresponding default selections
%
avail=[];
%merge supplied selections with default selections, and get available options
m.patch_num=[1:m.npatches]'; %make sure there is a field of numbers for each patch
for isel=1:length(sel_fields)
    m_field=m_fields{isel};
    if isfield(m,m_field)
        sel_field=sel_fields{isel};
        sel=filldefault(sel,sel_field,sel_defaults{isel});
        avail.(avail_fields{isel})=unique(m.(m_field)');
    end
end
%     
if (if_auto==0)
    patch_nums=getinp('patch selection (0 for all)','d',[0 max(avail.patch)],sel.patch_nums);
    if all(patch_nums==0)
        patch_nums=0;
    else
        patch_nums=intersect(patch_nums,avail.patch);
    end
    sel.patch_nums=patch_nums;
    %
    PID_nums=getinp('patient ID (image) selection (0 for all)','d',[0 max(avail.PID)],sel.PID_nums);
    if all(PID_nums==0)
        PID_nums=0;
    else
        PID_nums=intersect(PID_nums,avail.PID);
    end
    sel.PID_nums=PID_nums;
    %
    %     for view_ptr=1:length(avail.view)
    %         v=avail.view(view_ptr);
    %         disp(sprintf(' view %1.0f->%s',v,m.view_names{v}));
    %     end
    %     view_nums=getinp('view selection (0 for all)','d',[0 max(avail.view)],sel.view_nums);
    %     if all(view_nums==0)
    %         view_nums=0;
    %     else
    %         view_nums=intersect(view_nums,avail.view);
    %     end
    %     sel.view_nums=view_nums;
    sel_view_nums=1;
    %
    examp_nums=getinp('example [within each image] selection (0 for all)','d',[0 max(avail.examp)],sel.examp_nums);
    if all(examp_nums==0)
        examp_nums=0;
    else
        examp_nums=intersect(examp_nums,avail.examp);
    end
    sel.examp_nums=examp_nums;
    %
    gal_nums=getinp('gallery selection (0 for all)','d',[0 max(avail.gallery)],sel.gal_nums);
    if all(gal_nums==0)
        gal_nums=0;
    else
        gal_nums=intersect(gal_nums,avail.gallery);
    end
    sel.gal_nums=gal_nums;
end
%
%apply selections
%
patch_list=[1:m.npatches];
patch_list_sel.patch_num=avail.patch;
patch_list_sel.PID_num=avail.patch;
patch_list_sel.View_num=avail.patch;
patch_list_sel.Examp_num=avail.patch;
patch_list_sel.Gal_num=avail.patch;
if sel.patch_nums(1)>0
    patch_list_sel.patch_num=sel.patch_nums;
end
if sel.PID_nums(1)>0
    PID_sel=ismember(m.PID_num,sel.PID_nums);
    patch_list_sel.PID_num=find(PID_sel>0)';
end
if sel.view_nums(1)>0
    view_sel=ismember(m.View_num,sel.view_nums);
    patch_list_sel.View_num=find(view_sel>0)';
end
if sel.examp_nums(1)>0
    examp_sel=ismember(m.Examp_num,sel.examp_nums);
    patch_list_sel.Examp_num=find(examp_sel>0)';
end
if sel.gal_nums(1)>0
    gal_sel=ismember(m.Gal_num(m.PID_num),sel.gal_nums); %find the gallery number from the image (PID) number
    patch_list_sel.Gal_num=find(gal_sel>0);
end
patch_list=intersect(intersect(intersect(intersect(patch_list_sel.patch_num,patch_list_sel.PID_num),patch_list_sel.View_num),patch_list_sel.Examp_num),patch_list_sel.Gal_num);
%
sel_used=sel;
return
