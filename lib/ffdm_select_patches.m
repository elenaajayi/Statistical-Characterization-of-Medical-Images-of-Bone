function [patch_list,sel_used,avail,patch_list_sel]=ffdm_select_patches(m,sel,if_auto)
% [patch_list,sel_used,avail,patch_list_sel]=ffdm_select_patches(m,sel,if_auto) selects a list of
% patches according to patch number, PID_num (patient ID),  View_num (view), and Examp_num
%
% m: structure returned by ffdm_readmetadata
% sel: structure with fields corresponding to default choices for
%     patch_num (listed in patch_nums)
%     PID_num (listed in PID_nums)
%     view_num (listed in view_nums)
%     examp_num (listed in examp_nums)
%     Quality (listed in Quality_nums) -- this is ignored if m does not
%       have a Quality field, which means that 'Junk_In_Patch' was absent
%       from original csv file
%  In the above, empty means ask with no default, 0 means all
% if_auto: 1 to not ask, just uses defaults (or all), defaults to 0
%
% patch_list: list of selected patches
% sel_used: selection choices
% avail: structure of what is available
% patch_list_sel: structure of patch lists that meet each criterion
%
%  There is a lot of awkwardness here because of capitalization and naming conventions for header rows of original csv sheets
%   Names in input and output arguments, sel, sel_used, and avail are systematic (differ by "_nums"); fields in m are irregular.
%   Take care that sel_fields,m_fields,and avail_fields are properly
%   updated and aligned if additional selection criteria are added.
%
% 15Jan20: add option to select based on 'Quality', defaults to selecting Quality=1 (no junk)
%
% %See also:  FILLDEFAULT, FFDM_READ_SUBPATCH, FFDM_SUBPATCH_DEMO, NIMG_SELECT_PATCHES.
%
if (nargin<=1) sel=[]; end
if (nargin<=2) if_auto=0; end
sel_fields={'patch_nums','PID_nums','view_nums','examp_nums','quality_nums'}; %allowed field names for sel
m_fields={'patch_num','PID_num','View_num','Examp_num','Quality'}; %corresponding field names in m
avail_fields={'patch','PID','view','examp','quality'}; %corresponding field names in avail
sel_defaults={0,0,0,0,1}; %corresponding default selections
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
    PID_nums=getinp('patient ID selection (0 for all)','d',[0 max(avail.PID)],sel.PID_nums);
    if all(PID_nums==0)
        PID_nums=0;
    else
        PID_nums=intersect(PID_nums,avail.PID);
    end
    sel.PID_nums=PID_nums;
    %
    for view_ptr=1:length(avail.view)
        v=avail.view(view_ptr);
        disp(sprintf(' view %1.0f->%s',v,m.view_names{v}));
    end
    view_nums=getinp('view selection (0 for all)','d',[0 max(avail.view)],sel.view_nums);
    if all(view_nums==0)
        view_nums=0;
    else
        view_nums=intersect(view_nums,avail.view);
    end
    sel.view_nums=view_nums;
    %
    examp_nums=getinp('example selection (0 for all)','d',[0 max(avail.examp)],sel.examp_nums);
    if all(examp_nums==0)
        examp_nums=0;
    else
        examp_nums=intersect(examp_nums,avail.examp);
    end
    sel.examp_nums=examp_nums;
    %
end
%
%apply selections
%
patch_list=[1:m.npatches];
patch_list_sel.patch_num=avail.patch;
patch_list_sel.PID_num=avail.patch;
patch_list_sel.View_num=avail.patch;
patch_list_sel.Examp_num=avail.patch;
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
patch_list=intersect(intersect(intersect(patch_list_sel.patch_num,patch_list_sel.PID_num),patch_list_sel.View_num),patch_list_sel.Examp_num);
%apply quality selection if present
if isfield(m,'Quality')
    quality_nums=getinp('quality selection (0 for all)','d',[0 max(avail.quality)],sel.quality_nums);
    if all(quality_nums==0)
        quality_nums=0;
    else
        quality_nums=intersect(quality_nums,avail.quality);
    end
    sel.quality_nums=quality_nums;
    patch_list_sel.Quality=avail.patch;
    if sel.quality_nums(1)>0
        quality_sel=ismember(m.Quality,sel.quality_nums);
        patch_list_sel.Quality=find(quality_sel>0)';
    end
    patch_list=intersect(patch_list,patch_list_sel.Quality);
end
%
sel_used=sel;
return
