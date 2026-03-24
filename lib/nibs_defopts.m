function opts_use=nibs_defopts(opts)
%opts_use=nibs_defopts(opts) sets default options for natural image boundary statistics modules
%
% opts: optional overrides, may be omitted
%
% opts_use: full options structure
%
% 01Jun22: add other object
% 13Jun22: add obj_letter, obj_desc
%
%   See also:  NIBS_BPDATA_DEMO, NIBS_BPDATA_SHOW
%
if (nargin==0)
    opts=struct;
end
%mask_dict to use if omitted from bpdata
mask_dict_bpdata=struct;
mask_dict_bpdata.unassigned=0;
mask_dict_bpdata.in_figure_far_from_border=1;
mask_dict_bpdata.in_figure_close_to_border=2;
mask_dict_bpdata.in_ground_close_to_border=3;
mask_dict_bpdata.in_ground_far_from_border=4;
mask_dict_bpdata.other_object_near_border=5;
mask_dict_bpdata.other_object=6;
opts=filldefault(opts,'mask_dict_bpdata',mask_dict_bpdata);
%
%mask_dict for parsing composite mask
mask_dict_composite=struct;
mask_dict_composite.unassigned=0;
mask_dict_composite.figure_only=1;
mask_dict_composite.fig_near_border=2;
mask_dict_composite.ground_near_border=3;
mask_dict_composite.ground_only=4;
mask_dict_composite.other_object_near_border=5;
mask_dict_composite.other_object=6;
opts=filldefault(opts,'mask_dict_composite',mask_dict_composite);
%
% below added 13Jun22
opts.obj_letter={'a','m'};
opts.obj_desc={'animal','man-made'};
%
opts_use=opts;
return

