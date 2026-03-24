function meta_use=nibs_get_meta(bpdata,i_aspect,i_smooth)
% meta_use=nibs_get_meta(bpdata,i_aspect,i_smooth) gets a metadata structure from the border patch database
% and copies some fields from bpdata
%  works with original format (shared fields) or format after March 23, of individual fields
%
% bpdata: a data structure created by nibp_bdata_demo, calling segmentImagesBSDB.
% i_aspect: pointer to value in bpdata.patchProperties.heightWidthRatio
% i_smooth: pointer to value in bpdata.algoParams.smoothness
%
% meta_use: metadata structure, e.g., 
%              npatches: 682
%                   PID: {682×1 cell}
%               PID_num: [682×1 double]
%            image_file: {682×1 cell}
%              View_num: [682×1 double]
%                  View: {682×1 cell}
%                 XYPos: [682×2 double]
%                 XYLen: [682×2 double]
%     patch_count_table: [200×1 double]
%             Examp_num: [682×1 double]
%        nimg_meta_path: 'nibs-output'
%        nimg_meta_file: 'BorderPatchesBSDS.mat'
%               NR_orig: []
%                nsubjs: 200
%                 views: {'std'}
%      gal_names_sorted: {'BSDS'}
%          gal_img_list: {[1×200 double]}
%               Gal_num: [1×200 double]
%
%   See also: NIBS_BPDATA_DEMO, NIBS_BPDATA_SHOW, FFDM_BTC_CALC_GEN.
%
%subfields of bpdata shared by all values of  heightWidthRatio and  smoothness
shared_fields={'nimg_meta_path','nimg_meta_file','NR_orig','nsubjs','views','gal_names_sorted','gal_img_list','Gal_num'};
%
if isfield(bpdata,'patch_meta') %separate metadata stored for each aspect ratio and smoothness
    meta_use=bpdata.patch_meta{i_aspect,i_smooth};
    for ishared=1:length(shared_fields)
        meta_use.(shared_fields{ishared})=bpdata.(shared_fields{ishared});
    end
else %legacy: only one set of metadata
    meta_use=rmfield(bpdata,{'patchProperties','algoParams','segmentation'});
end
return

