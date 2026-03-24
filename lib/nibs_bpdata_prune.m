function bpdata_new=nibs_bpdata_prune(bpdata,keepfields,keep_aspect,keep_smoothness)
% function bpdata_new=nibs_bpdata_prune(bpdata,keepfields) is a utility that removes
% fields fro bpdata_segmentation that are not needed to recover patch pairs
%
% in nibs_get_patchpair, only ref to segmentation is
%  seg=bpdata.segmentation{PID_num}.patches{keep_aspect,keep_smoothness};
%
% input:
%    bpdata: full boundary patch database from segmentimagesbsdb
%    keepfields: list of fields to keep (defaults to {'patches'})
%    keep_aspect: if provided, only this value of keep_aspect are kept; []: keep all
%    keep_smoothness: if provided, only this value of keep_smoothness are kept; []: keep all
%
% output:
%   bpdata_new: bpdata, with bpdata.segmentation pruned to keep only fields listed as keepfields
%
%   See also:  NIBP_BPDATA_DEMO,  SEGMENTIMAGESBSDB. NIBS_GET_PATCHPAIR, NIBS_DEFOPTS.
%
if nargin<=1
    keepfields={'patches'};
end
if (nargin<=2)
    keep_aspect=[];
end
if isempty(keep_aspect)
    keep_aspect=[1:length(bpdata.patchProperties.heightWidthRatio)];
end
if (nargin<=3)
    keep_smoothness=[];
end
if isempty(keep_smoothness)
    keep_smoothness=[1:length(bpdata.algoParams.smoothness)];
end
bpdata_new=rmfield(bpdata,'segmentation');
%
bpdata_new.segmentation=cell(1,bpdata.nsubjs);
for PID_num=1:bpdata.nsubjs
    sold=bpdata.segmentation{PID_num};
    if isstruct(sold)
        snew=struct;
        for ik=1:length(keepfields)
            if isfield(sold,keepfields{ik})
                stemp=sold.(keepfields{ik});
                if (size(stemp,1)==length(bpdata.patchProperties.heightWidthRatio)) & (size(stemp,2)==length(bpdata.algoParams.smoothness))
                    stemp2=cell(size(stemp));
                    for i_aspect=1:size(stemp,1)
                        for i_smoothness=1:size(stemp,2)
                            if ismember(i_aspect,keep_aspect) & ismember(i_smoothness,keep_smoothness)
                                stemp2{i_aspect,i_smoothness}=stemp{i_aspect,i_smoothness};
                            end
                        end % i_smoothness
                    end %i_aspect
                    stemp=stemp2;
                end %proper size
                snew.(keepfields{ik})=stemp;
            end
        end
        bpdata_new.segmentation{PID_num}=snew;
    else
        bpdata_new.segmentation{PID_num}=sold;
    end
end
return
end
