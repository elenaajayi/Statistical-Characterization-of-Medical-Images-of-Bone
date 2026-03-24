function [opts_show_used,warnings]=nibs_bpdata_show(bpdata,opts_show,opts_nibs)
% [opts_show_used,warnings]=nibs_bpdata_show(bpdata,opts_show,opts_nibs) displays the contents of bpdata
%  and does some sanity checks
%
% bpdata: a data structure created by nibp_bdata_demo, calling segmentImagesBSDB.
% opts_show: options structure, may be omitted
%   opts_show.quantile_list: list of qusntiles of patch sizes to show
%   aspect_ptrs: pointers to aspect ratios to show
%   smoothness_ptrs: pointers to smoothness ratios to show
%   if_objlist: 1 to show object list for each image
% opts_nibs: global options for nibs modules, set by nibs_defopts if not supplied
%
%  opts_show_used: values used
%  warnings: warnings
%
% 25May22: changed access to patchMap so that it could also be a cell.
% 01Jun22: verify that region type exists, to enable backward compatibility
% 01Jun22: handle images with multiple objects and multiple patch maps
% 02Jun22: add sanity checks for multiple objects
% 06Jun22: documentation fix, add if_objlist
% 08Jun22: add tally of animal vs man-made objects
% 13Jun22: begin more details on object types
%
%   See also:  NIBS_BPDATA_DEMO, NIBS_DEFOPTS, NIBS_GET_META, NIBS_BPDATA_PARSEOBJ, FFDM_BTC_CALC_GEN.
%
quantile_list=[.1 .25 .5 .75 .9];
if (nargin<=1)
    opts_show=struct;
end
if (nargin<=2)
    opts_nibs=struct;
end
opts_nibs=nibs_defopts(opts_nibs);
nobj_types=length(opts_nibs.obj_letter);
%
opts_show=filldefault(opts_show,'quantile_list',quantile_list);
opts_show=filldefault(opts_show,'aspect_ptrs',[1:length(bpdata.patchProperties.heightWidthRatio)]);
opts_show=filldefault(opts_show,'smoothness_ptrs',[1:length(bpdata.algoParams.smoothness)]);
opts_show=filldefault(opts_show,'if_objlist',0);
opts_show_used=opts_show;
quantile_list=opts_show.quantile_list;
warnings=[];
%
for aspect_ptr=1:length(opts_show.aspect_ptrs)
    i_aspect=opts_show.aspect_ptrs(aspect_ptr);
    aspect=bpdata.patchProperties.heightWidthRatio(i_aspect);
    for smoothness_ptr=1:length(opts_show.smoothness_ptrs)
        i_smooth=opts_show.smoothness_ptrs(smoothness_ptr);
        smooth=bpdata.algoParams.smoothness(i_smooth);
        %
        disp(' ');
        %
        as_string=sprintf('aspect ratio %4.2f smoothness %4.2f',aspect,smooth);
        disp(sprintf(' survey for %s',as_string));
        %
        meta_use=nibs_get_meta(bpdata,i_aspect,i_smooth);
        PID_unique=unique(meta_use.PID_num);
        nPIDs_used=length(PID_unique);
        %
        disp(sprintf('total number of patch pairs: %5.0f',meta_use.npatches));
        %do statistics of object types of obj_ids is available
        if isfield(meta_use,'obj_ids')
            unique_objs=unique(meta_use.obj_ids);
            obj_table=nibs_bpdata_parseobj(meta_use.obj_ids);
            for iobj_type=1:nobj_types
                disp(sprintf('   %24s: %5.0f pairs',...
                    opts_nibs.obj_desc{iobj_type},sum(obj_table.type==opts_nibs.obj_letter{iobj_type})));
            end
            disp(sprintf('total number of objects: %5.0f',length(unique_objs)));
            obj_counts=zeros(nPIDs_used,nobj_types);
            for PID_ptr=1:nPIDs_used
                iPID=PID_unique(PID_ptr);
                for iobj_type=1:nobj_types
                    patchpair_pid_type=intersect(find(obj_table.pid_num==iPID),find(obj_table.type==opts_nibs.obj_letter{iobj_type})); %points to patch pairs with this PID and object type
                    obj_counts(PID_ptr,iobj_type)=length(unique(obj_table.obj_num(patchpair_pid_type))); %how many different objects are there?.
                end
            end
            disp(sprintf(cat(2,'          PID       tot obj',sprintf(' %12s',opts_nibs.obj_desc{:}))));
            if (opts_show.if_objlist)
                for PID_ptr=1:nPIDs_used
                    disp(sprintf(' %12.0f',PID_unique(PID_ptr),sum(obj_counts(PID_ptr,:)),obj_counts(PID_ptr,:)));
                end
                disp(' ');
            end
            disp(sprintf(cat(2,'      minimum',sprintf(' %12.0f',min(sum(obj_counts,2)),min(obj_counts,[],1)))));
            disp(sprintf(cat(2,'      maximum',sprintf(' %12.0f',max(sum(obj_counts,2)),max(obj_counts,[],1)))));
            disp(sprintf(cat(2,'        total',sprintf(' %12.0f',sum(obj_counts(:)),sum(obj_counts,1)))));
            if (opts_show.if_objlist)
                disp(' ');
            end
        end
        image_file_unique=unique(meta_use.image_file);
        disp(sprintf('PIDs         range from %14.0f to %14.0f; %14.0f unique PIDs out of %6.0f available',...
            min(meta_use.PID_num),max(meta_use.PID_num),nPIDs_used,meta_use.nsubjs));
        disp(sprintf('      image files range from %14s to %14s; %14.0f unique image files analyzed',...
            meta_use.image_file{1},meta_use.image_file{end},length(image_file_unique)));
        disp(sprintf('patch pair counts range from %14.0f to %14.0f (mean %6.1f) for      all (%6.0f) available images',...
            min(meta_use.patch_count_table),max(meta_use.patch_count_table),...
            mean(meta_use.patch_count_table),meta_use.nsubjs));
        disp(sprintf('patch pair counts range from %14.0f to %14.0f (mean %6.1f) for analyzed (%6.0f) PIDs',...
            min(meta_use.patch_count_table(PID_unique)),max(meta_use.patch_count_table(PID_unique)),...
            mean(meta_use.patch_count_table(PID_unique)),length(PID_unique)));
        %view and gallery fields are for compatibility; all should be 1
        disp(sprintf('view numbers range from %14.0f to %14.0f; %14.0f unique view numbers',...
            min(meta_use.View_num),max(meta_use.View_num),length(unique(meta_use.View_num))));
        disp(sprintf('gal  numbers range from %14.0f to %14.0f; %14.0f unique gal  numbers',...
            min(meta_use.Gal_num),max(meta_use.Gal_num),length(unique(meta_use.Gal_num))));
        %
        disp('patch statistics:      Position            Length            Area');
        XYArea=prod(meta_use.XYLen,2);
        disp('                      X       Y           X       Y');
        disp(sprintf(' min:           %7.0f %7.0f     %7.0f %7.0f   %12.0f',min(meta_use.XYPos),min(meta_use.XYLen),min(XYArea)));
        for iqr=1:length(quantile_list)
            disp(sprintf(' quant %6.3f:  %7.0f %7.0f     %7.0f %7.0f   %12.0f',quantile_list(iqr),...
                quantile(meta_use.XYPos,quantile_list(iqr)),...
                quantile(meta_use.XYLen,quantile_list(iqr)),...
                quantile(XYArea,quantile_list(iqr))));
        end
        disp(sprintf(' max:           %7.0f %7.0f     %7.0f %7.0f   %12.0f',max(meta_use.XYPos),max(meta_use.XYLen),max(XYArea)));
        disp(sprintf('    mean:        %7.1f  %7.1f    %7.1f  %7.1f   %12.1f',mean(meta_use.XYPos),mean(meta_use.XYLen),mean(XYArea)));
        disp(sprintf(' geomean:        %7.1f  %7.1f    %7.1f  %7.1f   %12.1f',geomean(meta_use.XYPos),geomean(meta_use.XYLen),geomean(XYArea)));
        disp(sprintf(' std dev:        %7.1f  %7.1f    %7.1f  %7.1f   %12.1f',std(meta_use.XYPos),std(meta_use.XYLen),std(XYArea)));
        %survey for each value of patchProperties and algoParams
        patches_per_PID=zeros(1,nPIDs_used);
        for PID_ptr=1:nPIDs_used
            iPID=PID_unique(PID_ptr);
            seg=bpdata.segmentation{iPID}.patches{i_aspect,i_smooth};
            patches_per_PID(PID_ptr)=size(seg.patch,2);
        end
        patches_total=sum(patches_per_PID);
        disp(sprintf(' total number of patches: %6.0f',patches_total));
        %sanity checks on total patches
        if patches_total~=2*meta_use.npatches 
            warnings=strvcat(warnings,sprintf('%s: patches from PIDs is %5.0f but npatches is %5.0f pairs',...
                as_string,patches_total,meta_use.npatches));
        end
        if patches_total~=2*sum(meta_use.patch_count_table)
            warnings=strvcat(warnings,sprintf('%s: patches from PIDs is %5.0f but sum of patch count table is %5.0f pairs',...
                as_string,patches_total,sum(meta_use.patch_count_table)));
        end
        region_types=fieldnames(opts_nibs.mask_dict_composite);
        nregion_types=length(region_types);
        region_sizes=zeros(patches_total,nregion_types);
        for PID_ptr=1:nPIDs_used
            iPID=PID_unique(PID_ptr);
            seg=bpdata.segmentation{iPID}.patches{i_aspect,i_smooth};
            npairs=size(seg.pairedPatches,2);
            if patches_per_PID(PID_ptr)~=2*npairs %number of patches must be twice patch pairs
                warnings=strvcat(warnings,sprintf('%s, PID %4.0f: %4.0f patches but %4.0f patch pairs',...
                    as_string,iPID,patches_per_PID(PID_ptr),npairs));
            end
            patchMaps=cell(0);
            if iscell(seg.patchMap) %old versions: patchMap is an array, now is a cell array of maps
                patchMaps=seg.patchMap;
            else
                patchMaps{1}=seg.patchMap;
            end
            nobjs=length(patchMaps); %may include maps with no patches
            patch_pairs_each=zeros(1,nobjs); %number of patch pairs in each object
            pair_string=[];
            map_unique=0;
            for iobj=1:nobjs
                if ~isempty(patchMaps{iobj})
                    patch_pairs_each(iobj)=(length(unique(patchMaps{iobj}(:)))-1)/2;
                    map_unique=unique([map_unique,2*sum(patch_pairs_each([1:iobj-1]))+unique(patchMaps{iobj}(:))']);
                end
            end
            if (opts_show.if_objlist)
                disp(sprintf('PID %5.0f has %3.0f objects, with patch pair counts: %s',iPID,nobjs,sprintf(' %5.0f ',patch_pairs_each)));
            end
            if meta_use.patch_count_table(iPID)~=sum(patch_pairs_each)
                warn_string=sprintf('%s, PID %4.0f: patch pairs in patch count table (%4.0f) and all patch maps (%4.0f) disagree',...
                    as_string,iPID,meta_use.patch_count_table(iPID),sum(patch_pairs_each));
                warnings=strvcat(warnings,warn_string);
            end
            %if object type fields are available, then check that objects are consistent
            if isfield(meta_use,'obj_ids') & isfield(meta_use,'obj_type')
                PID_patchpairs=find(meta_use.PID_num==iPID);
                PID_objs_unique=unique(meta_use.obj_ids(PID_patchpairs));
                if length(PID_objs_unique)~=sum(patch_pairs_each>=1)
                   warn_string=sprintf('%s, PID %4.0f: unique objects from obj_ids (%4.0f) and number of non-empty patch maps (%4.0f) disagree',...
                       as_string,iPID,length(PID_objs_unique),sum(patch_pairs_each>=1));
                   warnings=strvcat(warnings,warn_string);
                end
                if (opts_show.if_objlist)
                    disp(sprintf(' objects with >0 patch pairs: %s',sprintf('%8s',PID_objs_unique{:})));
                end
                %check that each map has patches of only one type
                for iobj=1:nobjs
                    PID_patchpairs_obj=min(PID_patchpairs)+sum(patch_pairs_each(1:iobj-1))+[0:patch_pairs_each(iobj)-1];
                    if ~isempty(PID_patchpairs_obj)
                        PID_obj_unique=unique(meta_use.obj_ids(PID_patchpairs_obj));
                        if length(PID_obj_unique)~=1
                           warn_string=sprintf('%s, PID %4.0f: map %2.0f has %4.0f patch pairs, which come from %2.0f different objects: %s',...
                               as_string,iPID,iobj,length(PID_patchpairs_obj),length(PID_obj_unique),sprintf('%8s',PID_obj_unique{:}));
                           warnings=strvcat(warnings,warn_string);
                        end
                    end
                end
            end
            %map_unique=unique(patchMaps{1}); %modified below added 01Jun22
            if length(map_unique)~=1+2*sum(patch_pairs_each)
                warn_string=sprintf('%s, PID %4.0f: total patch pairs across patch maps is %5.0f but with offsets there are %5.0f unique values in indiv maps',...
                    as_string,iPID,1+2*sum(patch_pairs_each),length(map_unique));
                warnings=strvcat(warnings,warn_string);
            end
            %number of patches must correspond to unique map entries -- 
            if (min(map_unique)~=0) | (max(map_unique)~=patches_per_PID(PID_ptr)) | (length(map_unique)~=(1+patches_per_PID(PID_ptr)))
                warn_string=sprintf('%s, PID %4.0f: patch map should have %5.0f unique values but has %5.0f',...
                    as_string,iPID,1+patches_per_PID(PID_ptr),length(map_unique));
                warn_string=cat(2,warn_string,';  missing from map:',sprintf(' %4.0f',setdiff([0:patches_per_PID(PID_ptr)],map_unique)));
                warn_string=cat(2,warn_string,';  extra   in   map:',sprintf(' %4.0f',setdiff(map_unique,[0:patches_per_PID(PID_ptr)])));
                warnings=strvcat(warnings,warn_string);
            end
            %now check each patch
            for ipair=1:npairs
                %look it up to find the overall patch number in the overall metadata
                meta_nums=intersect(find(iPID==meta_use.PID_num),find(ipair==meta_use.Examp_num));
                if length(meta_nums)~=1
                    warnings=strvcat(warnings,sprintf('%s, PID %4.0f pair %4.0f: found %2.0f matches',as_string,iPID,ipair,length(meta_nums)));
                else
                    metaPosLen=[meta_use.XYPos(meta_nums,:)   ,meta_use.XYLen(meta_nums,:)];
                    segPosLen =[seg.pairedPatches{ipair}.XYPos,seg.pairedPatches{ipair}.XYLen];
                    if ~all(metaPosLen==segPosLen)
                        warnings=strvcat(warnings,sprintf('%s, PID %4.0f pair %4.0f: XYPos, XYLen from meta=%4.0f %4.0f %4.0f %4.0f from pairedPatches=%4.0f %4.0f %4.0f %4.0f',...
                           as_string,iPID,ipair,metaPosLen,segPosLen));
                    end
                    region_sizes_check=zeros(1,nregion_types);
                    %compare region sizes in individual maps with region sizes in composite mask
                    for iregion=1:nregion_types
                        region_sizes(meta_nums,iregion)=sum(double(seg.pairedPatches{ipair}.mask.composite(:)==iregion-1));
                        if isfield(seg.pairedPatches{ipair}.mask,region_types{iregion})
                            region_sizes_check(iregion)=sum(double(seg.pairedPatches{ipair}.mask.(region_types{iregion})(:)));
                        end
                        if region_sizes(meta_nums,iregion)~=region_sizes_check(1,iregion)
                            warnings=strvcat(warnings,sprintf('%s, PID %4.0f pair %4.0f region sizes disagree for %20s: composite: %5.0f indiv %5.0f',...
                                as_string,iPID,ipair,region_types{iregion},region_sizes(meta_nums,iregion),region_sizes_check(1,iregion)));
                        end
                    end
                    %verify that total of region sizes is size of mask
                    if sum(region_sizes(meta_nums,:))~=prod(seg.pairedPatches{ipair}.XYLen+1)
                        warnings=strvcat(warnings,sprintf('%s, PID %4.0f pair %4.0f: sum of region sizes (%6.0f) and prod(XYLen+1) (%4.0f +1)*(%4.0f+1)  disagree',...
                            as_string,iPID,ipair,sum(region_sizes(meta_nums,:)),seg.pairedPatches{ipair}.XYLen));
                    end
                    %verify that all masks have the correct size
                    mask_fields=fieldnames(seg.pairedPatches{ipair}.mask);
                    for ifield=1:length(mask_fields)
                        mask_size=size(seg.pairedPatches{ipair}.mask.(mask_fields{ifield}));
                        if all(fliplr(mask_size)~=seg.pairedPatches{ipair}.XYLen+1) %note that masks have X and Y reversed w.r.t. XYLen and metadata
                            warnings=strvcat(warnings,sprintf('%s, PID %4.0f pair %4.0f: mask size for %20s ([%4.0f %4.0f]) and XYLen+1 ([%4.0f %4.0f]+1) disagree',...
                            as_string,iPID,ipair,mask_fields{ifield},mask_size,seg.pairedPatches{ipair}.XYLen));
                        end
                    end %ifield
                end %unique match found
            end %ipair
        end %PID_ptr
    end %smoothness_ptr
end %aspect_ptr
return

