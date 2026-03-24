%ffdm_btc_calc_gen: basic binary texture coordinate analysis of ffdm dataset
%
% Does analysis for a fixed value of the product N (downsampling) x R (downsampled blocks per patch)
% differs from ffdm_btc_calc in that it allows for non-constant values of NR.
%
% A list of log2(N) and log2(NR) is obtained, and the largest value of NR
% is used to subdivide into patches.  So analyses with smaller values of NR
% use exactly the same regions of the image.
%
% After this is run, can run in same (or saved) workspace:
%     ffdm_btc_scat_gen, to plot scattergrams of btc statistics and correlations with clinical metadata
%     ffdm_btc_plot_gen, to plot btc statistics means and standard devs
%     ffdm_btc_spec_gen: to plot and analyze statistics of spectra
%
% options:
%   padding 
%   to average according to number of samples per image, or per per subject
%   to use spectrum across entire dataset, across this view, or private
%   margin within each patch for calculating binary texture stats
%
%   02Apr18: option to subtract pixel-wise mean, and option to subtract within-ROI mean (latter had been obligate)
%   03Apr18: option for windowing with cosine bell
%   11Apr18: option for frozen randomization
%   11Apr18: added stdv (i.e., s.e.m) and confidence limits by jackknifing per subject
%   11Apr18: changed defaults to margin=1 (was 0) and cosine bell=1 (was no windowing)
%   18Apr18: allow for variable N*R
%   30Apr18: allow for limit on number of ROIs to analyze
%   30Apr18: allow for limit on number of individual patch spectra to plot
%   01May18: add expected contribution to standard deviation based on counting errors
%   19Nov19: begin option to get input from acmd database
%      data_source=1->original "pilot" ffdm database; they are cut into square ROI's here
%      data_source=2->image patches cut into square ROI's by CKA
%      data_source=3->image patches stats manipulated by mlis_alglib_pilot or similar, in a "sdb" database
%   22Nov19: fixed a bug in calculation of jackniife confidence limits on patch-weighted stdv across all views (btc_std_allv_drop)
%   24Dec19: begin option to get input from ultrasound database (data_source=4)
%   26Dec19: make plotting optional, fixed bug that skipped detailed plots of per-patch spectra 
%   26Dec19: retain ps_avg_v for ffdm_spec_stats, and fix bug in loglog plots
%   30Dec19: begin option to get input from MRI databases (data_source=5)
%   15Jan20: allow for an updated patch metadata file to be read for data_source=3
%   21Jan20: change default value of r_min_log2 to 4 from 5, to allow for regions as coarse as 16x16; add per-patch whitening option
%   07Feb20: call to getcorrs_p2x2 changed to compute only btc stats, to save time
%   13Apr20: begin option to get input from Penn natural image database (data_source=6)
%   08Aug20: option to not whiten before calculating btc stats (spectra still calculated), whiten_avg_type_use=7
%   10Aug20: begin option for image statistics to be calculated at power-of-2 scale after downsampling (S), and to step the start phase
%      * implemented via S_list{iNR} and Sstep_list{iNR}
%      * patches_binarized now a cell array of cell arrays -- outer dim for iNR, inner for start phase, and uint8 to save space
%      * binarization based on > median, rather than >= median, for consistency with mlis_alglib_pilot
%   16Sep20: use mlis_window_setup to create cosine bell.
%   23Sep20: add options for whitening via a spectrum of defined power law, whiten_avg_type_use=8
%   27Oct21: add documentation: re ffdm_upmc_nvms_calc, which uses data_source=7,8; and re patches*
%   20Jun22: begin addition of natural image database (nibs) with boundary patching, data_source=9
%      * modularize linear and log transform of nimg, nibs
%      * regtype_count: number of region types (e.g., figure, ground, all) -- for data_source=9, otherwise 1
%      * mask for region type (figure, ground), etc is retrieved three times:  for testing, at time of
%        binarization (to determine cutpoint) and at time of calculation of local image statistics
%      * can customize nibs_meta_path_def, nibs_meta_file_def, nibs_image_path_def
%      NOTE analysis by region type not yet done -- need to retrieve mask
%      and apply it, and use only masked region for block counts, and
%      evetnaully calculate covariances between figure and neighboring
%      ground
%   15Aug22: begin incorporation of bone radiographs and scintigraphs data -- data_source=10
%      * sources from xls sheets from Elena Ajayi: one subject may have several ROI files, each differnt sizes
%      *    roi_sizes_all (for each subject may have several files, with indiv ROIs)
%      *    patch_whichroi: list of rois for each patch
%      *    can define bxrs_path and bxrs_files to customize file location and xlsx databases
%      Also:
%      * fix a bug if NR_exclude is not defined
%      * fix a bug in the prompt for patches_usefrac
%      * fix a bug in patch_starts
%      * fix a bug in display of number of examples in statistical tables
%   24Sep22: use bone_scint_select
%   27Oct22: fix bug related to n_patches_eachdim when no patches are available
%   29Nov22: fix bug related to n_patches when no patches are available.
%   
%  See also:  FFDM_EXTRACT_DEMO, PREWHITEN_MRI_DEMO, BTC_DEFINE, GLIDER_MAPUPB, GETCORRS_P2X2, BTC_CORRS2VEC, BTC_VEC2LETCODE,
%  FFDM_BTC_SCAT,FFDM_BTC_CALC, FFDM_BTC_SCAT_GEN, JACK, TINV, GEOMOD_BTC_DEMO, FFDM_SELECT_PATCHES, FFDM_READMETADATA,
%  FFDM_READ_SUBPATCH, FFDM_BTC_SPEC_GEN, USIC_READ_ALL, MRIX_READ_ALL, MLIS_ALGLIB_FILENAME_AUG,
%  NIMG_READ_PATCHES, NIMG_READ_METADATA, NIMG_SELECT_PATCHES, MLIS_WINDOW_SETUP, FFDM_UPMC_NVMS_CALC.
%  NIBS_BPDATA_DEMO, NIBS_BPDATA_SHOW, NIBS_GET_PATCHPAIR, NIBS_GET_META, FFDM_BTC_CALC_LINLOG, NIBS_BPDATA_PRUNE,
%  BONE_READ_XLS, BONE_BTC_DEMO, BONE_SCINT_SELECT.
%
data_sources=cell(0);
%
data_sources{1}.label='pilot ffdm dataset';
data_sources{1}.roi_range=[0 4095]; %raw tiff
data_sources{1}.pixel_microns=75;
data_sources{1}.type='ffdm';
data_sources{1}.if_makepatches=1; %rois need to be cut into patches
%
data_sources{2}.label='full ffdm dataset for acmd, unprocessed';
data_sources{2}.roi_range=[0 1]; %rescaled tiff
data_sources{2}.pixel_microns=75;
data_sources{2}.type='ffdm';
data_sources{2}.if_makepatches=0; %rois DO NOT need to be cut into patches
%
data_sources{3}.label='subset of ffdm dataset for acmd, after mlis processing';
data_sources{3}.roi_range=[0 255]; %raw bmp
data_sources{3}.pixel_microns=75; 
data_sources{3}.type='ffdm';
data_sources{3}.if_makepatches=0; %rois DO NOT need to be cut into patches
%
data_sources{4}.label='carotid ultrasound';
data_sources{4}.roi_range=[0 1]; %rescaled bmp
data_sources{4}.pixel_microns=NaN; %unknown
data_sources{4}.type='usic';
data_sources{4}.if_makepatches=1; %rois need to be cut into patches
%
data_sources{5}.label='brain MRI (T1: ADNI, T1: OASIS, T2: TNS HV and Patients)';
data_sources{5}.roi_range=[NaN NaN]; %range depends on MRI database
data_sources{5}.pixel_microns=1000; %MRI slices are 1 x 1 mm
data_sources{5}.type='mrix';
data_sources{5}.if_makepatches=0; %rois DO NOT need to be cut into patches
%
data_sources{6}.label='Penn Natural Image Database';
data_sources{6}.roi_range_init=[0 2^16-1]; %16 bits of luminance resolution, will be set to [0 1]
data_sources{6}.roi_range=[0 1]; %rescaled after linear or log transform
data_sources{6}.pixel_microns=NaN; %unknown
data_sources{6}.type='nimg';
data_sources{6}.if_makepatches=0; %rois do not need to be cut into patches
%
%data_sources{7} and data_sources{8} used by ffdm_upmc_nvms_calc for upmc
%(mammograms) or nvms (Gaussian surrogates) -- btc statistics of patches calculated externally 
% and put into form for ffdm_ btc_plot_gen by ffdm_upmc_nvms_calc
%
data_sources{9}.label='Berkeley Img Dbase, boundary parsing';
data_sources{9}.roi_range_init=[0 2^8-1]; %8 bits of luminance resolution after extraction by nistats-process, will be set to [0 1]
data_sources{9}.roi_range=[0 1]; %rescaled after linear or log transform
data_sources{9}.pixel_microns=NaN; %unknown
data_sources{9}.type='nibs';
data_sources{9}.if_makepatches=0; %rois do not need to be cut into patches
data_sources{9}.imsize=[321 481]; %image size in pixels, can be landscape or portrait, needed to determine maximum patch size
data_sources{9}.if_regsel=1; %regions can be selected
%
%
data_sources{10}.label='bone xray or scintigram';
data_sources{10}.roi_range=[0 2^8-1]; %jpeg
data_sources{10}.pixel_microns=NaN; %unknown
data_sources{10}.type='bxrs';
data_sources{10}.if_makepatches=1; %rois need to be cut into patches
%

%defaults for nibs
if ~exist('opts_nibs') opts_nibs=struct; end
opts_nibs=nibs_defopts(opts_nibs);
if ~exist('opts_nibs_show') opts_nibs_show=struct; end
if ~exist('opts_nibs_patchpair') opts_nibs_patchpair=struct; end
if ~exist('nibs_meta_path_def') nibs_meta_path_def='C:\Users\jdvicto\Documents\jv\EY7977\nistats\nibs-output\';end
if ~exist('nibs_meta_file_def') nibs_meta_file_def='bpdata_BSDS_21Jun22.mat';end
if ~exist('nibs_image_path_def') nibs_image_path_def='C:\Users\jdvicto\Dropbox\ObjectBorderStats\Images'; end
if ~exist('nibs_NRpatch_log2def') nibs_NRpatch_log2def=7; end %default patch size for whitening
%
for ids=1:length(data_sources)
    if ~isempty(data_sources{ids})
        data_sources{ids}=filldefault(data_sources{ids},'if_regsel',0);
        disp(sprintf('%2.0f ->%6s data from %s, region selection flag %1.0f',...
            ids,data_sources{ids}.type,data_sources{ids}.label,data_sources{ids}.if_regsel));
    end
end
ifok=0;
while ifok==0
    data_source=getinp('choice','d',[1 length(data_sources)]);
    ifok=~isempty(data_sources{data_source});
end
if_regsel=data_sources{data_source}.if_regsel; %region selection?
if (if_regsel==0) %default values if region selection is inactive
    regtype_count=1;
    regtype_masks{1}=0;
    regtype_descs{1}='n/a';
end  
%
if_clear=getinp('1 to clear large arrays after use','d',[0 1],1);
if_norand=getinp('-1 for fixed randomization, 0 for full randomization, 1 for none','d',[-1 1],0);
if (if_norand==-1)
    rand('state',0);
end
%
if ~exist('if_plot_pipeline') if_plot_pipeline=0; end
if_plot_pipeline=getinp('1 to plot pipeline','d',[0 1],if_plot_pipeline);
if ~exist('if_plot_spectra') if_plot_spectra=0; end
if_plot_spectra=getinp('1 to plot spectra','d',[0 1],if_plot_spectra);
if ~exist('if_plot_details') if_plot_details=0; end
if_plot_details=getinp('1 to plot spectra anisotropy details','d',[0 1],if_plot_details);
%
%image data 
%
roi_range=data_sources{data_source}.roi_range;
if ~exist('pixel_microns') pixel_microns=data_sources{data_source}.pixel_microns;end
if ~exist('image_matfile') image_matfile='ffdm_extract_demo.mat'; end %alternatively: image_matfile='L:\AbbeyResources\ffdm_extract\demo.mat'
%
%analysis parameters -- fft
%
if ~exist('fft2_padfactor') fft2_padfactor=1; end %no padding
if ~exist('if_spec_subroimean') if_spec_subroimean=1; end %subtract the ROI mean before computiong spectra
if ~exist('if_spec_subpxlmean') if_spec_subpxlmean=0; end %subtract pixel mean before computiong spectra
if ~exist('if_spec_cosbell') if_spec_cosbell=1; end; %apply cosine bell window
%
%plotting parameters
%
if ~exist('view_colors') view_colors={'r','b','m','c'}; end %for LCC, LMLO, RCC, RMLO
if ~exist('ps_minfac') ps_minfac=10^-5; end %relative size of smallest power spectrum shown
ps_log10range=-log10(ps_minfac);
hv_labels={'horiz','vert'};
if ~exist('plot_spec_thinfac') plot_spec_thinfac=1; end %set to > 1 to plot only a fraction of spectra of individual patches
if ~exist('plot_spec_maxnum')  plot_spec_maxnum=2000; end %maximum numbe of spectra of individual patches to plot
%
%analysis parameters -- averaging methods
%
ps_avg_types_labels={'per patch','per subj'}; %ways to average the power spectra across patches
ps_avg_types=length(ps_avg_types_labels);
ps_avg_types_symbs={'-',':'};
%
btc_avg_types_labels={'per patch','per subj'}; %ways to average binary texture coords across patches
btc_avg_types=length(btc_avg_types_labels);
btc_avg_types_symbs={'-',':'};
%
whiten_avg_types_labels={'per subj and view','per view (patch wtd)','per view (subj wtd)','global (patch wtd)','global (subj wtd)','per patch','none','power law'};
whiten_avg_types_labels_orig=whiten_avg_types_labels; %so that 'power law' can be modified
whiten_avg_types=length(whiten_avg_types_labels);
%
% for binary texture coordinates
%
btc_checkdef=[0 0;0 1;1 0;1 1];
btc_ng=2;
btc_dict=btc_define;
btc_reshape=repmat(btc_ng,1,size(btc_checkdef,1)); %typically, [2 2 2 2];
btc_nconfigs=btc_ng.^size(btc_checkdef,1); %typically, 16=2^4 
btc_n=length(btc_dict.codel); %number of coords, typically 10
btc_order_JV=btc_dict.codel; %gbcdetuvwa
btc_order_HX='cbdevwtua'; %parameter order used in Hermundstad et al. and Xu et al.
%
if ~exist('whiten_avg_type_use') whiten_avg_type_use=3; end
if ~exist('R_min_log2') R_min_log2=4; end
if ~exist('margin_width') margin_width=1; end
R_min=2.^R_min_log2;
preprocess_downsample=1; %default; may be reset if images are preprocessed prior to the files that are reead, e.g., with data_source=3
%
% regtion types -- will be overridden for data_source=9
%
regtype_count=1; %number of region types
regtype_masks=cell(regtype_count,1);
regtype_descs=cell(regtype_count,1);
regtype_masks{1}=-1; %no masking
regtype_descs{1}='n/a';
regtype_binarize_useall=1; %1 to use the entire patch for binarizing
%
switch data_source
    case 1 %ffdm
        disp('image_matfile -- if not present, look on external drives')
        image_matfile
        load(image_matfile,'roi*','file_nos','view*');
        image_matfile_start=max([strfind(image_matfile,'\'),strfind(image_matfile,'/'),0]);
        image_matfile_short=image_matfile(image_matfile_start+1:end);
        image_matfile_short=strrep(image_matfile_short,'.mat','');
        count_subjs=size(roi_imgs,1);
        count_views=size(roi_imgs,2);
        datasource_string=image_matfile_short;
        %
    case {2,3} %ffdm, but rois already cut into squares and in mat-files
        %file locations for stimulus database
        [dos_err,user_profile]=dos('echo %USERPROFILE%');
        user_profile=deblank(user_profile);
        if ~exist('patch_metadata_filename') patch_metadata_filename='Patch_Metadata_Annotated_10Jan20'; end
        if ~exist('patch_metadata_path') patch_metadata_path='\Documents\JV\encl\Abbey\'; end
        %
        if (data_source==2) %acmd: 2 -> direct from patch database (read from patch_metadata), 
            %
            if ~exist('patch_path') patch_path='L:\AbbeyResources\Patches\'; end
            patch_path=strrep(getinp('patch path (for reading)','s',[],strrep(patch_path,'\','/')),'/',filesep);
            datasource_string=patch_metadata_filename;
            patch_metadata_filename=strrep(getinp('patch metadata file (for reading)','s',[],strrep(patch_metadata_filename,'\','/')),'/',filesep);
            patch_metadata_fullname=cat(2,user_profile,patch_metadata_path,filesep,patch_metadata_filename);
            patch_metadata=ffdm_readmetadata(patch_metadata_fullname);
        else  % 3-> from subpatches processed by mlis_alglib_pilot (read from sdb_subpatch files)
            if ~exist('sdb_path') sdb_path='\Documents\JV\encl\Abbey\'; end
            sdb_path=strrep(getinp('sdb path (for reading','s',[],strrep(sdb_path,'\','/')),'/',filesep);
            if ~exist('sdb_filename') sdb_filename='subpatch_sdb_07Nov19'; end
            sdb_filename=strrep(getinp('sdb file name(for reading)','s',[],strrep(sdb_filename,'\','/')),'/',filesep);
            sdb_fullname=cat(2,sdb_path,filesep,sdb_filename);
            sdb=getfield(load('subpatch_sdb_07Nov19'),'sdb');
            %
            disp('patch metadata extracted from stimulus database:');
            disp(sdb.patch_metadata)
            if getinp('1 to update patch metadata to include quality','d',[0 1],~isfield(sdb.patch_metadata,'Quality'))
                patch_metadata_filename=strrep(getinp('patch metadata file (for reading)','s',[],strrep(patch_metadata_filename,'\','/')),'/',filesep);
                patch_metadata_fullname=cat(2,user_profile,patch_metadata_path,filesep,patch_metadata_filename);
                patch_metadata=ffdm_readmetadata(patch_metadata_fullname);
                disp('patch metadata updated:')
            else
                patch_metadata=sdb.patch_metadata;
            end
            %
            if ~exist('prefix') prefix='Aux_s'; end
            subinfo=acmd_mlis_subinfo(sdb.stimulus_path,prefix);
            check_if_verbose=0;
            check_filenums='ends';
            check_if_ask=1; %if we check files, want option to change path
            [sdb,ifok,sdb_warnings,sdb_warnings_file,subinfo]=mlis_alglib_sdbcheck(sdb,check_filenums,check_if_ask,check_if_verbose,subinfo,prefix);
            %
            sdba=[];
            sdba.nstimuli=length(sdb.subpatch.alg_ptr);
            disp('algorithm list')
            for ialg=1:length(sdb.algs)
                disp(sprintf('alg %2.0f->%s',ialg,sdb.algs{ialg}.desc));
            end
            sdba.alg_ptr_avail=unique(sdb.subpatch.alg_ptr(:));
            sdba.alg_ptr_sel=getinp('algorithm pointer(s)','d',[min(sdba.alg_ptr_avail) max(sdba.alg_ptr_avail)]);
            sdba.alg_ptr_string=sprintf(' %1.0f',sdba.alg_ptr_sel);
            sdba.alg_ptr_sel=intersect(sdba.alg_ptr_sel,sdba.alg_ptr_avail);
            stimuli_alg_sel=find(ismember(sdb.subpatch.alg_ptr(:),sdba.alg_ptr_sel));
            disp(sprintf('Alg      choice selects %5.0f out of %5.0f stimuli',length(stimuli_alg_sel),sdba.nstimuli));
            %
            disp('target list');
            disp(sdb.btc_target_list(sdb.btc_target_ptrs,:));
            sdba.tvp_ptr_avail=unique(sdb.subpatch.target_vec_avail_ptr(:));
            sdba.tvp_ptr_sel=getinp('target vec pointer(s)','d',[min(sdba.tvp_ptr_avail) max(sdba.tvp_ptr_avail)]);
            sdba.tvp_ptr_string=sprintf(' %1.0f',sdba.tvp_ptr_sel);
            sdba.tvp_ptr_sel=intersect(sdba.tvp_ptr_sel,sdba.tvp_ptr_avail);
            stimuli_tvp_sel=find(ismember(sdb.subpatch.target_vec_avail_ptr(:),sdba.tvp_ptr_sel));
            disp(sprintf('Target  choice selects %5.0f out of %5.0f stimuli',length(stimuli_tvp_sel),sdba.nstimuli));
            stimuli_alg_tvp_sel=intersect(stimuli_alg_sel,stimuli_tvp_sel); %stimuli with the chosen algorithm and target vector
            disp(sprintf('Alg and target selects %5.0f out of %5.0f stimuli',length(stimuli_alg_tvp_sel),sdba.nstimuli));
            %
            datasource_string=cat(2,sdb_filename,' alg ptrs ',sdba.alg_ptr_string,' tgt vec ptrs',sdba.tvp_ptr_string);
            %
            %keep these variables so we can clear the large sdb
            sdba.preprocess=sdb.downsample;
            sdba.stimulus_path=sdb.stimulus_path;
            sdba.subpatch_patch_num=sdb.subpatch.patch_num;
            sdba.subpatch_filename=sdb.subpatch.filename;
            sdba.algs=sdb.algs;
            sdba.btc_target_ptrs=sdb.btc_target_ptrs;
            sdba.btc_target_list=sdb.btc_target_list;
            if if_clear
                clear sdb;
            end
            preprocess_downsample=sdba.preprocess;
        end
        %code common to data_source=2 or 3 to tally rois (patches or stimuli_ for each 
        %
        [patch_list,sel_used,avail,patch_list_sel]=ffdm_select_patches(patch_metadata);
        patch_PID_nums=patch_metadata.PID_num(patch_list);
        patch_View_nums=patch_metadata.View_num(patch_list);
        view_names=patch_metadata.view_names;
        file_nos=unique(patch_PID_nums(:)');
        view_nos=unique(patch_View_nums(:)');
        count_subjs=length(file_nos);
        count_views=length(view_nos);
        %rois are from patches (data_source=2) or stimuli (data_source=3_
        roi_count_sdb=zeros(count_subjs,count_views); %how many patches for each unique subject and view
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch sizes, should all be identical but NaN if no view for that subject
        for isubj=1:count_subjs
            for iview=1:count_views
                patch_sel=intersect(patch_list(patch_PID_nums==file_nos(isubj)),patch_list(patch_View_nums==view_nos(iview)));
                if (data_source==3)
                    stimuli_patch_sel=find(ismember(sdba.subpatch_patch_num,patch_sel)); %which stimuli come from a patch with the chosen subject and view
                    stimuli_sel=intersect(stimuli_patch_sel,stimuli_alg_tvp_sel);
                    roi_count_sdb(isubj,iview)=length(stimuli_sel);
                    roi_list_sdb{isubj,iview}=stimuli_sel;
                else
                    roi_count_sdb(isubj,iview)=length(patch_sel);
                    roi_list_sdb{isubj,iview}=patch_sel;
                end
                if length(patch_sel)>0
                    roi_sizes(isubj,iview,:)=reshape(min(patch_metadata.XYLen(patch_sel,:),[],1),[1 1 2])/preprocess_downsample; %roi size is minimum of roi sizes of all patches from this subject and view
                end
            end
        end
    case 4 %carotid ultrasound
        if ~exist('usic_path') usic_path=[]; end
        if ~exist('usic_file') usic_file=[]; end
        if ~exist('usic_preproc') usic_preproc=struct(); end
        if ~exist('usic_examp1_only') usic_examp1_only=1; end
        [usic_patches,usic_subpatches_all,usic_meta_all,usic_preproc_used]=usic_read_all(usic_path,usic_file,usic_preproc);
        disp('Ultrasound data read.');
        datasource_string=strrep(usic_meta_all.usic_file,'.mat.','');
        preprocess_downsample=usic_preproc_used.downsample;
        %
        usic_examp1_only=getinp('1 to keep example 1 only','d',[0 1],usic_examp1_only);
        usic_meta=struct();
        if (usic_examp1_only)
            usic_rois_sel=find(usic_meta_all.Examp_num==1);
            usic_meta.patch_count_table=min(usic_meta_all.patch_count_table,1);
        else
            usic_rois_sel=[1:length(usic_subpatches_all)];
            usic_meta.patch_count_table=usic_meta_all.patch_count_table;
        end
        %
        usic_meta.usic_path=usic_meta_all.usic_path;
        usic_meta.usic_file=usic_meta_all.usic_file;
        usic_meta.npatches=length(usic_rois_sel);
        usic_meta.views=usic_meta_all.views;
        usic_meta.PID_num=usic_meta_all.PID_num(usic_rois_sel);
        usic_meta.XYLen=usic_meta_all.XYLen(usic_rois_sel,:);
        usic_meta.XYPos=usic_meta_all.XYPos(usic_rois_sel,:);
        usic_meta.View_num=usic_meta_all.View_num(usic_rois_sel);
        usic_meta.Examp_num=usic_meta_all.Examp_num(usic_rois_sel);
        %
        usic_rois=cell(length(usic_rois_sel),1);
        for imptr=1:length(usic_rois_sel)
            isp=usic_rois_sel(imptr);
            usic_rois{imptr,1}=usic_subpatches_all{isp};
            usic_meta.PID_orig{imptr,1}=usic_meta_all.PID_orig{isp};
            usic_meta.PID{imptr,1}=usic_meta_all.PID{isp};
            usic_meta.View{imptr,1}=usic_meta_all.View{isp};
        end
        view_names=usic_meta.views;
        file_nos=usic_rois_sel;
        view_nos=unique(usic_meta.View_num);
        count_subjs=size(usic_meta.patch_count_table,1);
        count_views=length(view_nos);
        %
        roi_count_sdb=usic_meta.patch_count_table;
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch size, minimum for all patches for this view and this subject, should all be identical but NaN if no view for that subject
        for isubj=1:count_subjs
            for iview=1:count_views
                roi_sel=intersect(find(usic_meta.PID_num==isubj),find(usic_meta.View_num==iview));
                roi_count_sdb(isubj,iview)=length(roi_sel);
                roi_list_sdb{isubj,iview}=roi_sel;
                if length(roi_sel)>0
                    roi_sizes(isubj,iview,:)=usic_preproc_used.XYLen/preprocess_downsample; %all rois of same size, take downsampling into account
                end
            end
        end
        if (if_clear)
            clear usic_patches usic_subpatches_all
        end
    case 5 %MRI
        if ~exist('mrix_path') mrix_path=[]; end
        if ~exist('mrix_file') mrix_file=[]; end
        %call in special mode to get subdirectory options
        mrix_subdirs=mrix_read_all(mrix_path);
        for mrix_ds=1:length(mrix_subdirs)
            disp(sprintf(' %2.0f->%s',mrix_ds,mrix_subdirs{mrix_ds}));
        end
        mrix_ds=getinp('choice','d',[1 length(mrix_subdirs)]);
        [mrix_subpatches,mrix_meta,mrix_lims]=mrix_read_all(mrix_path,mrix_subdirs{mrix_ds},mrix_file);
        disp(sprintf('MRI data read.  %5.0f subjects, %5.0f slices, limits: [%8.3f %8.3f]',size(mrix_meta.patch_count_table,1),mrix_meta.npatches,mrix_lims));
        data_sources{5}.roi_range=mrix_lims; %range depends on MRI database
        roi_range=mrix_lims;
        datasource_string=strrep(mrix_meta.mrix_subdir,'.mat.','');
%        preprocess_downsample=1; set by default, so removed 27Jun22
        view_names=mrix_meta.views;
        file_nos=[1:length(mrix_meta.image_file)];
        view_nos=unique(mrix_meta.View_num);
        count_subjs=size(mrix_meta.patch_count_table,1);
        count_views=length(view_nos);
        %
        roi_count_sdb=mrix_meta.patch_count_table;
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch size, minimum for all patches for this view and this subject, should all be identical but NaN if no view for that subject
        for isubj=1:count_subjs
            for iview=1:count_views
                roi_sel=intersect(find(mrix_meta.PID_num==isubj),find(mrix_meta.View_num==iview));
                roi_count_sdb(isubj,iview)=length(roi_sel);
                roi_list_sdb{isubj,iview}=roi_sel;
                if length(roi_sel)>0
                    roi_sizes(isubj,iview,:)=min(mrix_meta.XYLen(roi_sel,:),[],1);
                end
            end
        end
    case 6 %natural images from Penn database
        if ~exist('nimg_meta_path') nimg_meta_path=[]; end % will default to 'C:\Users\jdvicto\Dropbox\Textures\RawNIstats\ImagePatchStats\data\'
        if ~exist('nimg_meta_file') nimg_meta_file=[]; end % will default to 'PPandNIstats_allAnalyses.mat'
        if ~exist('nimg_dir') nimg_dir=[]; end % will default to 'C:\Users\jdvicto\Desktop\PennNaturalImageDatabase\';
        if ~exist('NR_orig') NR_orig=[]; end %values of N and R used to determine focus; could choose NR_orig=[2 128]
        ifok=0;
        while (ifok==0) %N and R must be powers of 2
            nimg_meta_all=nimg_readmetadata(nimg_meta_path,nimg_meta_file,NR_orig);
            if all(2.^round(log(nimg_meta_all.NR_orig)/log(2))==nimg_meta_all.NR_orig)
                [nimg_patch_list,nimg_sel_used,nimg_avail,nimg_patch_list_sel]=nimg_select_patches(nimg_meta_all);
                if length(nimg_patch_list)>0
                    disp(sprintf(' %6.0f patches meet these criteria.',length(nimg_patch_list)));
                    ifok=getinp('1 if ok','d',[0 1]);
                else
                    disp(' No image patches meet these criteria.');
                end
            else
                disp(' N and R must be powers of 2.');
            end
        end
        [nimg_patches,nimg_dir_used,nimg_read_errs,nimg_imrange]=nimg_read_patches(nimg_meta_all,nimg_patch_list,nimg_dir);
        nimg_meta_all.nimg_dir=nimg_dir_used;
        disp(sprintf('Natural image data read (%3.0f of %3.0f galleries, %3.0f of %3.0f files): %8.0f patches of size %5.0f x %5.0f.',...
            length(unique(nimg_meta_all.Gal_num(nimg_meta_all.PID_num(nimg_patch_list)))),length(nimg_meta_all.gal_names_sorted),...
            length(unique(nimg_meta_all.PID_num(nimg_patch_list))),nimg_meta_all.nsubjs,length(nimg_patch_list),prod(nimg_meta_all.NR_orig),prod(nimg_meta_all.NR_orig)));
        %create a metadata structure just from the patches actually used, but retain full gallery information
        nimg_meta=struct();
        %fields to copy
        nimg_meta_fields={'nsubjs','nimg_meta_path','nimg_meta_file','nimg_dir','NR_orig','views','gal_names_sorted','gal_img_list','Gal_num'};
        for ifield=1:length(nimg_meta_fields)
            nimg_meta.(nimg_meta_fields{ifield})=nimg_meta_all.(nimg_meta_fields{ifield});
        end
        nimg_meta.patch_count_table=zeros(nimg_meta.nsubjs,1); %how many patches are used from each image (can be zero)
        for isubj=1:nimg_meta_all.nsubjs
            nimg_meta.patch_count_table(isubj)=sum(nimg_meta_all.PID_num(nimg_patch_list)==isubj);
        end
        nimg_meta.npatches=length(nimg_patch_list);
        nimg_meta.npatches_all=nimg_meta_all.npatches;
        nimg_meta.PID_num=nimg_meta_all.PID_num(nimg_patch_list);
        nimg_meta.XYLen=nimg_meta_all.XYLen(nimg_patch_list,:);
        nimg_meta.XYPos=nimg_meta_all.XYPos(nimg_patch_list,:);
        nimg_meta.View_num=nimg_meta_all.View_num(nimg_patch_list);
        nimg_meta.Examp_num=nimg_meta_all.Examp_num(nimg_patch_list);
        nimg_meta.patch_list=nimg_patch_list;
        nimg_meta.nimg_sel_used=nimg_sel_used;
        %
        for imptr=1:length(nimg_patch_list)
            isp=nimg_patch_list(imptr);
            nimg_meta.PID{imptr,1}=nimg_meta_all.PID{isp};
            nimg_meta.View{imptr,1}=nimg_meta_all.View{isp};
        end
        disp(sprintf(' Check:       data read (%3.0f of %3.0f galleries, %3.0f of %3.0f files): %8.0f patches of size %5.0f x %5.0f.',...
            length(unique(nimg_meta.Gal_num(nimg_meta.PID_num))),length(nimg_meta.gal_names_sorted),...
            sum(nimg_meta.patch_count_table>0),nimg_meta.nsubjs,nimg_meta.npatches,prod(nimg_meta.NR_orig),prod(nimg_meta.NR_orig)));
        %preprocess: log transform if requested
        [nimg_patches,nimg_xform_label]=ffdm_btc_calc_linlog(nimg_patches,data_sources{data_source});
        %
        nimg_focus_label=sprintf('NRorig=[%2.0f %3.0f]',nimg_meta_all.NR_orig);
        datasource_string=sprintf('Penn Natural Images %s %s',nimg_xform_label,nimg_focus_label);
%        preprocess_downsample=1;% set by default, so removed 27Jun22
        view_names=nimg_meta.views;
        file_nos=[1:nimg_meta_all.nsubjs];
        view_nos=unique(nimg_meta.View_num);
        count_subjs=size(nimg_meta.patch_count_table,1); %this is a count of all files, including files that don't have a patch that is used -- each "subject" is an image
        count_views=length(view_nos);
        roi_count_sdb=nimg_meta.patch_count_table;
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch size, minimum for all patches for this view and this subject, should all be identical but NaN if no view for that subject
        for isubj=1:count_subjs
            for iview=1:count_views
                roi_sel=intersect(find(nimg_meta.PID_num==isubj),find(nimg_meta.View_num==iview));
                roi_count_sdb(isubj,iview)=length(roi_sel);
                roi_list_sdb{isubj,iview}=roi_sel;
                if length(roi_sel)>0
                    roi_sizes(isubj,iview,:)=reshape(min(nimg_meta.XYLen(roi_sel,:),[],1),[1 1 2])/preprocess_downsample; %all rois of same size, take downsampling into account
                end
            end
        end
        %
    case 9 %border-parsed patches
        %
        %read the metadata
        %
        nibs_meta_path=strrep(nibs_meta_path_def,'\','/');
        nibs_meta_file=strrep(nibs_meta_file_def,'\','/');
        nibs_image_path=strrep(nibs_image_path_def,'\','/');
        ifok=0;
        while (ifok==0)
            nibs_meta_path=getinp('path to nibs metadata','s',[],nibs_meta_path);
            nibs_meta_file=getinp('file name for nibs metadata','s',[],nibs_meta_file);
            nibs_meta_fullname=cat(2,nibs_meta_path,filesep,nibs_meta_file);
            nibs_image_path=getinp('path to nibs images','s',[],nibs_image_path);
            ifok=1;
            if ~exist(nibs_meta_fullname,'file')
                disp(sprintf('metadata file %s not found.',strrep(nibs_meta_fullname,'\','/')));
                ifok=0;
            end
            if ~exist(nibs_image_path,'file')
                disp(sprintf('image path %s not found.',strrep(nibs_image_path,'\','/')));
                ifok=0;
            end
            ifok=getinp('1 if ok to read metadata','d',[0 1],ifok);
        end
        nibs_bpdata=getfield(load(nibs_meta_fullname),'bpdata');
        disp(sprintf('bpdata loaded from %s',nibs_meta_fullname));
        nibs_bpdata=filldefault(nibs_bpdata,'mask_dict',opts_nibs.mask_dict_bpdata);
        disp(sprintf('summary of %s',strrep(nibs_meta_fullname,'\','/')));
        disp(nibs_bpdata);
        disp(sprintf('typical pixels per patch: %5.0f',nibs_bpdata.patchProperties.pixelsPerPatch));
        %
        %select an aspect ratio and smoothness parameter and retrieve
        %corresponding boundary patch metadata
        %
        ifok=0;
        while (ifok==0)
            nibs_aspects=length(nibs_bpdata.patchProperties.heightWidthRatio);
            for i_aspect=1:nibs_aspects
                disp(sprintf('%2.0f->aspect ratio %5.2f',i_aspect,nibs_bpdata.patchProperties.heightWidthRatio(i_aspect)));
            end
            nibs_aspect_ptr=getinp('choice','d',[1 nibs_aspects]);
            %
            nibs_smooths=length(nibs_bpdata.algoParams.smoothness);
            for i_smooth=1:nibs_smooths
                disp(sprintf('%2.0f->smoothness param value %5.2f',i_smooth,nibs_bpdata.algoParams.smoothness(i_smooth)));
            end
            nibs_smoothness_ptr=getinp('choice','d',[1 nibs_smooths]);
            %    
            opts_nibs_show.if_objlist=getinp('1 to show object list for each image','d',[0 1],0);
            opts_nibs_show.aspect_ptrs=nibs_aspect_ptr;
            opts_nibs_show.smoothness_ptrs=nibs_smoothness_ptr;
            [opts_nibs_show_used,nibs_warnings]=nibs_bpdata_show(nibs_bpdata,opts_nibs_show);
            disp(sprintf('warnings encountered: %4.0f',size(nibs_warnings,1)));
            if (length(nibs_warnings)>0)
                disp('first and last warnings:')
                disp(nibs_warnings(1,:));
                disp(nibs_warnings(end,:));
            end
            ifok=getinp('1 if ok to retrieve patches (-1 to retain bpdata for all aspect ratios and smoothness)','d',[-1 1]);
            if (ifok==1)
                nibs_bpdata=nibs_bpdata_prune(nibs_bpdata,{'patches'},nibs_aspect_ptr,nibs_smoothness_ptr);
            end
        end
        nibs_NRpatch_maxposs=min(data_sources{data_source}.imsize);
        nibs_NRpatch_log2maxposs=floor(log(nibs_NRpatch_maxposs)/log(2));
        nibs_NRpatch_log2def=min(nibs_NRpatch_log2def,nibs_NRpatch_log2maxposs);
        nibs_NRpatch_log2max=getinp('log(2) NRmax (size of region surrounding figure and ground used for whitening)','d',...
            [5 nibs_NRpatch_log2maxposs],nibs_NRpatch_log2def);
        nibs_NRpatch_max=2.^nibs_NRpatch_log2max;
        %
        opts_nibs_patchpair.rect_size=nibs_NRpatch_max;
        opts_nibs_patchpair.pathImagesDir=nibs_image_path;
        opts_nibs_patchpair.if_showimages=0;
        opts_nibs_patchpair.if_noread=0;
        %
        %get the metadata for the chosen aspect and smoothness params, and some other fields from bpdata
        nibs_meta_use=nibs_get_meta(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr);
        %preprocess_downsample=1; % set by default
        %
        %PID is a unique designator for each "subject" (was a "patient ID" in ffdm)
        %each "subject" corresponds to a single image file
        %but there can be multiple examples (patch pairs) per subject
        %above is as in nimg.  But in ffdm, each patient subject has multiple files, one for each view.
        %
        view_names=nibs_meta_use.views;
        file_nos=[1:nibs_meta_use.nsubjs];
        view_nos=unique(nibs_meta_use.View_num);
        count_subjs=nibs_meta_use.nsubjs; %this is a count of all files, including files that don't have a patch that is used -- each "subject" is an image
        count_views=length(view_nos);
        roi_count_sdb=zeros(count_subjs,count_views); %how many pat6chpairs for each subject
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch size, minimum for all patches for this view and this subject, should all be identical but NaN if no view for that subject
        nibs_patches=zeros(nibs_NRpatch_max,nibs_NRpatch_max,nibs_meta_use.npatches);
        %nibs_composite_masks=zeros(nibs_NRpatch_max,nibs_NRpatch_max,nibs_meta_use.npatches);
        nibs_patches_img_size=zeros(nibs_meta_use.npatches,2);
        %retrieve all patch pairs and composite mask
        for patch_pair_ID=1:nibs_meta_use.npatches
            isubj=nibs_meta_use.PID_num(patch_pair_ID);
            iview=nibs_meta_use.View_num(patch_pair_ID); %always 1
            roi_count_sdb(isubj,iview)=roi_count_sdb(isubj,iview)+1;
            roi_list_sdb{isubj,iview}(end+1,1)=patch_pair_ID;
            [nibs_patchpair,hfig,opts_nibs_patchpair_used]=nibs_get_patchpair(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,opts_nibs_patchpair);
            nibs_patches_img_size(patch_pair_ID,:)=opts_nibs_patchpair_used.img_size; %save image size for later retrieval of masks
            nibs_patches(:,:,patch_pair_ID)=nibs_patchpair.rect_image;
            %nibs_composite_masks(:,:,patch_pair_ID)=nibs_patchpair.rect_mask.composite;
            roi_sizes(isubj,iview,:)=nibs_NRpatch_max;
        end
        disp(sprintf(' %5.0f patch pairs and composite masks retrieved',nibs_meta_use.npatches));
        ifok=0; %show some example patches if requested
        while (ifok==0)
            patch_pair_ID=getinp('patch pair to show (0 to end)','d',[0 nibs_meta_use.npatches]);
            if (patch_pair_ID==0)
                ifok=1;
            else
                [nibs_patchpair,hfig,opts_nibs_patchpair_show]=nibs_get_patchpair(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,...
                    setfield(opts_nibs_patchpair,'if_showimages',1));
            end
        end
        ifok=0; %show some stats on masks if requested
        while (ifok==0)
            patch_pair_ID=getinp('patch pair ID to test retrieval of masks (0 to end)','d',[0 nibs_meta_use.npatches]);
            if (patch_pair_ID==0)
                ifok=1;
            else
                [nibs_patchpair,hfig,opts_nibs_patchpair_retrieve]=nibs_get_patchpair(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,...
                    setfields(opts_nibs_patchpair,{'if_noread','img_size'},{1,nibs_patches_img_size(patch_pair_ID,:)}));
                %disp(nibs_patchpair.rect_mask);
                nibs_rect_mask_fields=fieldnames(nibs_patchpair.rect_mask);
                %verification with nibs_get_mask
                rect_mask_composite=nibs_get_mask(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,nibs_patches_img_size(patch_pair_ID,:),nibs_NRpatch_max);
                for im=1:length(nibs_rect_mask_fields)
                    nibs_mask_name=nibs_rect_mask_fields{im};
                    nibs_mask=nibs_patchpair.rect_mask.(nibs_mask_name);
                    if isfield(opts_nibs.mask_dict_composite,nibs_mask_name) & ~strcmp(nibs_mask_name,'composite')
                        nibs_composite_key=opts_nibs.mask_dict_composite.(nibs_mask_name);
                        nibs_composite_count=sum(rect_mask_composite(:)==nibs_composite_key);
                        nibs_composite_string=sprintf(' via get_mask (key=%2.0f): %4.0f',nibs_composite_key,nibs_composite_count);
                    else
                        nibs_composite_string=[];
                    end
                    if all(size(nibs_mask)==nibs_NRpatch_max) %for each mask, show how many entries of each value
                        nibs_mask_string=sprintf('%30s',nibs_mask_name);
                        nu=unique(nibs_mask(:));
                        for inu=1:length(nu)
                            nibs_mask_string=cat(2,nibs_mask_string,sprintf(' %2.0f->%4.0f',nu(inu),sum(nibs_mask(:)==nu(inu))));
                        end
                        disp(cat(2,nibs_mask_string,nibs_composite_string));                           
                    end %mask?
                end %fields
            end
        end
        %ask about region types
        disp('options for region selection');
        nibs_mask_types=fieldnames(opts_nibs.mask_dict_bpdata);
        nibs_mask_type_vals=[];
        for im=1:length(nibs_mask_types)
            nibs_mask_type_vals=[nibs_mask_type_vals,opts_nibs.mask_dict_bpdata.(nibs_mask_types{im})];
            disp(sprintf('%1.0f->%s',im-1,nibs_mask_types{im}));
        end
        ifok=0;
        while (ifok==0)
            regtype_count=getinp('number of region types','d',[1 length(nibs_mask_types)+2]); %+2 to also allow for all, or repeats
            regtype_masks=cell(regtype_count,1);
            regtype_descs=cell(regtype_count,1);
            for iregtype=1:regtype_count
                regtype_descs{iregtype}=[];
                regtype_masks{iregtype}=getinp(sprintf('choice(s) for region type %1.0f (-1: use all of patch)',iregtype),'d',[-1 max(nibs_mask_type_vals)]);
                if any(regtype_masks{iregtype}==-1)
                    regtype_masks{iregtype}=nibs_mask_type_vals;
                    regtype_descs{iregtype}='all';
                else
                    regtype_masks{iregtype}=intersect([regtype_masks{iregtype}],nibs_mask_type_vals);
                    for im=1:length(regtype_masks{iregtype})
                        regtype_descs{iregtype}=cat(2,regtype_descs{iregtype},nibs_mask_types{find(nibs_mask_type_vals==regtype_masks{iregtype}(im))},'+');
                    end
                    regtype_descs{iregtype}=regtype_descs{iregtype}(1:end-1);
                end
            end
            disp('chosen region selections:')
            for iregtype=1:regtype_count
                disp(sprintf('region type %2.0f: %s',iregtype,regtype_descs{iregtype}));
            end
            regtype_binarize_useall=getinp('1 to use the entire patch for binarizing, 0 to use just selected area','d',[0 1],regtype_binarize_useall);
            ifok=getinp('1 if ok','d',[0 1]);
        end
        %ask about log transform and do it
        [nibs_patches,nibs_xform_label]=ffdm_btc_calc_linlog(nibs_patches,data_sources{data_source});
        %
        datasource_string=sprintf('%s %s asp rat: %4.2f, smooth: %4.2f, pxls: %4.0f, NRmax: %6.0f',...
            data_sources{data_source}.label,nibs_xform_label,...
            nibs_bpdata.patchProperties.heightWidthRatio(nibs_aspect_ptr),...
            nibs_bpdata.algoParams.smoothness(nibs_smoothness_ptr),...
            nibs_bpdata.patchProperties.pixelsPerPatch,...
            nibs_NRpatch_max);
    case 10 %bxrs: bone xray, scintigrams
        %the patches in ffdm are the rectangular ROIs in bone_btc_demo;
        if ~exist('bxrs_path') bxrs_path='..\..\Teaching\Ajayi\bone\'; end
        [bxrs_file,bxrs_desc,opts_bxrs]=bone_scint_select();
        %
        opts_bxrs=filldefault(opts_bxrs,'if_log',0);
        opts_bxrs.file_path=bxrs_path;
        opts_bxrs.im_path=bxrs_path;
        opts_bxrs.roi_path=bxrs_path;
        %
        [bxrs_im_data_all,bxrs_roi_data_all,bxrs_roi_sizes,bxrs_roi_origimg,bxrs_im_files,opts_bxrs_used]=bone_read_xls(bxrs_file,opts_bxrs);
        bxrs_n_entries=length(bxrs_im_data_all);
        bxrs_n_rois=length(bxrs_roi_data_all);
        bxrs_minsize=min(bxrs_roi_sizes(:));
        bxrs_minsize_pwr2=2^floor(log2(bxrs_minsize));
        disp(sprintf('Read bxrs data: %3.0f entries with %5.0f separate ROI files',bxrs_n_entries,bxrs_n_rois));
        disp(sprintf('minimum roi size, minimum power-of-2: %5.0f %5.0f',bxrs_minsize,bxrs_minsize_pwr2));
        datasource_string=strrep(strrep(bxrs_file,'.xlsx',''),'.xls','');
        preprocess_downsample=1;
        %
        bxrs_meta=struct();
        bxrs_meta.bxrs_path=bxrs_path;
        bxrs_meta.bxrs_file=bxrs_file;
        bxrs_meta.npatches=bxrs_n_rois;
        bxrs_meta.views={'std'};
        bxrs_meta.PID_num=bxrs_roi_origimg(:);
        bxrs_meta.XYLen=bxrs_roi_sizes;
        bxrs_meta.XYPos=ones(bxrs_n_rois,2);
        bxrs_meta.View_num=ones(bxrs_n_rois,1);
        %
        bxrs_meta.PID=cell(bxrs_n_rois,1);
        bxrs_meta.Examp_num=zeros(bxrs_n_rois,1);
        bxrs_meta.patch_count_table=zeros(bxrs_n_entries,1);
        for ib=1:bxrs_n_rois
            bxrs_meta.PID{ib}=bxrs_im_files{bxrs_roi_origimg(ib)};
            bxrs_meta.Examp_num(ib)=sum(bxrs_roi_origimg(1:ib)==bxrs_roi_origimg(ib));
            bxrs_meta.patch_count_table(bxrs_roi_origimg(ib))=1+bxrs_meta.patch_count_table(bxrs_roi_origimg(ib));
        end
        %
        %may want to ask about excluding subpatches with no change, and usethis later
        %
        view_names=bxrs_meta.views;
        %file_nos=bxrs_meta.PID_num;
        file_nos=[1:bxrs_n_entries]'; %file_nos is indexed by subject ID, and there is one entry per subject ID
        view_nos=unique(bxrs_meta.View_num);
        count_subjs=size(bxrs_meta.patch_count_table,1);
        count_views=length(view_nos);
        %
        roi_count_sdb=bxrs_meta.patch_count_table;
        roi_list_sdb=cell(count_subjs,count_views); %list of patch numbers
        roi_sizes=NaN(count_subjs,count_views,2); %patch size, minimum for all patches for this view and this subject, NaN if no view for that subject
        roi_sizes_all=cell(count_subjs,count_views);
        for isubj=1:count_subjs
            for iview=1:count_views
                roi_sel=intersect(find(bxrs_meta.PID_num==isubj),find(bxrs_meta.View_num==iview));
                roi_count_sdb(isubj,iview)=length(roi_sel);
                roi_list_sdb{isubj,iview}=roi_sel;
                if length(roi_sel)>0
                    roi_sizes(isubj,iview,:)=reshape(min(bxrs_meta.XYLen(roi_sel,:),[],1),[1 1 2]);
                    roi_sizes_all{isubj,iview}=bxrs_meta.XYLen(roi_sel,:);
                end
            end
        end
        %
        if if_clear
            clear bxrs_im_data_all
        end
end %data_source
if ~exist('roi_sizes_all') %added 16Aug: roi_sizes_all separately calculated for data_source=10 since each subject may have several rois of different sizes
    roi_sizes_all=cell(count_subjs,count_views);
    for isubj=1:count_subjs
        for iview=1:count_views
            roi_sizes_all{isubj,iview}=reshape(roi_sizes(isubj,iview,:),[1 2]);
        end
    end
end
disp(sprintf('loaded %s with %2.0f subjects (files %2.0f to %2.0f) and %2.0f views (%2.0f->%s to %2.0f->%s)',...
    datasource_string,...
    count_subjs,min(file_nos),max(file_nos),...
    count_views,min(view_nos),view_names{min(view_nos)},max(view_nos),view_names{max(view_nos)}));
disp('all view names')
disp(view_names)
disp('Options for determining whitening spectra');
for iavg_type=1:whiten_avg_types
    disp(sprintf('%1.0f->%s',iavg_type,whiten_avg_types_labels_orig{iavg_type}));
end
whiten_avg_type_use=getinp('choice','d',[1 whiten_avg_types],whiten_avg_type_use);
if strcmp(whiten_avg_types_labels{whiten_avg_type_use},'power law')
    whiten_power_law=getinp('power law slope (-2 for "1/f^2" natural images)','f',[-10 10]);
    whiten_avg_types_labels{whiten_avg_type_use}=sprintf('power law %4.1f',whiten_power_law);
end
%
fft2_padfactor=getinp('pad factor','d',[1 8],fft2_padfactor);
if_spec_subpxlmean=getinp('1 to subtract pixel mean','d',[0 1],if_spec_subpxlmean);
if_spec_subroimean=getinp('1 to subtract ROI mean','d',[0 1],if_spec_subroimean);
if_spec_cosbell=getinp('1 or 2 for cosine bell window (1: den=R-1, 2: den=R)','d',[0 2],if_spec_cosbell);
%
roi_minsize=[min(min(roi_sizes(:,:,1))),min(min(roi_sizes(:,:,2)))];
disp(sprintf('Minimum roi size in files: %4.0f x %4.0f (i.e., all subjects and all views have an ROI this size or larger)',roi_minsize));
% override? (added 16Aug22 for data_source==10)
NR_max_allowed_log2=12;
NR_max_poss_log2=floor(log(min(roi_minsize))/log(2));
if (data_source==10)
    NR_max_poss_log2=getinp('minimum roi size to require (log 2(NR))','d',[NR_max_poss_log2 NR_max_allowed_log2],NR_max_poss_log2);
end
NR_max_poss=2.^NR_max_poss_log2;
N_list_log2_def=[1:NR_max_poss_log2-R_min_log2];
if_NRS_ok=0;
Sstep_char={' ','+'}; %+ will indicate that phases are stepped
while (if_NRS_ok==0)
    N_list_log2=getinp('values of blocking (log 2(N))','d',[0 max(N_list_log2_def)],[0:max(N_list_log2_def)]);
    N_list=2.^N_list_log2;
    NR_log2_list=getinp('values of ROI size in pixels (log 2(NR))','d',[N_list_log2(1)+2 NR_max_poss_log2],repmat(NR_max_poss_log2,1,length(N_list)));
    NR_list=2.^NR_log2_list;
    if (data_source==1)
        NR_max_log2=max(NR_log2_list);
    else
        NR_max_log2=NR_max_poss_log2; %since NR_max determines the size of the roi to be read, and these are fixed by the sdb database
    end
    NR_max=2.^NR_max_log2;
    R_list=NR_list./N_list;
    %
    % 10 Aug 20:  S is the scale for calculating image statistics (after downsampling by N)
    %
    disp('analysis parameters so far:')
    for iNR=1:length(N_list)
        fprintf(' %2.0f-> blocking (N) =%4.0f    blocks per patch edge (R) =%4.0f    patch size (NR)=%4.0f\n',...
            iNR,N_list(iNR),R_list(iNR),N_list(iNR)*R_list(iNR));
    end
    if_S_ok=0;
    while if_S_ok==0
        S_list_log2=getinp('scales for calculating image statistics (log2(S), must be less than correpsponding values of R/4)',...
            'd',[0 round(log(max(R_list))/log(2))],zeros(1,length(N_list)));
        S_list=2.^S_list_log2;
        if_S_ok=all(S_list<=(R_list/4));
    end
    Sstep_list=getinp('1 to step the starting phase for scaling','d',[0 1],zeros(1,length(N_list)));
    %
    disp(sprintf('analysis parameters including scaling, maximum total patch size (N*R)=%5.0f',NR_max));
    for iNR=1:length(N_list)
        disp(sprintf(' %2.0f-> blocking (N) =%4.0f    blocks per patch edge (R) =%4.0f    patch size (NR)=%4.0f   scale for btc statistics, in blocks(S)=%4.0f step phase: %s',...
            iNR,N_list(iNR),R_list(iNR),N_list(iNR)*R_list(iNR),S_list(iNR),Sstep_char{Sstep_list(iNR)+1}));
    end
    if_NRS_ok=getinp('1 if ok','d',[0 1]);
end %if_NRS_ok
%
N_list_log2_psa=getinp('values of blocking (log 2) for detailed spectra, -1 for none','d',[-1  max(N_list_log2)],N_list_log2);
N_list_psa=2.^N_list_log2_psa;
%
margin_width=getinp('margin size for calculating binary texture statistics','d',[0 floor(min(R_list)/2)-1],margin_width);
%
n_patches_old=floor(roi_sizes(:,:,1)/NR_max).*floor(roi_sizes(:,:,2)/NR_max);
n_patches=zeros(count_subjs,count_views);
for isubj=1:count_subjs
    for iview=1:count_views
        n_patches_eachdim=floor(roi_sizes_all{isubj,iview}/NR_max);
        if ~isempty(n_patches_eachdim)
            n_patches(isubj,iview)=sum(n_patches_eachdim(:,1).*n_patches_eachdim(:,2));
        end
    end
end
%note that n_patches may have some NaN's if roi_sizes_all has NaN's, which
%could happen if there are no images with the specified subject and view
n_patches_min=min(n_patches(:));
n_patches_max=max(n_patches(:));
disp(sprintf('There are %4.0f to %4.0f patches of size %4.0f x %4.0f pixels per image, total of %6.0f patches',...
    n_patches_min,n_patches_max,NR_max,NR_max,sum(n_patches(~isnan(n_patches(:))))));
if ~exist('patches_usefrac') patches_usefrac=1; end
patches_usefrac=getinp('1/(fraction of patches to use) (1: use all)','f',[1 n_patches_max],1); %"1/" added 15Aug22
max_patches=getinp('maximum number of big patches per image to use','d',[1 ceil(n_patches_max/patches_usefrac)],ceil(n_patches_max/patches_usefrac));
%
if ~exist('cl_prob') cl_prob=0.95; end
cl_prob=getinp('confidence limit probability (0 for none)','f',[0 0.99999],cl_prob);
%
%for data_source=1  (pilot ffdm), there is only one roi for each subject and view, and we need to divide it up into patches
%for data_source=4  (usic), there may be more than one roi for each subject (only one view), and we need to divide it up into patches
%for data_source=10 (bxrs), there may be more than one roi for each subject, from different files, and we need to divide it up into patches
%for data_source=2, 3, 5, 6, 9: there are multiple patches for each subject and view already, and no need to divide up
%
patch_whichroi=cell(count_subjs,count_views);
if data_sources{data_source}.if_makepatches %data_source = 1 or 4 or 10
    if_exclude_constant=getinp('1 to exclude patches with constant subpatches','d',[0 1],double(data_source==10));
    if (if_exclude_constant)
        NR_exconst_log2=getinp('size of subpatch (log 2)','d',[R_min_log2 NR_max_poss_log2],R_min_log2);
        NR_exconst=2.^NR_exconst_log2;
    end
    %
    %choose image patches by tesselating them into size NR_max x NR_max
    %with a random start point for each row and column of patches (given by patch_strips)
    %choosing up to max_patches systematically or randomly from them
    %
    patch_starts=cell(count_subjs,count_views); %low x and y of each patch and examp number (0 for data_source=1)
    patch_strips=cell(count_subjs,count_views,2,1); %start of each strip of possible patch positions, for each dimension, each example number
    patch_count=zeros(count_subjs,count_views);
    patch_exconst_lists=cell(count_subjs,count_views); %list of excluded patches because they are constant
    %
    patch_exconst_tot=0;
    patch_count_avail_tot=0;
    %
    for isubj=1:count_subjs
        for iview=1:count_views
            %randomly divide excess into patch_each+1 non-negative quantities
            if (data_source==1)
                nexamps=1;
            else %more than one
                nexamps=roi_count_sdb(isubj,iview);
            end
            for iexamp=1:nexamps
                patch_exconst_lists{isubj,iview}=cell(1,nexamps);
                if (data_source~=10) %same size roi for each example
                    roi_size=reshape(roi_sizes(isubj,iview,:),1,2);
                else %different roi sizes
                    roi_size=roi_sizes_all{isubj,iview}(iexamp,:);
                end
                excess=mod(roi_size,NR_max);
                patch_each=floor(roi_size/NR_max); %number of patches along each dimension
                part_rands=cell(1,2);
                for ixy=1:2
                    parts=excess(ixy)+patch_each(ixy);
                    if (if_norand==1)
                        part_segs=([1:parts]>patch_each(ixy));
                    else
                        part_segs=randperm(parts)>(patch_each(ixy));
                    end
                    part_segs_cs=cumsum(part_segs);
                    part_rands{ixy}=diff([0 part_segs_cs(part_segs==0)]);
                    patch_strips{isubj,iview,ixy,iexamp}=cumsum(part_rands{ixy})+NR_max*[0:patch_each(ixy)-1];
                end
                %patch_count_thisexamp=min(ceil(n_patches(isubj,iview)/patches_usefrac),max_patches); %logic changed 16Aug22
                patch_count_avail=floor(prod(patch_each)); %patches available in this example
                patch_count_avail_tot=patch_count_avail_tot+patch_count_avail;
                patch_avail_list=[1:patch_count_avail]; %list of all available patches
                %
                %optionally (esp. data_source==10) remove patches that are constant
                %
                if (if_exclude_constant)
                    switch data_source
                        case 1
                            orig_img=roi_imgs{isubj,iview};
                        case 4
                            imptr=roi_list_sdb{isubj,iview}(iexamp);
                            orig_img=usic_rois{imptr};
                        case 10
                            imptr=roi_list_sdb{isubj,iview}(iexamp);
                            orig_img=bxrs_roi_data_all{imptr};
                        otherwise
                            error('cannot retrieve images to exclude constant patches')
                    end
                    patch_exconst_list=[]; %patches to be excluded
%                    disp('[isubj iview iexamp ipatch patch_x patch_y size(patch_to_check) patch_strips{isubj,iview,1,iexamp}(patch_x) patch_strips{isubj,iview,2,iexamp}(patch_y)]');
                    for ipatch=1:patch_count_avail
                        patch_x=1+mod(ipatch-1,patch_each(1));
                        patch_y=1+floor((ipatch-1)/patch_each(1));
                        patch_to_check=orig_img([1:NR_max]+patch_strips{isubj,iview,1,iexamp}(patch_x),[1:NR_max]+patch_strips{isubj,iview,2,iexamp}(patch_y));
 %                       disp([isubj iview iexamp ipatch patch_x patch_y size(patch_to_check) patch_strips{isubj,iview,1,iexamp}(patch_x) patch_strips{isubj,iview,2,iexamp}(patch_y)]);
                        %code adapted from bone_btc_demo to see if any subpatch of size NR_exconst is constant
                        patch_data=permute(reshape(patch_to_check,[NR_exconst,NR_max/NR_exconst,NR_exconst,NR_max/NR_exconst]),[1 3 2 4]);
                        patch_data_min=min(min(patch_data,[],2),[],1);% minimum within each region of size NR_exconst
                        patch_data_max=max(max(patch_data,[],2),[],1);% maximum within each region of size NR_exconst
                        if any(patch_data_min(:)==patch_data_max(:)) 
                            patch_exconst_list=[patch_exconst_list,ipatch];
                        end
                    end
                    patch_exconst_lists{isubj,iview}{iexamp}=patch_exconst_list;
                    patch_avail_list=setdiff(patch_avail_list,patch_exconst_list); %remove constant patches
                    patch_exconst_tot=patch_exconst_tot+length(patch_exconst_list);
                end
                patch_count_thisexamp=min(floor(length(patch_avail_list)/patches_usefrac),max_patches); %number of patches to use
                if (if_norand==1)
                    patch_choices=patch_avail_list(1:patch_count_thisexamp); %was patch_choices=find([1:n_patches(isubj,iview)]<=patch_count_thisexamp);
                else
                    patch_choices=patch_avail_list(randperm(length(patch_avail_list))<=patch_count_thisexamp); % was patch_choices=find(randperm(n_patches(isubj,iview))<=patch_count_thisexamp);
                end
                patch_count(isubj,iview)=patch_count(isubj,iview)+patch_count_thisexamp;
                if patch_count_thisexamp~=length(patch_choices)
                    error('logic error in patch count')
                end
                if length(patch_choices)>0
                    patch_x=1+mod(patch_choices-1,patch_each(1));
                    patch_y=1+floor((patch_choices-1)/patch_each(1));
                    patch_starts{isubj,iview}=cat(1,patch_starts{isubj,iview},...
                        [patch_strips{isubj,iview,1,iexamp}(patch_x)',patch_strips{isubj,iview,2,iexamp}(patch_y)' repmat(iexamp,patch_count_thisexamp,1)]); %fourth dim (iexamp) added 16Aug22
                    patch_whichroi{isubj,iview}=[patch_whichroi{isubj,iview},repmat(iexamp,1,patch_count_thisexamp)]; %added 16Aug22 to keep track of which ROI each patch comes from
                end
            end %iexamp
        end %iview
    end %isubj
    disp(sprintf(' %6.0f big patches available in images, %6.0f excluded as constant, %6.0f big patches made (1/patches_usefrac=%3.0f, max patches per image =%5.0f)',...
        patch_count_avail_tot,patch_exconst_tot,sum(patch_count(:)),patches_usefrac,max_patches));
else %datasource=2, 3, 5, 6, 9
    patch_count=roi_count_sdb; %row: each subject, col: each view, one patch per ROI
    for isubj=1:count_subjs
        for iview=1:count_views
            patch_starts{isubj,iview}=zeros(patch_count(isubj,iview),2);
        end
    end
end
s_havedata=cell(1,count_views); %which subjects have data for each view
for iview=1:count_views
    s_havedata{iview}=find(patch_count(:,iview)>0)';
end
for iview=1:count_views
    disp(sprintf('view %1.0f (%10s): %4.0f big patches to be used in %4.0f subjects',...
        iview,view_names{iview},sum(patch_count(:,iview)),sum(patch_count(:,iview)>0)));
end
if any(sum(patch_count(:,iview))==0)
    disp('cannot proceed.');
end       
s_havedata_allv=find(sum(patch_count,2)>0)'; %which subjects have data for any view
%make a list of all patches
patch_count_sum_big=sum(patch_count(:));
patches_orig=zeros(NR_max,NR_max,patch_count_sum_big); %before downsampling
patches_subj_big=zeros(1,patch_count_sum_big);  %for grouping and averaging later
patches_view_big=zeros(1,patch_count_sum_big);  %for grouping and averaging later
patches_numb_big=zeros(1,patch_count_sum_big);  %unique patch number within roi
%
patches_subj=cell(1,length(N_list)); %subject number for each patch
patches_view=cell(1,length(N_list)); %view number for each patch
patches_numb=cell(1,length(N_list)); %serial number, beginning with 1, to distinguish patches with same subject number and view number
%
patch_mult_side=(NR_max./NR_list);
patch_mult=patch_mult_side.^2; %number of patches that a "big" patch (maximum NR) can be broken into for smaller NR's
%
for iNR=1:length(N_list)
    patches_subj{iNR}=zeros(1,patch_count_sum_big*patch_mult(iNR));
    patches_view{iNR}=zeros(1,patch_count_sum_big*patch_mult(iNR));
    patches_numb{iNR}=zeros(1,patch_count_sum_big*patch_mult(iNR));
end
index=0;
for isubj=1:count_subjs
    for iview=1:count_views
        patches_subj_big(index+[1:patch_count(isubj,iview)])=isubj;
        patches_view_big(index+[1:patch_count(isubj,iview)])=iview;
        patches_numb_big(index+[1:patch_count(isubj,iview)])=[1:patch_count(isubj,iview)];
        for iNR=1:length(N_list)
            patches_subj{iNR}(index*patch_mult(iNR)+[1:patch_count(isubj,iview)*patch_mult(iNR)])=isubj;
            patches_view{iNR}(index*patch_mult(iNR)+[1:patch_count(isubj,iview)*patch_mult(iNR)])=iview;
            patches_numb{iNR}(index*patch_mult(iNR)+[1:patch_count(isubj,iview)*patch_mult(iNR)])=[1:patch_count(isubj,iview)*patch_mult(iNR)];
        end
        for ipatch=1:patch_count(isubj,iview)
            index=index+1;
            switch data_source
                case 1
                    patches_orig(:,:,index)=roi_imgs{isubj,iview}([1:NR_max]+patch_starts{isubj,iview}(ipatch,1),[1:NR_max]+patch_starts{isubj,iview}(ipatch,2));
                case 2
                    %preproc=[] means downsample=1, full patch read
                    patches_orig(:,:,index)=ffdm_read_subpatch(roi_list_sdb{isubj,iview}(ipatch),patch_metadata,patch_path,[]);
                case 3 %code modified from acmd_mlis_load
                    stim_filename=sdba.subpatch_filename{roi_list_sdb{isubj,iview}(ipatch)};
                    [stim_filename_aug,stim_path_file]=mlis_alglib_filename_aug(stim_filename,sdba,subinfo);
%                     if (length(subinfo.numeric)>0)
%                         stim_number=str2num(stim_filename(2:min(find(stim_filename=='_'))-1));
%                         subdirectory_num=sum(double(subinfo.numeric<=stim_number));
%                         subdirectory_name=deblank(subinfo.names(subdirectory_num,:));
%                         stim_filename_aug=cat(2,subdirectory_name,filesep,stim_filename);
%                         stim_path_file=cat(2,sdba.stimulus_path,filesep,stim_filename_aug);
%                     else
%                         stim_path_file=cat(2,sdba.stimulus_path,filesep,stim_filename);
%                     end
                    patches_orig(:,:,index)=imread(stim_path_file,'bmp');
                case 4
                    iexamp=patch_starts{isubj,iview}(ipatch,3);
                    imptr=roi_list_sdb{isubj,iview}(iexamp);
                    patches_orig(:,:,index)=usic_rois{imptr}([1:NR_max]+patch_starts{isubj,iview}(ipatch,1),[1:NR_max]+patch_starts{isubj,iview}(ipatch,2));
                case 5
                    imptr=roi_list_sdb{isubj,iview}(ipatch);
                    patches_orig(:,:,index)=mrix_subpatches{imptr}([1:NR_max],[1:NR_max]);
                case 6
                    imptr=roi_list_sdb{isubj,iview}(ipatch);
                    patches_orig(:,:,index)=nimg_patches([1:NR_max],[1:NR_max],imptr);
                case 9
                    imptr=roi_list_sdb{isubj,iview}(ipatch);
                    patches_orig(:,:,index)=nibs_patches([1:NR_max],[1:NR_max],imptr);
                case 10
                    iexamp=patch_starts{isubj,iview}(ipatch,3);
                    imptr=roi_list_sdb{isubj,iview}(iexamp);
%                     if (isubj==1) | isubj==count_subjs
%                         disp('isubj iexamp imptr ipatch patch_starts');
%                         disp([isubj iexamp imptr ipatch patch_starts{isubj,iview}(ipatch,:)]);
%                     end
                    patches_orig(:,:,index)=bxrs_roi_data_all{imptr}([1:NR_max]+patch_starts{isubj,iview}(ipatch,1),[1:NR_max]+patch_starts{isubj,iview}(ipatch,2));
            end
        end
    end
end
disp(sprintf('Extracted %4.0f big patches of size %4.0f x %4.0f',patch_count_sum_big,NR_max,NR_max));
%
%show original roi and patches and then start processing
%
patches_downsampled=cell(1,length(N_list));
ps_patches=cell(1,length(N_list));
ps_avg_v=cell(1,length(N_list));
ps_avg_vs=cell(1,length(N_list)); %power, per view and per subject
patches_whitened=cell(1,length(N_list));
patches_binarized=cell(1,length(N_list));
patches_blockcounts=cell(1,length(N_list));
%
%examples to show
if data_sources{data_source}.if_makepatches
    % example 1 is first patch in first view of first subject
    % last example is last patch in last view of last subject
    % if there are 2 or more views,then other examples are intermediate views of intermediate subjects
    if (count_views==1)
        examp_subj=[1 count_subjs];
        examp_view=[1 count_views];
        examp_numb=[1 patch_count(count_subjs,count_views)];
        examp_roi=[1 1];
        if (data_source==10) %special logic since some subjects may not have any patches (if NR_max is set larger than ROI size)
            examp_subj=[min(find(sum(patch_count,2)>0)) max(find(sum(patch_count,2)>0))];
            examp_numb=[1 patch_count(examp_subj(end),count_views)];
            examp_roi=[min(patch_whichroi{1,1}) max(patch_whichroi{examp_subj(end),count_views})];
        end
    else
        examp_subj=1+round((count_subjs-1)*[0:count_views-1]/(count_views-1));
        examp_view=[1:count_views];
        for iexamp=1:length(examp_subj)
            examp_numb(iexamp)=1+round((patch_count(examp_subj(iexamp),examp_view(iexamp))-1)*(iexamp-1)/(count_views-1));
        end
    end
    count_examps=length(examp_subj);
else %take examples equally spaced in list of patches
    count_examps=max(2,count_views);
    examp_indices=1+round((patch_count_sum_big-1)*[0:count_examps-1]/(count_examps-1));
    examp_subj=patches_subj_big(examp_indices);
    examp_view=patches_view_big(examp_indices);
    examp_numb=patches_numb_big(examp_indices);
    examp_roi=ones(1,count_examps);
end
roi_img_examp=cell(1,count_examps);
if ~exist('NR_exclude') %added 15Aug22
    NR_exclude(examp_subj(:))=examp_subj(:);
end
for iexamp=1:count_examps
    if (data_source==1) 
        roi_img_examp{iexamp}=roi_imgs{examp_subj(iexamp),examp_view(iexamp)};
    elseif (data_source==4)
        roi_img_examp{iexamp}=usic_rois{roi_list_sdb{examp_subj(iexamp),examp_view(iexamp)}(examp_numb(iexamp))};
    elseif (data_source==10)
        roi_img_examp{iexamp}=bxrs_roi_data_all{roi_list_sdb{examp_subj(iexamp),examp_view(iexamp)}(examp_roi(iexamp))};
    else %data from patches or stimuli
        %find index
        subj_sel=find(patches_subj_big==examp_subj(iexamp));
        view_sel=find(patches_view_big==examp_view(iexamp));
        numb_sel=find(patches_numb_big==examp_numb(iexamp));
        examp_index=intersect(intersect(subj_sel,view_sel),numb_sel);
        if length(examp_index)~=1
            subj_sel
            view_sel
            numb_sel
            examp_index
            error('Cannot find example index.  Above *_sel quantities should have a unique intersection.');
        else
            disp(sprintf(' example is index %5.0f, taken from subj ptr %3.0f (orig subj %3.0f), view ptr %2.0f (orig view %2.0f), patch %3.0f',...
                examp_index,examp_subj(iexamp),...
                NR_exclude(examp_subj(iexamp)),examp_view(iexamp),....
                view_nos(examp_view(iexamp)),examp_numb(iexamp)));
            roi_img_examp{iexamp}=patches_orig(:,:,examp_index);
        end
    end
    disp(sprintf(' pipeline example %2.0f: isubj %3.0f iview %3.0f examp numb (patch) %3.0f roi %3.0f',...
        iexamp,examp_subj(iexamp),examp_view(iexamp),examp_numb(iexamp),examp_roi(iexamp)));
end
nrows_patch=count_examps;
if (if_clear)
    clear roi_imgs;
    clear mrix_subpatches usic_rois nimg_patches nibp_patches nibp_data;
end
%
%loop to show pipeline, downsample, calculate spectra, whiten, and calculate block counts
%
for iNR=1:length(N_list)
    N=N_list(iNR);
    NR=NR_list(iNR);
    R=NR/N;
    S=S_list(iNR);
    Sstep=Sstep_list(iNR);
    Splus=Sstep_char{Sstep_list(iNR)+1};
    %
    tstring=sprintf('%s N=%1.0f R=%1.0f S=%1.0f%s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f  gen usefrac=%1.0f',...
        datasource_string,N,R,S,Splus,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
        if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell,patches_usefrac);
    if (if_norand==1)
        tstring=cat(2,tstring,' NO randomization');
    end
    if (if_norand==-1)
        tstring=cat(2,tstring,' Frozen randomization');
    end
    if (if_plot_pipeline)
        %
        %show first and last samples
        %
        hf_patch=figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'pipeline: ',tstring));
        set(gcf,'Position',[50 50 1200 750]);
        ncols_patch=5; %roi, selected patch, selected patch downsampled, whitened, binarized
        for irow=1:nrows_patch
            isubj=examp_subj(irow);
            iview=examp_view(irow);
            ipatch_show=examp_numb(irow);
            patch_index_show_big(irow)=find((patches_subj_big==isubj) & (patches_view_big==iview) & (patches_numb_big==ipatch_show));
            patch_index_show(irow)=1+(patch_index_show_big(irow)-1)*patch_mult(iNR);
            %show the ROI with bigpatch selection
            icol=1;
            subplot(nrows_patch,ncols_patch,icol+(irow-1)*ncols_patch);
            imagesc(roi_img_examp{irow},roi_range);
            colormap gray;
            axis equal;
            axis tight;
            set(gca,'XTick',[1 roi_sizes(isubj,iview,2)]); %because imagesc plots the transpose
            set(gca,'YTick',[1 roi_sizes(isubj,iview,1)]);
            title(sprintf('subj %1.0f v%1.0f->%s',isubj,iview,view_names{iview}));
            hold on;
            if isempty(patch_whichroi{isubj,iview}) %logic added 16Aug22 to only show the patches from the chosen roi
                patch_showlist=[1:patch_count(isubj,iview)];
            else
                patch_showlist=find(patch_whichroi{isubj,iview}==examp_roi(irow));
            end
            for ipatch=patch_showlist
                patch_start=patch_starts{isubj,iview}(ipatch,:);
                if (ipatch==ipatch_show)
                    if patch_mult_side(iNR)>1 %draw subdivisions
                        for isub=1:patch_mult_side(iNR)-1
                            plot(repmat(patch_start(2)+NR_max*(isub/patch_mult_side(iNR)),1,2),patch_start(1)+[0 NR_max],'c');
                            plot(patch_start(2)+[0 NR_max],repmat(patch_start(1)+NR_max*(isub/patch_mult_side(iNR)),1,2),'c');
                        end
                        plot(repmat(patch_start(2)+NR_max/patch_mult_side(iNR),1,2),patch_start(1)+[0 NR_max/patch_mult_side(iNR)],'m');
                        plot(patch_start(2)+[0 NR_max/patch_mult_side(iNR)],repmat(patch_start(1)+NR_max/patch_mult_side(iNR),1,2),'m');
                    end
                    color='g';
                else
                    color='y';
                end
                plot(patch_start(2)-0.5+[1 NR_max NR_max 1 1],patch_start(1)-0.5+[1 1 NR_max NR_max 1],color);
            end
            %show one big patch
            icol=2;
            subplot(nrows_patch,ncols_patch,icol+(irow-1)*ncols_patch);
            imagesc(patches_orig(:,:,patch_index_show_big(irow)),roi_range);
            colormap gray;
            axis square;
            axis tight;
            set(gca,'XTick',[1 NR_max]);
            set(gca,'YTick',[1 NR_max]);
            hold on;
            if patch_mult_side(iNR)>1 %draw subdivisions
                for isub=1:patch_mult_side(iNR)-1
                    plot(repmat(NR_max*(isub/patch_mult_side(iNR)),1,2),[0 NR_max],'c');
                    plot([0 NR_max],repmat(NR_max*(isub/patch_mult_side(iNR)),1,2),'c');
                end                   
                plot(repmat(NR_max/patch_mult_side(iNR),1,2),[0 NR_max/patch_mult_side(iNR)],'m');
                plot([0 NR_max/patch_mult_side(iNR)],repmat(NR_max/patch_mult_side(iNR),1,2),'m');
            end
            title(sprintf('patch %1.0f of %1.0f',ipatch_show,patch_count(isubj,iview)));
            hold on;
        end
        %
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end %if_plot_pipeline
    %
    %downsample
    %
    disp(sprintf('downsampling for N=%3.0f R=%3.0f',N,R));
    patch_count_sum=patch_count_sum_big*patch_mult(iNR);
    patches_downsampled{iNR}=zeros(R,R,patch_count_sum);
    %patches_downsampled{iNR}=reshape(mean(mean(reshape(patches_orig,[N R N R patch_count_sum]),3),1),[R R patch_count_sum]);
    Z=mean(mean(reshape(patches_orig,[N R patch_mult_side(iNR) N R patch_mult_side(iNR) patch_count_sum_big]),4),1); %block-average 
    Zp=permute(Z,[2 5 3 6 7 1 4]); %now dimensions are Xblock Yblock Xsubpatch Ysubpatch bigpatch 1 1 
    patches_downsampled{iNR}=reshape(Zp,[R R patch_count_sum]);
    clear Z Zp;
    if (if_plot_pipeline)
        %show one patch
        icol=3;
        for irow=1:nrows_patch
            subplot(nrows_patch,ncols_patch,icol+(irow-1)*ncols_patch);
            imagesc(patches_downsampled{iNR}(:,:,patch_index_show(irow)),roi_range);
            colormap gray;
            axis square;
            axis tight;
            set(gca,'XTick',[1 R]);
            set(gca,'YTick',[1 R]);
            title(sprintf('downsampled, N=%1.0f',N));
            hold on;
        end %irow
    end %if_plot_pipeline
    %
    %compute 2dFFT's in each patch, both for spectrum and to do whitening
    %
    disp(sprintf('doing fft''s for %4.0f patches of size %4.0f x %4.0f with pad factor %1.0f',patch_count_sum,R,R,fft2_padfactor));
    fft_length=fft2_padfactor*R;
    fft_half=fft_length/2;
    ps_patches{iNR}=zeros(fft_length,fft_length,patch_count_sum);
    sfreqs_per_pixel=min([0:fft_length-1],fliplr([1:fft_length]))/fft_length/N;
    ps_avg_v{iNR}=zeros(fft_length,fft_length,count_views,ps_avg_types);
    ps_avg_vs{iNR}=zeros(fft_length,fft_length,count_views,count_subjs);
    patches_fft=zeros(fft_length,fft_length,patch_count_sum);
%
%     cosbell=ones(1,R);
%     if (if_spec_cosbell==1)
%         cosbell=(1-cos(2*pi*[0:R-1]/(R-1)))/2;
%     end
%     if (if_spec_cosbell==2)
%         cosbell=(1-cos(2*pi*[0:R-1]/R))/2;
%     end
%     cosbell_2d=cosbell'*cosbell;
    cosbell_2d=mlis_window_setup(R,if_spec_cosbell); %16Sep20
    for ipatch=1:patch_count_sum
         patch=patches_downsampled{iNR}(:,:,ipatch);
         %subtract pixel mean if requested
         if (if_spec_subpxlmean)
             patch=patch-mean(patches_downsampled{iNR},3);
         end
         %subtract ROI mean if requested
         if (if_spec_subroimean)
             patch=patch-mean(patch(:));
         end
         %cosine bell window if requested
         if (if_spec_cosbell>0)
             patch=patch.*cosbell_2d;
         end
         patch_fft=fft2(patch,R*fft2_padfactor,R*fft2_padfactor);
         patches_fft(:,:,ipatch)=patch_fft;
         ps_patches{iNR}(:,:,ipatch)=patch_fft.*conj(patch_fft)/(R.^4); %area is R^2, and power squares it again
    end
    %
    for iview=1:count_views
        for isubj=1:count_subjs
            ps_avg_vs{iNR}(:,:,iview,isubj)=mean(ps_patches{iNR}(:,:,find((patches_view{iNR}==iview)&(patches_subj{iNR}==isubj))),3);
        end
        %average, weight by patch
        ps_avg_v{iNR}(:,:,iview,1)=mean(ps_patches{iNR}(:,:,find(patches_view{iNR}==iview)),3);
        %average, weight by subject
        ps_avg_v{iNR}(:,:,iview,2)=mean(ps_avg_vs{iNR}(:,:,iview,s_havedata{iview}),4);
    end
    %
    %whiten
    %
    disp(sprintf('whitening %4.0f patches of size %4.0f x %4.0f using average type %s for whitening and binarizing blocks of %3.0f x %3.0f, stepping %1.0f',...
        patch_count_sum,R,R,whiten_avg_types_labels{whiten_avg_type_use},S,S,Sstep));
    patches_whitened{iNR}=zeros(R,R,patch_count_sum);
    %patches_binarized{iNR}=zeros(R,R,patch_count_sum); %10Aug
    %set aside space for binarized images; if starting phase is stepped then array length is one smaller
    nstep=(Sstep==0)+(Sstep==1)*S;
    patches_binarized{iNR}=cell(nstep,nstep);
    Rsample=cell(1,nstep);
    for ixs=1:nstep
        Rsample{ixs}=[ixs:(R+ixs-1)-S*(double(ixs>1))]; %how to sample an R x R array of blocks, starting on block ixs, in groups of S
    end
    for ixs=1:nstep
        for iys=1:nstep
            patches_binarized{iNR}{ixs,iys}=uint8(zeros(R/S-double(ixs>1),R/S-double(iys>1),patch_count_sum,regtype_count));
        end
    end
    for ipatch=1:patch_count_sum
        patch=patches_downsampled{iNR}(:,:,ipatch);
        isubj=patches_subj{iNR}(ipatch);
        iview=patches_view{iNR}(ipatch);
        switch whiten_avg_types_labels_orig{whiten_avg_type_use} %changed 23Oct20 so that 'power law' in whiten_avg_types_labels can be changed
            case 'per subj and view'
                ps_patch=ps_avg_vs{iNR}(:,:,iview,isubj);
            case 'per view (patch wtd)'
                ps_patch=ps_avg_v{iNR}(:,:,iview,1);
            case 'per view (subj wtd)'
                ps_patch=ps_avg_v{iNR}(:,:,iview,2);
            case 'global (patch wtd)'
                ps_patch=mean(ps_avg_v{iNR}(:,:,:,1),3);
            case 'global (subj wtd)'
                ps_patch=mean(ps_avg_v{iNR}(:,:,:,2),3);
            case 'per patch'
                ps_patch=ps_patches{iNR}(:,:,ipatch);
            case 'none'
                ps_patch=ones(fft_length,fft_length);
            case 'power law' %added 23Oct20
                %compute power law spectrum
                sf=min([0:fft_length-1],fft_length-[0:fft_length-1]);
                sfsq=repmat(sf.^2,fft_length,1);
                sfgrid=sqrt(sfsq+sfsq');
                sfgrid(1,1)=1; %avoid zero-divide
                ps_patch=sfgrid.^whiten_power_law;
        end
        %following Ann's code, as in prewhiten_mri_demo
        filt=1./sqrt(real(ps_patch));
        filt(1,1)=0; %ps_patch(1,1) is 0 so its recip is infinite
        patch_whitened_padded=real(ifft2(patches_fft(:,:,ipatch).*filt));
        patch_whitened=patch_whitened_padded(1:R,1:R);
        patches_whitened{iNR}(:,:,ipatch)=patch_whitened;
        %use the scale S (10Aug20)
        %patches_binarized{iNR}(:,:,ipatch)=double(patch_whitened>=median(patch_whitened(:)));
        %
        for ixs=1:nstep
            for iys=1:nstep
                patch_sampled=patch_whitened(Rsample{ixs},Rsample{iys}); %choose the subset, stepping the sampling phase
                patch_scaled=squeeze(mean(mean(reshape(patch_sampled,[S,size(patch_sampled,1)/S,S,size(patch_sampled,2)/S]),3),1));
                %here need to calculate patches_binarized for each kind of region
                for iregtype=1:regtype_count
%                     rect_mask_scaled=zeros(size(patch_scaled)); %zeros indicate pixels to be used
                      if (regtype_binarize_useall==0) &  (regtype_masks{iregtype}(1)>0)
                        %retrieve composite mask
                        patch_pair_ID=ceil(ipatch/patch_mult(iNR));
                        [rect_mask_composite,rect_mask_sampled,rect_mask_scaled]=...
                            nibs_get_mask(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,nibs_patches_img_size(patch_pair_ID,:),nibs_NRpatch_max,...
                            regtype_masks{iregtype},setfields([],{'xs','ys','S'},{Rsample{ixs},Rsample{iys},S}));
                        median_cut=median(patch_scaled(rect_mask_scaled(:)==0)); %find cutpoint only on patch pixels that are selected by mask
                    else
                        median_cut=median(patch_scaled(:));
                    end
 %                   patches_binarized{iNR}{ixs,iys}(:,:,ipatch,iregtype)=uint8(patch_scaled>median(patch_scaled(:)));
                    patches_binarized{iNR}{ixs,iys}(:,:,ipatch,iregtype)=uint8(patch_scaled>median_cut);
                 end
            end
        end
    end
    if (if_clear)
        clear patches_downsampled
    end
    %show one patch
    if (if_plot_pipeline)
        icol=4;
        for irow=1:nrows_patch
            subplot(nrows_patch,ncols_patch,icol+(irow-1)*ncols_patch);
            imagesc(patches_whitened{iNR}(1:R,1:R,patch_index_show(irow)));
            colormap gray;
            axis square;
            axis tight;
            set(gca,'XTick',[1 R]);
            set(gca,'YTick',[1 R]);
            title(sprintf('whitened, %s',whiten_avg_types_labels{whiten_avg_type_use}));
            hold on;
        end %irow
    end %if_plot_pipeline
    if (if_clear)
        clear patches_whitened
    end
    %show binarized
    if (if_plot_pipeline)
        icol=5;
        for irow=1:nrows_patch
            subplot(nrows_patch,ncols_patch,icol+(irow-1)*ncols_patch);
            imagesc(patches_binarized{iNR}{1,1}(:,:,patch_index_show(irow)),[0 1]); %10Aug29=0: {1,1} added to accommodate stepping
            colormap gray;
            axis square;
            axis tight;
            set(gca,'XTick',[1 R]);
            set(gca,'YTick',[1 R]);
            title(sprintf('binarized'));
            hold on;
        end %irow
    end %if_plot_pipeline
    ps_max=max(ps_patches{iNR}(:));
    if (if_plot_spectra)
        %
        %a figure with power spectra
        %
        hf_ps=figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'spectra: ',tstring));
        set(gcf,'Position',[50 50 1200 750]);
        nrows_ps=2;
        ncols_ps=3;
        icol=1; %spectra per patch, horiz or vert
        for irow=1:nrows_ps
            hline=[];
            htext=[];
            subplot(nrows_ps,ncols_ps,icol+(irow-1)*ncols_ps)
            plot_spec_thinfac_use=max(plot_spec_thinfac,floor(patch_count_sum/plot_spec_maxnum));
            for ipatch=1:patch_count_sum
                ps=fftshift(ps_patches{iNR}(:,:,ipatch));
                color=view_colors{view_nos(patches_view{iNR}(ipatch))};
                %thin out the plot but always plot the first patch of each view
                if_plot_indiv=(ipatch==min(find(patches_view{iNR}==patches_view{iNR}(ipatch))));
                if (if_plot_indiv==0)
                    if_plot_indiv=(mod(ipatch-1,plot_spec_thinfac_use)==0);
                end
                if (if_plot_indiv)
                    switch irow
                        case 1
                            ps_vals=ps(fft_half+1+[0:fft_half-1],fft_half+1)'; %starts at zero spatial frequency
                        case 2
                            ps_vals=ps(fft_half+1,fft_half+1+[0:fft_half-1]);
                    end
                    hplot=loglog(sfreqs_per_pixel(2:fft_half),ps_vals(2:end),color); %plot nonzero spatial frequency
                    hold on;
                    if ipatch==min(find(patches_view{iNR}==patches_view{iNR}(ipatch)))
                        hline=[hline;hplot];
                        htext=strvcat(htext,view_names{patches_view{iNR}(ipatch)});
                    end
                end
            end
            hleg=legend(hline,htext,'Location','NorthEast');
            set(gca,'XLim',[1/NR,0.5]);
            set(hleg,'FontSize',7);
            set(gca,'YLim',ps_max*[ps_minfac,1]);
            title(sprintf('patch (%s) thinfac %1.0f',hv_labels{irow},plot_spec_thinfac_use));
            xlabel('cy/pixel');
            ylabel('pwr');
        end %irow
        icol=2; %spectra per subj, horiz or vert
        for irow=1:nrows_ps
            hline=[];
            htext=[];
            subplot(nrows_ps,ncols_ps,icol+(irow-1)*ncols_ps)
            for isubj=1:count_subjs
                for iview=1:count_views
                    ps=fftshift(ps_avg_vs{iNR}(:,:,iview,isubj));
                    color=view_colors{view_nos(iview)};
                    switch irow
                        case 1
                            ps_vals=ps(fft_half+1+[0:fft_half-1],fft_half+1)'; %starts at zero spatial frequency
                        case 2
                            ps_vals=ps(fft_half+1,fft_half+1+[0:fft_half-1]);
                    end
                    hplot=loglog(sfreqs_per_pixel(2:fft_half),ps_vals(2:end),color); %plot nonzero spatial frequency
                    hold on;
                    if (isubj==1)
                        hline=[hline;hplot];
                        htext=strvcat(htext,view_names{view_nos(iview)});
                    end
                end
            end
            hleg=legend(hline,htext,'Location','NorthEast');
            set(gca,'XLim',[1/NR_max,0.5]);
            set(hleg,'FontSize',7);
            set(gca,'YLim',ps_max*[ps_minfac,1]);
            title(sprintf('subj (%s)',hv_labels{irow}));
            xlabel('cy/pixel');
            ylabel('pwr');
        end %irow
        icol=3; %averaged spectra
        for irow=1:nrows_ps
            hline=[];
            htext=[];
            subplot(nrows_ps,ncols_ps,icol+(irow-1)*ncols_ps)
            for ps_avg_type=1:ps_avg_types
                for iview=1:count_views
                    ps=fftshift(ps_avg_v{iNR}(:,:,iview,ps_avg_type));
                    color=view_colors{view_nos(iview)};
                    switch irow
                        case 1
                            ps_vals=ps(fft_half+1+[0:fft_half-1],fft_half+1)'; %starts at zero spatial frequency
                        case 2
                            ps_vals=ps(fft_half+1,fft_half+1+[0:fft_half-1]);
                    end
                    hplot=loglog(sfreqs_per_pixel(2:fft_half),ps_vals(2:end),cat(2,color,ps_avg_types_symbs{ps_avg_type}));
                    hold on;
                    hline=[hline;hplot];
                    htext=strvcat(htext,cat(2,view_names{view_nos(iview)},' ',ps_avg_types_labels{ps_avg_type}));
                end
            end
            h2=loglog(1/NR*[1 10],ps_max*[1 10^-2],'k:');
            h3=loglog(1/NR*[1 10],ps_max*[1 10^-3],'k-');
            hline=[hline;h2;h3];
            htext=strvcat(htext,'-2','-3');
            hleg=legend(hline,htext,'Location','NorthEast');
            set(gca,'XLim',[1/NR,0.5]);
            set(hleg,'FontSize',7);
            set(gca,'YLim',ps_max*[ps_minfac,1]);
            title(sprintf('avgs (%s)',hv_labels{irow}));
            xlabel('cy/pixel');
            ylabel('pwr');
        end %irow
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end %if_plot_spectra
    if (if_clear)
        clear ps_patches
    end
    %
    %a figure to look at anisotropies in power spectra
    %
    if (if_plot_details) & ismember(N,N_list_psa)
       for ps_avg_type=1:ps_avg_types
            hf_psa=figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'spectra details: ',tstring));
            set(gcf,'Position',[50 50 1200 750]);
            nrows_psa=2;
            ncols_psa=count_views;
            for iview=1:count_views
                icol=iview;
                irow=1; %heatmap
                subplot(nrows_psa,ncols_psa,icol+(irow-1)*ncols_psa)
                imagesc(max(log10(ps_max)-ps_log10range,log10(fftshift(ps_avg_v{iNR}(:,:,iview,ps_avg_type)))),log10(ps_max)+[-ps_log10range 0]);
                title(sprintf('%s (avg %s)',view_names{view_nos(iview)},ps_avg_types_labels{ps_avg_type}));
                set(gca,'XTick',0.5+[0 fft_half fft_length]);
                set(gca,'XTickLabel',[-1/(2*N) 0 1/(2*N)]); 
                set(gca,'XLim',0.5+[0 fft_length]);
                set(gca,'YTick',0.5+[0 fft_half fft_length]);
                set(gca,'YTickLabel',[-1/(2*N) 0 1/(2*N)]); 
                set(gca,'YLim',0.5+[0 fft_length]);
                set(gca,'YDir','normal') %so that positive Y is plotted at the top, in contrast to Matlab'e standard imagesc convention
                axis square;
                %
                irow=2; %compare H with V, and oblique axes
                subplot(nrows_psa,ncols_psa,icol+(irow-1)*ncols_psa)
                ps=fftshift(ps_avg_v{iNR}(:,:,iview,ps_avg_type));
                ps_v=ps(fft_half+1+[0:fft_half-1],fft_half+1)'; %starts at zero spatial frequency
                ps_h=ps(fft_half+1,fft_half+1+[0:fft_half-1]);
                loglog(sfreqs_per_pixel(2:fft_half),ps_v(2:end)./ps_h(2:end),view_colors{view_nos(iview)}); %plot starting at nonzero spatial frequency
                hold on;
                ps_ssdiag=diag(ps);
                ps_osdiag=diag(flipud(ps));
                ps_ss=ps_ssdiag(1:fft_half);
                ps_os=ps_osdiag(1:fft_half);
                loglog(sfreqs_per_pixel(2:fft_half)*sqrt(2),ps_ss(2:end)./ps_os(2:end),cat(2,view_colors{view_nos(iview)},':'));
                loglog([sfreqs_per_pixel(2) sfreqs_per_pixel(fft_half)*sqrt(2)],[1 1],'k');
                hold on;
                legend('vert/horiz','same-sign/opp-sign','unity');
                xlabel('cy/pixel');
                ylabel(cat(2,'power ratio',' (avg ',ps_avg_types_labels{ps_avg_type},')'));
                set(gca,'YLim',10.^[-1 1]);
            end %irow
            axes('Position',[0.02,0.02,0.01,0.01]); %for text
            text(0,0,tstring,'Interpreter','none');
            axis off;
       end %avg type
    end %if_plot_details & ismember(N,N_list_psa)
    %
    %calculate block counts
    %July 2022: this is done separately for each region type (regtype_count)
    %
    patches_blockcounts{iNR}=zeros(btc_nconfigs,patch_count_sum,regtype_count);
    patches_btc{iNR}=zeros(btc_n,patch_count_sum,regtype_count);
    for iregtype=1:regtype_count
        disp(sprintf('calculating block counts for %4.0f patches of size %4.0f x %4.0f (%2.0f x %2.0f downsampled), margin %1.0f, derived from %1.0f big patches of size %4.0f x %4.0f pixels',...
            patch_count_sum,R,R,N,N,margin_width,patch_count_sum_big,NR_max,NR_max));
        disp(datasource_string);       
        if if_regsel
            disp(sprintf('region type %2.0f: %s',iregtype,regtype_descs{iregtype}));
        end
        for ipatch=1:patch_count_sum
            %07Sep22: add masking if if_regsel=1
            %10Aug20: add stepping loop, accumulating block counts
            blockcounts=zeros(1,btc_nconfigs);
            for ixs=1:nstep
                for iys=1:nstep
                    patch=double(patches_binarized{iNR}{ixs,iys}(:,:,ipatch,iregtype)); %needed so that NaN's can be set
                    if (if_regsel)
                        %rect_mask_scaled is 0 to select, 1 to deselect, takes S into account
                            patch_pair_ID=ceil(ipatch/patch_mult(iNR));
                            [rect_mask_composite,rect_mask_sampled,rect_mask_scaled]=...
                            nibs_get_mask(nibs_bpdata,nibs_aspect_ptr,nibs_smoothness_ptr,patch_pair_ID,nibs_patches_img_size(patch_pair_ID,:),nibs_NRpatch_max,...
                            regtype_masks{iregtype},setfields([],{'xs','ys','S'},{Rsample{ixs},Rsample{iys},S}));
                        patch(rect_mask_scaled==1)=NaN;                       
                    end                   
                    patch_marg=patch((1+margin_width):(end-margin_width),(1+margin_width):(end-margin_width));
                    blockcounts=blockcounts+glider_mapubi(patch_marg,btc_checkdef,btc_ng,setfield([],'mapubi_bc',0)); %NON-periodic boundary conditions
                end
            end
            patches_blockcounts{iNR}(:,ipatch,iregtype)=blockcounts;
            if sum(blockcounts)>0 %the selected region could be empty
                p2x2=reshape(blockcounts,btc_reshape)/sum(blockcounts(:));
                patches_btc{iNR}(:,ipatch,iregtype)=btc_corrs2vec(getcorrs_p2x2(p2x2,0,1),btc_dict); %compute only btc stats (07Feb20)
            end
        end
    end %iregtype
    if (if_clear)
        clear patches_binarized
    end
end % iNR
if (if_clear)
    clear patches_orig
    clear patches_fft patch_fft
    clear ps_avg_vs  %ps_avg_v %ps_avg_v not cleared since we will use this in ffdm_spec_stats
    clear roi_img_examp
end
%
% calculate binary image statistics
%
btc_avg_vs=zeros(length(N_list),btc_n,count_views,count_subjs,regtype_count); %average, each view and subject
btc_avg_v=zeros(length(N_list),btc_n,count_views,btc_avg_types,regtype_count);
btc_std_v=zeros(length(N_list),btc_n,count_views,btc_avg_types,regtype_count);
btc_avg_allv=zeros(length(N_list),btc_n,btc_avg_types,regtype_count);
btc_std_allv=zeros(length(N_list),btc_n,btc_avg_types,regtype_count);
btc_nus_v=zeros(length(N_list),count_views,btc_avg_types,regtype_count);
tstring_NR=sprintf('%s, NR_max=%1.0f pad=%1.0f whiten=%s margin=%1.0f gen usefrac=%1.0f',...
    datasource_string,NR_max,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,patches_usefrac);
if (if_norand==1)
    tstring_NR=cat(2,tstring_NR,' NO randomization');
end
if (if_norand==-1)
    tstring_NR=cat(2,tstring_NR,' Frozen randomization');
end
%for jackknifing based on subjects
btc_avg_v_drop=zeros(length(N_list),btc_n,count_views,btc_avg_types,count_subjs,regtype_count);
btc_std_v_drop=zeros(length(N_list),btc_n,count_views,btc_avg_types,count_subjs,regtype_count);
btc_avg_v_jsem=zeros(length(N_list),btc_n,count_views,btc_avg_types,regtype_count);
btc_std_v_jsem=zeros(length(N_list),btc_n,count_views,btc_avg_types,regtype_count);
btc_avg_allv_drop=zeros(length(N_list),btc_n,btc_avg_types,count_subjs,regtype_count);
btc_std_allv_drop=zeros(length(N_list),btc_n,btc_avg_types,count_subjs,regtype_count);
btc_avg_allv_jsem=zeros(length(N_list),btc_n,btc_avg_types,regtype_count);
btc_std_allv_jsem=zeros(length(N_list),btc_n,btc_avg_types,regtype_count);
%
tmult=tinv((1+cl_prob)/2,count_subjs); %two-tailed t-test crit value 
counts_per_patch=zeros(1,length(N_list)); %total counts in a full patch
std_counts=zeros(length(N_list),count_views,btc_avg_types,regtype_count); %expected standard dev based on counts
%need to do this for each region type (masks already applied) and use 
for iNR=1:length(N_list)
    N=N_list(iNR);
    NR=NR_list(iNR);
    R=NR/N;
    S=S_list(iNR);
    patch_count_sum=patch_count_sum_big*patch_mult(iNR); %total number of patches, all views, all subjects
    counts_per_patch(iNR)=(R/S-2*margin_width-1).^2; %number of block counts per patch %modified 10Aug20 to take into account blocking
    std_counts_perpatch=sqrt(1./counts_per_patch(iNR));  %standard dev of a btc statistic calculated from a patch without masking, was std_counts(iNR,iview,1)
    have_patches=zeros(count_views,1); %how many patches with each view
    have_btc_patchlist=cell(count_views,regtype_count); %list of patches with btc data for each view (only includes patches with sufficient pixels in region)
    have_btc_subjlist=cell(count_views,regtype_count); %list of subjects with btc data for each view and region
    for iview=1:count_views
        have_patches(iview)=sum(patches_view{iNR}==iview);
        std_counts_persubj=sqrt(1./counts_per_patch(iNR)*mean(1./(patch_count(s_havedata{iview},iview)*patch_mult(iNR)))); %patch_mult is patch_mult_side.^2; was std_counts(iNR,iview,2)
        %patch weighted -- when region type is present, this takes into account the number of pixels in the selected regions
        for iregtype=1:regtype_count
            have_btc_patchlist{iview,iregtype}=intersect(find(patches_view{iNR}==iview),find(sum(patches_blockcounts{iNR}(:,:,iregtype),1)'>0)); %which patches have at least some counts of btc stats
            have_btc_subjlist{iview,iregtype}=unique(patches_subj{iNR}(have_btc_patchlist{iview,iregtype}));% which subjects have data for this view and region; same as s_havedata{iview} for if_regsel=0
            disp(sprintf(' iNR=%3.0f, [N,R,S]=[%3.0f %5.0f %3.0f], view %2.0f region type %2.0f (%s)',iNR,N,R,S,iview,iregtype,regtype_descs{iregtype}));
            disp(sprintf('    number of contributing  patches: %5.0f of %5.0f (this view) and %5.0f total',length(have_btc_patchlist{iview,iregtype}),have_patches(iview),patch_count_sum));
            disp(sprintf('    number of contributing subjects: %5.0f of %5.0f (this view) and %5.0f total',length(have_btc_subjlist{iview,iregtype}),length(s_havedata{iview}),count_subjs));
            %
            btc_nus_v(iNR,iview,1,iregtype)=length(have_btc_patchlist{iview,iregtype}); % was btc_nus_v(iNR,iview,1)=sum(patches_view{iNR}==iview);
            btc_avg_v(iNR,:,iview,1,iregtype)=mean(patches_btc{iNR}(:,have_btc_patchlist{iview,iregtype},iregtype),2); % was mean(patches_btc{iNR}(:,find(patches_view{iNR}==iview)),2);
            btc_std_v(iNR,:,iview,1,iregtype)=std(patches_btc{iNR}(:,have_btc_patchlist{iview,iregtype},iregtype),0,2); % was std(patches_btc{iNR}(:,find(patches_view{iNR}==iview)),0,2);
            %
            if length(have_btc_patchlist{iview,iregtype})>0
                std_counts(iNR,iview,1,iregtype)=sqrt(mean(1./sum(patches_blockcounts{iNR}(:,have_btc_patchlist{iview,iregtype},iregtype),1))); % was std_counts(iNR,iview,1,1)=sqrt(1./counts_per_patch(iNR));
                %for avg_types=2, sum up all of the pixel counts for each subject
                subj_list=unique(patches_subj{iNR});
                counts_svr=zeros(1,length(subj_list)); %counts for this subject, view, and region
                for subj_ptr=1:length(subj_list)
                    patch_list=intersect(have_btc_patchlist{iview,iregtype},find(patches_subj{iNR}==subj_list(subj_ptr)));
                    counts_svr(subj_ptr)=sum(sum(patches_blockcounts{iNR}(:,patch_list,iregtype))); %all the blockcounts in all patches for this view, subject, and region
                end
                counts_svr_nz=find(counts_svr>0);
                if (length(counts_svr_nz)>0)
                    std_counts(iNR,iview,2,iregtype)=sqrt(mean(1./counts_svr(counts_svr_nz)));
                end
            end
            disp(sprintf('    estimated standard deviation based on counting errors for per-patch stats: %7.5f (without region masking: %7.5f)',...
                std_counts(iNR,iview,1,iregtype),std_counts_perpatch));
            disp(sprintf('    estimated standard deviation based on counting errors for per-subj  stats: %7.5f (without region masking: %7.5f)',...
                std_counts(iNR,iview,2,iregtype),std_counts_persubj));
        end %regtype      
        for isubj=1:count_subjs
            btc_avg_vs(iNR,:,iview,isubj,:)=mean(patches_btc{iNR}(:,find((patches_view{iNR}==iview)&(patches_subj{iNR}==isubj)),:),2); %regions calculated in parallel
        end       
        for iregtype=1:regtype_count
            %btc_nus_v(iNR,iview,2,iregtype)=sum(s_havedata{iview}>0); %changed from count_subjs, 04Dec19
            btc_nus_v(iNR,iview,2,iregtype)=length(have_btc_subjlist{iview,iregtype}); % was btc_nus_v(iNR,iview,2,iregtype)=sum(s_havedata{iview}>0); %changed from count_subjs, 04Dec19
            if length(have_btc_subjlist{iview,iregtype})>=1
                btc_avg_v(iNR,:,iview,2,iregtype)=mean(btc_avg_vs(iNR,:,iview,have_btc_subjlist{iview,iregtype},iregtype),4); %was sum(s_havedata{iview}>0); %changed from count_subjs, 04Dec19
            end
            if length(have_btc_subjlist{iview,iregtype})>=2
                btc_std_v(iNR,:,iview,2,iregtype)=std(btc_avg_vs(iNR,:,iview,have_btc_subjlist{iview,iregtype},iregtype),0,4); %std(btc_avg_vs(iNR,:,iview,s_havedata{iview}),0,4) normalize by count-1
            end
        end
        %drop-one quantities
        for iregtype=1:regtype_count %section modified to take into account regtype
        %note that this is computed for all subjects so that all of the "drops" are available for a cross-view average
        %jackknifing within views (called by jack) only considers the dropped subjects that have data
            for isubj=1:count_subjs
                %weight patches equally
                patchlist_drop=setdiff(have_btc_patchlist{iview,iregtype},find(patches_subj{iNR}==isubj)); %only the patches that are not from this subject
                btc_avg_v_drop(iNR,:,iview,1,isubj,iregtype)=mean(patches_btc{iNR}(:,patchlist_drop,iregtype),2);
                btc_std_v_drop(iNR,:,iview,1,isubj,iregtype)=std(patches_btc{iNR}(:,patchlist_drop,iregtype),0,2);
                %weight subjects equally
                subjlist_drop=setdiff(have_btc_subjlist{iview,iregtype},isubj);
                btc_avg_v_drop(iNR,:,iview,2,isubj,iregtype)=mean(btc_avg_vs(iNR,:,iview,subjlist_drop,iregtype),4);
                btc_std_v_drop(iNR,:,iview,2,isubj,iregtype)=std(btc_avg_vs(iNR,:,iview,subjlist_drop,iregtype),0,4); %normalize by count-1          
            end %isubj
        end %iregtype
        %
    end %iview
    %across views
    btc_avg_allv(iNR,:,:,:)=reshape(mean(btc_avg_v(iNR,:,:,:),3),[1 btc_n btc_avg_types regtype_count]); %average across views
    btc_std_allv(iNR,:,1,:)=reshape(std(patches_btc{iNR},0,2),[1 btc_n 1 regtype_count]); %patches_btc{iNR} is [btc_n, patch_count_sum regtype_count], includes all views
    %
    %across views, per-subject, taking into account that not every subject has every view with every region
    have_btc_patchlist_allv=cell(1,regtype_count); %list of patches that have data, this view, for any view 
    have_btc_subjlist_allv=cell(1,regtype_count);
    for iregtype=1:regtype_count %section modified to take into account regtype
        btc_avg_vs_allv=zeros(1,btc_n,1,0);
        for iview=1:count_views
            have_btc_patchlist_allv{iregtype}=union(have_btc_patchlist_allv{iregtype},have_btc_patchlist{iview,iregtype});
        end
        have_btc_subjlist_allv{iregtype}=unique(patches_subj{iNR}(have_btc_patchlist_allv{iregtype}));
        %
        %btc_avg_vs=zeros(length(N_list),btc_n,count_views,count_subjs,regtype_count);
        for iview=1:count_views
            subj_list=have_btc_subjlist{iview,iregtype};
            btc_avg_vs_allv=cat(4,btc_avg_vs_allv,btc_avg_vs(iNR,:,iview,subj_list,iregtype)); %concatenate data from all subjects that have data with this view
        end
        btc_std_allv(iNR,:,2,iregtype)=std(btc_avg_vs_allv,0,4);
        %%%TO HERE
        %drop-one quantities
        btc_avg_allv_drop(iNR,:,:,:)=mean(btc_avg_v_drop(iNR,:,:,:,:),3);
        for isubj=1:count_subjs
            subj_list=setdiff([1:count_subjs],isubj);
            btc_std_allv_drop(iNR,:,1,isubj)=std(patches_btc{iNR}(:,find(patches_subj{iNR}~=isubj)),0,2); %bug fixed 22 Nov 19:  LHS did not have isubj
            %take into account that not every subject has every view
            %previous: btc_std_allv_drop(iNR,:,2,isubj)=std(reshape(btc_avg_vs(iNR,:,:,subj_list),[1 btc_n,count_views*(count_subjs-1)]),0,3);
            btc_avg_vs_allv_drop=zeros(1,btc_n,1,0);
            for iview=1:count_views
                subj_list=setdiff(s_havedata{iview},isubj);
                btc_avg_vs_allv_drop=cat(4,btc_avg_vs_allv_drop,btc_avg_vs(iNR,:,iview,subj_list));
            end
            btc_std_allv_drop(iNR,:,2,isubj)=std(btc_avg_vs_allv_drop,0,4);
        end
    end %iregtype
    %
end
%jackknife standard error of measurement
%modified so that we don't jackknife the subjects that don't have any data for a view
for iview=1:count_views
    if length(s_havedata{iview})>1
        [jbias,jdebiased,jvar,btc_avg_v_jsem(:,:,iview,:)]=jack(btc_avg_v(:,:,iview,:),btc_avg_v_drop(:,:,iview,:,s_havedata{iview}));
        [jbias,jdebiased,jvar,btc_std_v_jsem(:,:,iview,:)]=jack(btc_std_v(:,:,iview,:),btc_std_v_drop(:,:,iview,:,s_havedata{iview}));
    end
end
if length(s_havedata_allv)>1
    [jbias,jdebiased,jvar,btc_avg_allv_jsem]=jack(btc_avg_allv,btc_avg_allv_drop(:,:,:,s_havedata_allv));
    [jbias,jdebiased,jvar,btc_std_allv_jsem]=jack(btc_std_allv,btc_std_allv_drop(:,:,:,s_havedata_allv));
end
%
btc_header='  N  R   S';
btc_fmt='%3.0f %3.0f %2.0f%s';
btc_stat_fmt='%8s  ';
for icoord=1:btc_n
    btc_header=cat(2,btc_header,'       ',btc_order_JV(icoord));
    btc_fmt=cat(2,btc_fmt,'%7.4f ');
    btc_stat_fmt=cat(2,btc_stat_fmt,' %7.4f');
end
disp(' ');
disp(tstring_NR);
%NEED TO DO THIS SEPARATELY FOR EACH REGION TYPE
for btc_avg_type=1:btc_avg_types %1->patch weighted, 2->subject weighted
    disp(' ');
    disp(tstring_NR);
    disp(sprintf('Image statistics calculated with weighting %s; sems and CLs (%7.5f) by subject -based jackknife',...
        btc_avg_types_labels{btc_avg_type},cl_prob));
    disp(sprintf('texture coordinate averages'));
    for iview=1:count_views
        disp(sprintf(' view %5s',view_names{view_nos(iview)}));
        disp(btc_header);
        for iNR=1:length(N_list);
            N=N_list(iNR);
            NR=NR_list(iNR);
            R=NR/N;
            S=S_list(iNR);
            Splus=Sstep_char{Sstep_list(iNR)+1};
            disp(sprintf(' from %6.0f examples:',btc_nus_v(iNR,iview,btc_avg_type)));
            disp(sprintf(btc_fmt,N,R,S,Splus,btc_avg_v(iNR,:,iview,btc_avg_type)));
            disp(sprintf(btc_stat_fmt,'jk_sem',btc_avg_v_jsem(iNR,:,iview,btc_avg_type)));
            if (cl_prob>0)
                disp(sprintf(btc_stat_fmt,'cl_lo',btc_avg_v(iNR,:,iview,btc_avg_type)-tmult*btc_avg_v_jsem(iNR,:,iview,btc_avg_type)));
                disp(sprintf(btc_stat_fmt,'cl_hi',btc_avg_v(iNR,:,iview,btc_avg_type)+tmult*btc_avg_v_jsem(iNR,:,iview,btc_avg_type)));
            end
        end %iNR
    end %iview
    disp(sprintf(' average across all views'));
    disp(btc_header);
    for iNR=1:length(N_list);
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        disp(sprintf(' from %6.0f examples',sum(btc_nus_v(iNR,:,btc_avg_type))));
        disp(sprintf(btc_fmt,N,R,S,Splus,btc_avg_allv(iNR,:,btc_avg_type)));
        disp(sprintf(btc_stat_fmt,'jk_sem',btc_avg_allv_jsem(iNR,:,btc_avg_type)));
            if (cl_prob>0)
                disp(sprintf(btc_stat_fmt,'cl_lo',btc_avg_allv(iNR,:,btc_avg_type)-tmult*btc_avg_allv_jsem(iNR,:,btc_avg_type)));
                disp(sprintf(btc_stat_fmt,'cl_hi',btc_avg_allv(iNR,:,btc_avg_type)+tmult*btc_avg_allv_jsem(iNR,:,btc_avg_type)));
            end
    end %iNR
    disp(sprintf('texture coordinate standard devs'));
    for iview=1:count_views
        disp(sprintf(' view %5s',view_names{view_nos(iview)}));
        disp(btc_header);
        for iNR=1:length(N_list);
            N=N_list(iNR);
            NR=NR_list(iNR);
            R=NR/N;
            S=S_list(iNR);
            Splus=Sstep_char{Sstep_list(iNR)+1};
            disp(sprintf(' from %6.0f examples',btc_nus_v(iNR,:,btc_avg_type)));
            disp(sprintf(btc_fmt,N,R,S,Splus,btc_std_v(iNR,:,iview,btc_avg_type)));
            disp(sprintf(btc_stat_fmt,'jk_sem',btc_std_v_jsem(iNR,:,iview,btc_avg_type)));
            if (cl_prob>0)
                disp(sprintf(btc_stat_fmt,'cl_lo',btc_std_v(iNR,:,iview,btc_avg_type)-tmult*btc_std_v_jsem(iNR,:,iview,btc_avg_type)));
                disp(sprintf(btc_stat_fmt,'cl_hi',btc_std_v(iNR,:,iview,btc_avg_type)+tmult*btc_std_v_jsem(iNR,:,iview,btc_avg_type)));
            end
            disp(sprintf('expected contribution to standard dev based on counting errors within patches: %7.4f (patches have %5.0f counts)',...
                std_counts(iNR,iview,btc_avg_type),counts_per_patch(iNR)));
        end %iNR
    end %iview
    disp(sprintf(' standard dev across all views'));
    disp(btc_header);
    for iNR=1:length(N_list);
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        disp(sprintf(' from %6.0f examples',sum(btc_nus_v(iNR,:,btc_avg_type))));
        disp(sprintf(btc_fmt,N,R,S,Splus,btc_std_allv(iNR,:,btc_avg_type)));
        disp(sprintf(btc_stat_fmt,'jk_sem',btc_std_allv_jsem(iNR,:,btc_avg_type)));
        if (cl_prob>0)
            disp(sprintf(btc_stat_fmt,'cl_lo',btc_std_allv(iNR,:,btc_avg_type)-tmult*btc_std_allv_jsem(iNR,:,btc_avg_type)));
            disp(sprintf(btc_stat_fmt,'cl_hi',btc_std_allv(iNR,:,btc_avg_type)+tmult*btc_std_allv_jsem(iNR,:,btc_avg_type)));
        end
    end %iNR
end %btc_avg_type
