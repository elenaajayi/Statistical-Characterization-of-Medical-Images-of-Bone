%ffdm_btc_plot_nistats: basic plots for binary texture coordinate analysis of natural images
%  computed by Hermundstad et al. pipeline, and saved in NIstats.mat
%
%  Plots will NOT match those of Figure 3 in Hermundstad et al. since here we plot raw statistics, and in that plot,
%   each set of image stats is multiplied by a scale factor to best match image stats
%
% NIstats.mat has the local image statistics from  each analyzed patch
%  This is data_source=6
%
%  In comparison with images processed through ffdm_btc_calc_gen:
%  * whitening has already been done based on global spectra of patches
%  * whitening used no cosine bell, no padding, no subtraction of pixel or patch average
%  * binarization was at median (as standard), but 2x2 neighborhoods are
%     non-overlapping rather than staggered
%  * there is no grouping into "subjects"
%  * confidence limits are not calculated since it is unclear how the individual patches are grouped into images
%  *
%  Source designation: here (natural images) it is 6:
%  * data_source= 1,2,3: full-field digital mammography (ffdm), for pilot data, original patches, and processed patches
%  * data_source=4: carotid ultrasound (usic, from Amanda Simon
%  * data_source=5: MRI (mrix), from Sam Xu
%
%  See also:  BTC_DEFINE, NISTATS_PLOT_DEMO, FFDM_NRPLOT, FFDM_BTC_PLOT_GEN.
%
if ~exist('btc_dict') btc_dict=btc_define([]); end
btc_n=length(btc_dict.codel); %10
btc_eLife_order='gcbdevwtua';
if ~exist('dropbox_path') dropbox_path='../../../Dropbox'; end %also  ../../../../Dropbox
if ~exist('dropbox_dir') dropbox_dir='/Textures/RawNIstats/ImagePatchStats'; end %also /Textures/RawNIstats/ImagePatchStats/Data;
if ~exist('filename') filename='NIstats.mat'; end
dropbox_path=getinp('dropbox path','s',[],dropbox_path);
dropbox_dir=getinp('dropbox dir','s',[],dropbox_dir);
filename=getinp('file containing NI data','s',[],filename);
fullname=cat(2,dropbox_path,filesep,dropbox_dir,filesep,filename);
fullname=strrep(fullname,'/',filesep);
fullname=strrep(fullname,cat(2,filesep,filesep),filesep);
fullname_show=strrep(fullname,filesep,'/');
NIstats=getfield(load(fullname),'NIstats')

if ~exist('data_source') | ~exist('N_list') | ~exist('if_allv_reorg') ...
        | ~exist('fft2_padfactor') | ~exist('whiten_avg_type_use') | ~exist('if_norand') ...
        | ~exist('preprocess_downsample') | ~exist('btc_order_HX') | ~exist('view_names') | ~exist('btc_avg_types')
    %should not be set, as they are usually set by ffdm_btc_calc_gen
    %data source
    data_source=6;
    datasource_string='NIstats';
    % whitening details, from Hermundstad et al. paper 
    % and ../Textures/Texture_Analysis_II/localFourierWhiten.m
    % and ../Textures/Texture_Analysis_II/processBlock.m
    fft2_padfactor=1; %no padding in Hermundstad et al. code
    whiten_avg_types_labels={'per subj and view','per view (patch wtd)','per view (subj wtd)','global (patch wtd)','global (subj wtd)','per patch'};
    whiten_avg_type_use=4; %global specified in Hermundstad et al.
    margin_width=0; %no margin .
    if_spec_subpxlmean=0; %from code
    if_spec_subroimean=0; %from code
    if_spec_cosbell=0; %from code
    if_norand=1; %no randomization mentioned
    preprocess_downsample=1;
    %
    N_list=NIstats.N;
    R_list=NIstats.R;
    NR_list=N_list.*R_list;
    %
    btc_order_HX=btc_eLife_order(2:end);
    %
    view_names={'std'};
    count_views=length(view_names);
    count_subjs=1; %"subjects" is irrelevant; don't know how patches are grouped into images
    if_allv_reorg=0;
    %
    btc_avg_types_labels={'per patch','per subj'}; %ways to average binary texture coords across patches
    btc_avg_types=length(btc_avg_types_labels);
    btc_avg_types_symbs={'-',':'};
else
    warning('Some quantities have already been unexpectedly set')
end
%
%calculate statistics from patch values in NIstata
%
%needed for ffdm_btc_plot_gen:
%
btc_avg_v=zeros(length(N_list),btc_n,count_views,btc_avg_types);
btc_std_v=zeros(length(N_list),btc_n,count_views,btc_avg_types);
btc_avg_v_jsem=NaN(length(N_list),btc_n,count_views,btc_avg_types);
btc_std_v_jsem=NaN(length(N_list),btc_n,count_views,btc_avg_types);
%
btc_avg_allv=zeros(length(N_list),btc_n,btc_avg_types);
btc_std_allv=zeros(length(N_list),btc_n,btc_avg_types);
btc_avg_allv_jsem=NaN(length(N_list),btc_n,btc_avg_types);
btc_std_allv_jsem=NaN(length(N_list),btc_n,btc_avg_types);
%
patches_view=cell(1,length(N_list));
patches_btc=cell(1,length(N_list));
npatches=zeros(1,length(N_list));
btc_eLife_posit=zeros(1,btc_n);
for ilet=1:btc_n
    btc_eLife_posit(ilet)=find(btc_eLife_order==btc_dict.codel(ilet));
end
%extract data from NIstats and reformat into the quantitites
%(btc_[avg|std]*) expected by ffdm_btc_plot_gen
for iNR=1:length(N_list)
    if (NIstats.set(iNR).N ~= N_list(iNR))
        warning('mismatch of N for iNR=%3.0f',iNR);
    end
    if (NIstats.set(iNR).R ~= R_list(iNR))
        warning('mismatch of R for iNR=%3.0f',iNR);
    end
    NIdata=NIstats.set(iNR).stats;
    %dim 1 of NIdata are the individual patches
    %dim 2 of NIdata has the statistics in btc_eLife_order, i.e., [gcbdevwtua]
    npatches(iNR)=size(NIdata,1);
    patches_view{iNR}=ones(1,npatches(iNR));
    patches_btc{iNR}=NIdata(:,btc_eLife_posit)';
    for iview=1:count_views %average for each view (same as all views, if count_views=1)
        patches_btc_viewsel=patches_btc{iNR}(:,patches_view{iNR}==iview); %patches in this view
        btc_avg_v(iNR,:,iview,:)=repmat(mean(patches_btc_viewsel',1),[1 1 1 btc_avg_types]);
        btc_std_v(iNR,:,iview,:)=repmat(std(patches_btc_viewsel',0,1),[1 1 1 btc_avg_types]);
    end
    btc_avg_allv(iNR,:,:)=repmat(mean(patches_btc{iNR}',1),[1 1 btc_avg_types]);
    btc_std_allv(iNR,:,:)=repmat(std(patches_btc{iNR}',0,1),[1 1 btc_avg_types]);
end
%
ffdm_btc_plot_gen;
