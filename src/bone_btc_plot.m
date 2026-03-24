%bone_btc_plot: basic plots for binary texture coordinate analysis of bone and scintigraphic images
%  analyzed by Elena Ajayi
%
% derived from ffdm_btc_plot_nistats, to use the same input format and then 
% call ffdm_btc_plot_gen to do the plotting
%
% to be run after bone_btc_demo, after clearing extraneous variables
% clear patch_fft patches_binarized patches_blockcounts patches_downsampled patches_whitened roi_data_all patches windowed im_data_all
% clear big*
%
% limitations:
%  * does not compute error bars (patches from the same ROI or same study are not tracked)
%  * only some metadata are properly captured
%
%  To fix: data_source, metadata
%
%  Source designation: here (natural images) it is 6:
%  * data_source= 1,2,3: full-field digital mammography (ffdm), for pilot data, original patches, and processed patches
%  * data_source=4: carotid ultrasound (usic, from Amanda Simon
%  * data_source=5: MRI (mrix), from Sam Xu
%
%  See also:  BTC_DEFINE, NISTATS_PLOT_DEMO, FFDM_NRPLOT, FFDM_BTC_PLOT_GEN, FFDM_BTC_PLOT_NISSTATS.
%
if ~exist('btc_dict') btc_dict=btc_define([]); end
btc_n=length(btc_dict.codel); %10
btc_eLife_order='gcbdevwtua';
if ~exist('bone_data_path') bone_data_path='./'; end
if ~exist('bone_dir') bone_dir='./'; end 
if ~exist('filename') filename='bone_standard.mat'; end
bone_data_path=getinp('data path','s',[],bone_data_path);
bone_dir=getinp('data dir within path','s',[],bone_dir);
filename=getinp('file containing data','s',[],filename);
fullname=cat(2,bone_data_path,filesep,bone_dir,filesep,filename);
fullname=strrep(fullname,'/',filesep);
fullname=strrep(fullname,cat(2,filesep,filesep),filesep);
fullname_show=strrep(fullname,filesep,'/');
bs=load(fullname);
nsetups=length(bs.patches_btc);
BSstats.N=bs.N_list;
BSstats.R=bs.R_list;
%convert output of bone_btc_demo to Hermundstad form
for isetup=1:nsetups
    BSstats.set(isetup).stats=bs.patches_btc{isetup}';
    BSstats.set(isetup).N=bs.N_list(isetup);
    BSstats.set(isetup).R=bs.R_list(isetup);
end
%capture metadata
%data source
data_source=NaN; %not yet specified
datasource_string=filename;
%
fft2_padfactor=bs.fft2_padfactor;%  as specified in bone_btc_demo (this is 1 for natural images, no padding in Hermundstad et al. code)
whiten_avg_types_labels={'per subj and view','per view (patch wtd)','per view (subj wtd)','global (patch wtd)','global (subj wtd)','per patch'};
whiten_avg_type_use=4; %global, patch-weighted, also specified in Hermundstad et al.
margin_width=bs.margin_width; %no margin .
if_spec_subpxlmean=0; %never done in bone_btc_demo
if_spec_subroimean=1; %always done in bone_btc_demo
if_spec_cosbell=bs.wintype; %from code
%if_norand=1; %no randomization mentioned
if_norand=1-2*bs.if_frozen; %if_frozen=1 ->if_norand=-1, if_frozen=0->if_norand=0
preprocess_downsample=1;
%
N_list=BSstats.N;
R_list=BSstats.R;
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
%extract data from BSstats and reformat into the quantitites
%(btc_[avg|std]*) expected by ffdm_btc_plot_gen
for iNR=1:length(N_list)
    if (BSstats.set(iNR).N ~= N_list(iNR))
        warning('mismatch of N for iNR=%3.0f',iNR);
    end
    if (BSstats.set(iNR).R ~= R_list(iNR))
        warning('mismatch of R for iNR=%3.0f',iNR);
    end
    NIdata=BSstats.set(iNR).stats;
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
