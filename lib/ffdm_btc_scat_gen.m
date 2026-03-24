%ffdm_btc_scat_gen: scattergrams and correlations with clinical data for binary texture coordinate analysis
%
% run this after ffdm_btc_calc_gen, or on saved workspace.
% analysis parameters inherited from ffdm_btc_calc_gen
%
%   11Apr18: added stdv (i.e., s.e.m.) and confidence limits by jackknifing per subject, after Fisher transform (atanh)
%   18Apr18: allowed for variable N*R
%   30Apr18: added heatmaps and PCA
%   12May18: added PCA of mean and diff of L and R views
%   13May18: added correlation with clinical data, and force PCA's to be positive if projected against pca_upos=[0 1 1 1 1 0 0 0 0 0]
%   15May18: added coloring of pc weights in scattergrams according to clinical data
%   17May18: added 3d plot of above
%   03Dec19: begin options use other image databases, from fftm_btc_calc_gen
%      data_source=1->original "pilot" ffdm database; they are cut into square ROI's here
%      data_source=2->image patches cut into square ROI's by CKA
%      data_source=3->image patches stats manipulated by mlis_alglib_pilot or similar, in a "sdb" database
%   16Dec19: changes in calculation of jackknife stats because of missing data, i.e., ,length of drop-one lists may depend on view
%   17Jan20: fix bug related to display of pca vs VGF for data_source=1
%   20Jan20: bring in metadata for definitive ffdm dataset (data_source=2 or 3), prevent jackknife stats for correlation coefs if less than 4 samples
%   10Aug20: allow for subtraction of mean of statistics, default is only to do so if there is no whitening
%   12Aug20: add labeling corresponding to image statistic scale (S) from ffdm_btc_calc_gen.m
%   18Aug22:  option to call ffdm_btc_s2nr to take s and n into account -- only affects labeling of plots and output log
%
%  scattergram of btc coords in one view vs another will be empty if the two views are not available in the same subject
%  correlations with clinicla variables (BIRADS, VGF) will be empty if clinical data are not available
%
% Need to fix correlations within views for missing views
% May need to fix correlations with clin for missing views
%
%  See also:  FFDM_EXTRACT_DEMO, PREWHITEN_MRI_DEMO, BTC_DEFINE,
%  GLIDER_MAPUBI, GETCORRS_P2X2, BTC_CORRS2VEC, BTC_VEC2LETCODE, FFDM_BTC_CALC_GEN, BTC_VFLIP,
%  FFDM_BTC_SCAT, JACK, FFDM_GETCLIN, FFDM_BTC_PLOT_GEN, FFDM_BTC_S2NR.
%
if_s2nr=0; %section added 18Aug22
if ~exist('S_list')
    S_list=ones(size(N_list));
end
N_list_orig=N_list;
R_list_orig=R_list;
S_list_orig=S_list;
if any(S_list~=1)
    disp('Some values of S_list are not 1.');
    if_s2nr=getinp('1 to replace S_list, R_list, N_list by equivalents for S=1','d',[0 1]);
    if (if_s2nr)
        [N_list,R_list,S_list]=ffdm_btc_s2nr(N_list_orig,R_list_orig,S_list_orig);
        datasource_string=cat(2,datasource_string,'+s2nr');
    end
end
%
if ~exist('scat_scales')
    scat_scales=[0.1 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.1 0.1];
end
if ~exist('gray_min')
    gray_min=0.25; %minimum gray value for scatter plots -- set to 1 so that all are black
end
if ~exist('coord_pairs')
    coord_pairs={'bc','bd','be','de','cd','ce'};
end
if ~exist('pca_nplot') %number of principal components to plot
    pca_nplot=5;
end
if ~exist('pca_nkeep') %number of principal components to keep for corrrelation analysis
    pca_nkeep=5;
end
if ~exist('color_pcs')
    color_pcs='rmbcgy';
end
if ~exist('clin_colors') 
    clin_colors=[1 0 0;0 1 0];
end
if ~exist('data_source')
    data_source=1;
end
%for backward compatibility of S_list doesn't exist
if ~exist('S_list')
    S_list=ones(1,length(NR_list));
end
if ~exist('Sstep_list')
    Sstep_list=zeros(1,length(NR_list));
end
Sstep_char={' ','+'}; %+ will indicate that phases are stepped
%
data_source_strings={'ffdm: pilot','ffdm: unprocessed','ffdm: processed','usic: carotid ultrasound','mrix: brain MRI'};
%
switch data_source
    case 1
        clin_data_corr_labels={'BIRADS','VGF_L','VGF_R','VGF_(L+R)/2','VGF_L-R'};
        if ~exist('clin_filename') clin_filename=[]; end
        if ~exist('clin_columns') clin_columns=[]; end %location of clinical data for data_source=1 clin_columns=[];
        [clin_data,clin_read_info]=ffdm_getclin(clin_filename,clin_columns);
        disp('clin_read_info');
        disp(clin_read_info);
        clin_data_corr=[clin_data.BIRADS,clin_data.VGF_L,clin_data.VGF_R];
        clin_data_corr=clin_data_corr(1:count_subjs,:);
        clin_data_corr(:,[4 5])=clin_data_corr(:,[2:3])*[0.5 1;0.5 -1];
        tstring_clindata='Read clinical data.';
        reorg_string='(L)'; %to indicate that btc coordinates from R views are reorganized to correspond to L
    case {2,3}
        clin_data_corr_labels={'BIRADS'};
        %retrieve metdata from aeach suject; subject numbers are the unique values of patch_PID_nums
        PID_nums=unique(patch_PID_nums);
        clin_data_corr=zeros(count_subjs,1);
        for isubj=1:count_subjs
            PID_num=PID_nums(isubj);
            patches=find(patch_metadata.PID_num==PID_num); %all patches for this patient
            if isempty(patches)
                error(sprintf(' cannot find any patches for subj %4.0f -> PID %4.0f',isubj,PID_num));
            else
                BIRADS_vals=patch_metadata.BIRADS(patches);
                if min(BIRADS_vals)==max(BIRADS_vals)
                    clin_data_corr(isubj,1)=BIRADS_vals(1);
                    disp(sprintf(' subj %2.0f -> PID %4.0f: %5.0f patches found (%3.0f to %3.0f), BIRADS for all is %2.0f',...
                        isubj,PID_num,length(patches),min(patches),max(patches),BIRADS_vals(1)));
                else
                    error(sprintf(' inconsistent BIRADS values for subj %4.0f -> PID %4.0f, %2.0f to %2.0f',isubj,PID_num,min(BIRADS_vals),max(BIRADS_vals)));
                end
            end
        end
        tstring_clindata='Clinical data (BIRADS only) extracted from patch_metadata.';
        reorg_string='(L)'; %to indicate that btc coordinates from R views are reorganized to correspond to L
    case {4,5}
        clin_data_corr_labels=cell(0);
        tstring_clindata='No clinical data available.  All clinical metadata set to 1';
        reorg_string=''; %to indicate that there is no reorganization of btc coordinates 
    otherwise
        error(sprintf('data_source %3.0f unrecognized.',data_source));
end
clin_params=length(clin_data_corr_labels);
have_clin=double(clin_params>0);
data_source_string=data_source_strings{data_source};
%
disp(' ');
disp(sprintf(' data source: %s',data_source_string));
disp(tstring_clindata);
disp(' ');
disp('views');
disp(view_names);
disp(' ');
if_ffdm=double(~isempty(strmatch('ffdm',data_source_string)));
if (if_ffdm)
    corr_view_xlab={'LCC','LMLO','avCC'};
    corr_view_ylab={'RCC','RMLO','avMLO'};
    pca_views={1,2,3,4,[1 3],[1 3],[2 4],[2 4]};
    pca_views_howcomb=[0 0 0 0 1 -1 1 -1];
    %view_names: 'LCC', 'LMLO','RCC','RMLO' (LCC and LMLO  should be displayed against VGF_L (clin col 2), RCC and RMLO should be displayed against VGF_R (clin col 3), 
    clin_cols=[2 2 3 3 4 5 4 5]; %VGF: L, L, R, R, L+R, L-R, L+R, L-R, fixed 17Jan20, was clin_cols=[2 3 2 3 4 5 4 5] 
else %usic, mrix:  no correlations across views
    corr_view_xlab=cell(0);
    corr_view_ylab=cell(0);
    pca_views={1};
    pca_views_howcomb=0;
    clin_cols=0; %no clinical correlations
end
%
%reorganize to choose the vert-mirror-flipped btc coordinates for the R view
btc_avg_vs_reorg=zeros(length(N_list),btc_n,count_views,count_subjs);
for let=1:btc_n
    codel(let)=btc_order_JV(let);
    codel_flip(let)=btc_vflip(codel(let),btc_dict);
    let_flip=find(btc_order_JV==codel_flip(let));
    for iview=1:count_views
        let_use=let;
        if (view_names{view_nos(iview)}(1)=='R')
            let_use=let_flip;
        end
        btc_avg_vs_reorg(:,let,iview,:)=btc_avg_vs(:,let_use,iview,:);
        if let~=let_use
            disp(sprintf(' view %5s: coord, when analyzed here, %s taken from original coord %s',view_names{iview},codel(let),btc_order_JV(let_use)));
        end
    end %iview
end %let
disp(' ');
if (count_views>1)
    % set up scattergrams of a single coord between views, and determine how much data is available
    corr_view=zeros(length(N_list),length(corr_view_xlab),btc_n);
    corr_view_atanh_jsem=zeros(length(N_list),length(corr_view_xlab),btc_n);
    corr_view_pval=zeros(length(N_list),length(corr_view_xlab),btc_n);
    %
    nrows=length(corr_view_xlab); %LCC vs RCC, LMLO vs RMLO, av(CC) vs av (MLO)
    ncols=btc_n;
    %
    corr_view_xsel=zeros(length(corr_view_xlab),count_views);
    corr_view_ysel=zeros(length(corr_view_xlab),count_views);
    for irow=1:2
        corr_view_xsel(irow,strmatch(corr_view_xlab{irow},view_names,'exact'))=1;
        corr_view_ysel(irow,strmatch(corr_view_ylab{irow},view_names,'exact'))=1;
    end
    x_havedata=double((patch_count*corr_view_xsel')>0);
    y_havedata=double((patch_count*corr_view_ysel')>0);
    corr_view_xsel(3,:)=(corr_view_xsel(1,:)+corr_view_ysel(1,:))/2;
    corr_view_ysel(3,:)=(corr_view_xsel(2,:)+corr_view_ysel(2,:))/2;
    havedata_viewcorr=cell(1,nrows);
    havedata_viewcorr{1}=intersect(find(x_havedata(:,1)>0),find(y_havedata(:,1)>0));
    havedata_viewcorr{2}=intersect(find(x_havedata(:,2)>0),find(y_havedata(:,2)>0));
    havedata_viewcorr{3}=intersect(havedata_viewcorr{1},havedata_viewcorr{2});
    ndata_viewcorr=zeros(1,nrows);
    for irow=1:nrows
        ndata_viewcorr(irow)=length(havedata_viewcorr{irow});
        disp(sprintf('scatter plots between views:  %6s vs %6s has data from %4.0f subjects of %4.0f',...
            corr_view_xlab{irow},corr_view_ylab{irow},ndata_viewcorr(irow),size(patch_count,1)));
    end
    disp(' ');
else
    ndata_viewcorr=0;
end %count_views
%
%set up correlations between coords in a single view, and determine how much data is available
%
havedata_cpair=cell(1,count_views+1);
for iview=1:count_views
    havedata_cpair{iview}=find(patch_count(:,iview)>0);
end
havedata_cpair{count_views+1}=find(all(patch_count>0,2));
for iview=1:count_views+1
    ndata_cpair(iview)=length(havedata_cpair{iview});
    if (iview<=count_views)
        disp(sprintf('scatter plots within views:  %10s has data from %4.0f subjects of %4.0f',...
            view_names{iview},ndata_cpair(iview),size(patch_count,1)));
    else
        disp(sprintf('scatter plots within views:  all views have data from %4.0f subjects of %4.0f',...
            ndata_cpair(iview),size(patch_count,1)));
    end
end
%
if ~exist('whiten_avg_type_use') | ~exist('whiten_avg_types_labels')
    if_submean=0;
else
    if_submean=double(strcmp(whiten_avg_types_labels(whiten_avg_type_use),'none'));
    if_submean=getinp('1 to subtract mean of btc statistics prior to further processing','d',[0 1],if_submean);
end
% btc_avg_vs_reorg: dim 1 is iNR, dim 2 is btc coord, dim 3 is view, dim 4 is subj
% patch_count: dim 1 is views, dim 2 is subj
if (if_submean) %subtract mean, weighting all subjects equally (as long as at least one patch)
    for iNR=1:length(N_list)
        for iview=1:size(patch_count,2)
            if sum(patch_count(:,iview)>0)
                N=N_list(iNR);
                NR=NR_list(iNR);
                R=NR/N;
                S=S_list(iNR);
                Splus=Sstep_char{Sstep_list(iNR)+1};
                btc_mean=mean(btc_avg_vs_reorg(iNR,:,iview,patch_count(:,iview)>0),4);
                btc_avg_vs_reorg(iNR,:,iview,:)=btc_avg_vs_reorg(iNR,:,iview,:)-repmat(btc_mean,[1 1 1 size(btc_avg_vs_reorg,4)]);
                %
                disp(sprintf('for N=%1.0f R=%1.0f S=%1.0f%s view %s, mean subtracted from %4.0f subjects:',N,R,S,Splus,view_names{iview},sum(patch_count(:,iview)>0)))
                disp(btc_mean);
            end %sum(patch_count)
        end
    end
end
%
getinp('1 to proceed','d',[0 1],1);
%
% calculate and plot scattergrams of a single coord between views
%
if any(ndata_viewcorr>0)
    for iNR=1:length(N_list)
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        %
        tstring=sprintf('%s N=%1.0f R=%1.0f S=%1.0f%s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f submean %1.0f gen',...
            datasource_string,N,R,S,Splus,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
            if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell,if_submean);
        if (if_norand==1)
            tstring=cat(2,tstring,' NO rand');
        end
        if (if_norand==-1)
            tstring=cat(2,tstring,' Frozen rand');
        end
        % btc_avg_vs=zeros(length(N_list),btc_n,count_views,count_subjs);
        % btc_order_JV(icoord)
        hf_scat=figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'scattergrams betw views: ',tstring));
        set(gcf,'Position',[50 50 1400 750]);
        for icol=1:ncols
            %choose the flipped coordinate on appropriate views
            btc_data=reshape(btc_avg_vs_reorg(iNR,icol,:,:),count_views,count_subjs);
            btc_data(isnan(btc_data))=0; %eliminate these with havedata_pca
            for irow=1:nrows
                subplot(nrows,ncols,icol+(irow-1)*ncols);
                patch_count_view=patch_count*(sign(corr_view_xsel(irow,:))+sign(corr_view_ysel(irow,:)))';
                patch_count_min=min(patch_count_view(:));
                patch_count_max=max(patch_count_view(:));
                xcv=zeros(ndata_viewcorr(irow),1);
                ycv=zeros(ndata_viewcorr(irow),1);
                for isubj_ptr=1:ndata_viewcorr(irow)
                    isubj=havedata_viewcorr{irow}(isubj_ptr);
                    patch_frac=(patch_count_view(isubj)-patch_count_min)/(patch_count_max-patch_count_min);
                    graylev=gray_min+(1-gray_min)*patch_frac;
                    grayRGB=(1-graylev)*[1 1 1];
                    xcv(isubj_ptr)=corr_view_xsel(irow,:)*btc_data(:,isubj);
                    ycv(isubj_ptr)=corr_view_ysel(irow,:)*btc_data(:,isubj);
                    hp=plot(xcv(isubj_ptr),ycv(isubj_ptr),'k.');
                    set(hp,'MarkerEdgeColor',grayRGB);
                    set(hp,'MarkerFaceColor',grayRGB);
                    hold on;
                end           
                %modified to take into account missing data
                if ndata_viewcorr(irow)>0
                    [corr_view(iNR,irow,icol),corr_view_pval(iNR,irow,icol)]=corr(xcv,ycv);
                    %drop-one values
                    %corr_view_drop=zeros(length(N_list),length(corr_view_xlab),btc_n,count_subjs);
                    corr_view_drop=zeros(ndata_viewcorr(irow),1);
                    if ndata_viewcorr(irow)>1
                        for isubj_ptr=1:ndata_viewcorr(irow);
                            subj_ptr_list=setdiff([1:ndata_viewcorr(irow)],isubj_ptr);
                            corr_view_drop(isubj_ptr)=corr(xcv(subj_ptr_list),ycv(subj_ptr_list));
                        end          
                        [jbias,jdebiased,jvar,corr_view_atanh_jsem(iNR,irow,icol)]=jack(atanh(corr_view(iNR,irow,icol)),atanh(corr_view_drop));
                    end
                end
                %jackknife standard error of measurement
                %
                set(gca,'XLim',scat_scales(icol)*[-1 1]);
                set(gca,'YLim',scat_scales(icol)*[-1 1]);
                plot(scat_scales(icol)*[-1 1],[0 0],'k');
                plot([0 0],scat_scales(icol)*[-1 1],'k');
                plot(scat_scales(icol)*[-1 1],scat_scales(icol)*[-1 1],'k:');
                set(gca,'XTick',scat_scales(icol)*[-1 0 1]);
                set(gca,'YTick',scat_scales(icol)*[-1 0 1]);
                set(gca,'FontSize',7);
                axis square;
                if (irow==1)
                    title(sprintf('%s (L) and %s (R)',codel(icol),codel_flip(icol)));
                end
                xlabel(corr_view_xlab{irow},'FontSize',7);
                if (icol==1)
                    ylabel(corr_view_ylab{irow},'FontSize',7);
                end
            end %icol
        end %irow
        axes('Position',[0.02,0.07,0.01,0.01]); %for text
        text(0,0,'darkness: number of ROIs','Interpreter','none');
        axis off;
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end
    %jackknife standard error of measurement
    %[jbias,jdebiased,jvar,corr_view_atanh_jsem]=jack(atanh(corr_view),atanh(corr_view_drop));
    %make a table
    disp(' ');
    disp('Correlations between views');
    disp(tstring_NR);
    disp(sprintf('Image statistics calculated with weighting %s; sems and CLs (%7.5f) by subject -based jackknife',...
        btc_avg_types_labels{btc_avg_type},cl_prob));
    corr_view_header='  N  R   S';
    corr_view_NRfmt='%3.0f %3.0f %2.0f%s';
    corr_view_statfmt='%8s   ';
    for ilet=1:btc_n
        corr_view_header=cat(2,corr_view_header,'    ',btc_dict.codel(ilet),'     pval ');
        corr_view_NRfmt=cat(2,corr_view_NRfmt,' %7.4f %6.3f');
        corr_view_statfmt=cat(2,corr_view_statfmt,' %7.4f       ');
    end
    for ileg=1:length(corr_view_xlab)
        disp(' ');
        disp(sprintf('%s vs %s',corr_view_xlab{ileg},corr_view_ylab{ileg}));
        disp(corr_view_header);
        for iNR=1:length(N_list)
            N=N_list(iNR);
            NR=NR_list(iNR);
            R=NR/N;
            S=S_list(iNR);
            Splus=Sstep_char{Sstep_list(iNR)+1};
            disp(sprintf(corr_view_NRfmt,N,R,S,Splus,squeeze(cat(2,corr_view(iNR,ileg,:),corr_view_pval(iNR,ileg,:)))));
            disp(sprintf(corr_view_statfmt,'jk_sem',tanh(corr_view_atanh_jsem(iNR,ileg,:))));
            if (cl_prob>0)
                disp(sprintf(corr_view_statfmt,'cl_lo',tanh(atanh(corr_view(iNR,ileg,:))-tmult*corr_view_atanh_jsem(iNR,ileg,:))));
                disp(sprintf(corr_view_statfmt,'cl_hi',tanh(atanh(corr_view(iNR,ileg,:))+tmult*corr_view_atanh_jsem(iNR,ileg,:))));
            end
        end 
    end %ileg
end %data available for correlations between views?
%
% calculate and plot scattergrams of a one coordinate vs another, and tabulate, with significance
%
corr_pair=zeros(length(N_list),length(coord_pairs),count_views+1);
corr_pair_atanh_jsem=NaN(length(N_list),length(coord_pairs),count_views+1);
corr_pair_pval=zeros(length(N_list),length(coord_pairs),count_views+1);
%
havedata_cpair=cell(1,count_views+1);
for iview=1:count_views
    havedata_cpair{iview}=find(patch_count(:,iview)>0);
end
havedata_cpair{count_views+1}=find(all(patch_count>0,2));
for iview=1:count_views+1
    ndata_cpair(iview)=length(havedata_cpair{iview});
    if (iview<=count_views)
        disp(sprintf('scatter plots within views:  %10s has data from %4.0f subjects of %4.0f',...
            view_names{iview},ndata_cpair(iview),size(patch_count,1)));
    else
        disp(sprintf('scatter plots within views:  all views have data from %4.0f subjects of %4.0f',...
            ndata_cpair(iview),size(patch_count,1)));
    end
end
for iNR=1:length(N_list)
    N=N_list(iNR);
    NR=NR_list(iNR);
    R=NR/N;
    S=S_list(iNR);
    Splus=Sstep_char{Sstep_list(iNR)+1};
    %
    tstring=sprintf('%s N=%1.0f R=%1.0f S=%1.0f%s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f submean %1.0f gen',...
        datasource_string,N,R,S,Splus,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
        if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell,if_submean);
    if (if_norand==1)
        tstring=cat(2,tstring,' NO rand');
    end
    if (if_norand==-1)
        tstring=cat(2,tstring,' Frozen rand');
    end
    % btc_avg_vs=zeros(length(N_list),btc_n,count_views,count_subjs);
    % btc_order_JV(icoord)
    hf_scatpair=figure;
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'scattergrams betw coords: ',tstring));
    set(gcf,'Position',[50 50 1200 750]);
    %
    [nrows,ncols]=nicesubp(length(coord_pairs),0.7);
    for ipair=1:length(coord_pairs)
        subplot(nrows,ncols,ipair)
        coord_pair=coord_pairs{ipair};
        %
        xlet=find(coord_pair(1)==btc_order_JV);
        ylet=find(coord_pair(2)==btc_order_JV);
        hl=[];
        ht=[];
        for iview=1:count_views
            if ndata_cpair(iview)>0
                xvals=squeeze(btc_avg_vs_reorg(iNR,xlet,iview,havedata_cpair{iview}));
                yvals=squeeze(btc_avg_vs_reorg(iNR,ylet,iview,havedata_cpair{iview}));
                [corr_pair(iNR,ipair,iview),corr_pair_pval(iNR,ipair,iview)]=corr(xvals(:),yvals(:));
                hp=plot(xvals,yvals,'k.');
                set(hp,'MarkerEdgeColor',view_colors{view_nos(iview)});
                set(hp,'MarkerFaceColor',view_colors{view_nos(iview)});
                hl=[hl;hp];
                ht=strvcat(ht,cat(2,'view: ',view_names{view_nos(iview)}));
                hold on;
            end
        end
        if ndata_cpair(count_views+1)>0
            xvals=squeeze(mean(btc_avg_vs_reorg(iNR,xlet,:,havedata_cpair{count_views+1}),3));
            yvals=squeeze(mean(btc_avg_vs_reorg(iNR,ylet,:,havedata_cpair{count_views+1}),3));
            [corr_pair(iNR,ipair,count_views+1),corr_pair_pval(iNR,ipair,count_views+1)]=corr(xvals(:),yvals(:));
            hp=plot(xvals,yvals,'k.');
            hl=[hl;hp];
            ht=strvcat(ht,'mean across views');
        end
        hleg=legend(hl,ht,'Location','SouthEast');
        set(hleg,'FontSize',7);
        %
        set(gca,'XLim',scat_scales(xlet)*[-1 1]);
        set(gca,'YLim',scat_scales(ylet)*[-1 1]);
        plot(scat_scales(xlet)*[-1 1],[0 0],'k');
        plot([0 0],scat_scales(ylet)*[-1 1],'k');
        scat_scale_min=min(scat_scales([xlet ylet]));
        plot(scat_scale_min*[-1 1],scat_scale_min*[-1 1],'k:');
        set(gca,'XTick',scat_scales(xlet)*[-1 0 1]);
        set(gca,'YTick',scat_scales(ylet)*[-1 0 1]);
        xlabel(sprintf('%s %s',coord_pair(1),reorg_string));
        ylabel(sprintf('%s %s',coord_pair(2),reorg_string));
        title(coord_pair);
        %drop-one values
        for iview=1:count_views+1
            %corr_pair_drop=zeros(length(N_list),length(coord_pairs),length(count_views)+1,count_subjs);
            corr_pair_drop=zeros(ndata_cpair(iview),1);
            if (ndata_cpair(iview)>=4) %need at least 4 points otherwise the drop-one values will have only 2 points, and have a correlation of 1
                for isubj_ptr=1:ndata_cpair(iview)
                subj_list=setdiff(havedata_cpair{iview},havedata_cpair{iview}(isubj_ptr));
                    if (iview<=count_views)
                        xvals=squeeze(btc_avg_vs_reorg(iNR,xlet,iview,subj_list));
                        yvals=squeeze(btc_avg_vs_reorg(iNR,ylet,iview,subj_list));
                    else
                        xvals=squeeze(mean(btc_avg_vs_reorg(iNR,xlet,:,subj_list),3));
                        yvals=squeeze(mean(btc_avg_vs_reorg(iNR,ylet,:,subj_list),3));
                    end
                    corr_pair_drop(isubj_ptr)=corr(xvals(:),yvals(:));
                end %isubj_ptr
                %jackknife standard error of measurement
                [jbias,jdebiased,jvar,corr_pair_atanh_jsem(iNR,ipair,iview)]=jack(atanh(corr_pair(iNR,ipair,iview)),atanh(corr_pair_drop));
            end %ndata_cpair{iview}>1
        end %iview
    end %ipair
    %
    axes('Position',[0.02,0.02,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
    axis off;
end %iNR
%make a table of correlations between coordinates in same view
disp(' ');
disp('Correlations between image statistics in same view');
disp(tstring_NR);
disp(sprintf('Image statistics calculated with weighting %s; sems and CLs (%7.5f) by subject -based jackknife',...
    btc_avg_types_labels{btc_avg_type},cl_prob));
corr_pair_header='  N  R   S';
corr_pair_NRfmt='%3.0f %3.0f %2.0f%s';
corr_pair_statfmt='%8s   ';
for ipair=1:size(corr_pair,2)
    corr_pair_header=cat(2,corr_pair_header,'  ',coord_pairs{ipair},'      pval ');
    corr_pair_NRfmt=cat(2,corr_pair_NRfmt,' %7.4f %6.3f');
    corr_pair_statfmt=cat(2,corr_pair_statfmt,' %7.4f       ');
end
pca_upos=[0 1 1 1 1 0 0 0 0 0]; %projection of every u on this must be positive
for iview=1:count_views+1
    disp(' ');
    if iview<=count_views
        vn=cat(2,'view: ',view_names{view_nos(iview)},' ',reorg_string);
    else
        vn=cat(2,'mean across views',' ',reorg_string);
    end
    disp(sprintf(' correlations between image statistics for: %s',vn));
    disp(corr_pair_header)
    for iNR=1:length(N_list)
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        disp(sprintf(corr_pair_NRfmt,N,R,S,Splus,[corr_pair(iNR,:,iview);corr_pair_pval(iNR,:,iview)]));
        disp(sprintf(corr_pair_statfmt,'jk_sem',tanh(corr_pair_atanh_jsem(iNR,:,iview))));
        if (cl_prob>0)
            disp(sprintf(corr_pair_statfmt,'cl_lo',tanh(atanh(corr_pair(iNR,:,iview))-tmult*corr_pair_atanh_jsem(iNR,:,iview))));
            disp(sprintf(corr_pair_statfmt,'cl_hi',tanh(atanh(corr_pair(iNR,:,iview))+tmult*corr_pair_atanh_jsem(iNR,:,iview))));
        end
    end 
end %iview
%
% do PCA for each view and each scale
%
hf_pca=cell(0);
hf_pca_scat=cell(0);
nr_pca=5;
nr_pca_scat=2;
nc_pca=length(N_list);
pca_vars=length(pca_views);
pca_wts=zeros(count_subjs,pca_nkeep,pca_vars,length(N_list));
pca_var_name=cell(0);
pca_corr_hdfmt='                 ';
pca_corr_dafmt='     %12s';
%determine number of ways to color pca plots according to clinical data
if (clin_params>1)
    n_clin_colorings=2; %use BIRADS and VGF
elseif clin_params==1
    n_clin_colorings=1;
else
    n_clin_colorings=0;
end
havedata_pca=cell(length(N_list),pca_vars);
for pca_var=1:pca_vars
    iviews=pca_views{pca_var}; 
    disp(' ');
    switch pca_views_howcomb(pca_var)
        case 0
            vn=view_names{view_nos(iviews)};
        case 1
            vn=cat(2,'(',view_names{view_nos(iviews(1))},'+',view_names{view_nos(iviews(2))},')/2');
        case -1
            vn=cat(2,view_names{view_nos(iviews(1))},'-',view_names{view_nos(iviews(2))});
    end
    pca_var_name{pca_var}=vn;
    pca_corr_hdfmt=cat(2,pca_corr_hdfmt,sprintf('%15s',vn));
    pca_corr_dafmt=cat(2,pca_corr_dafmt,'       %8.5f');
    tstring=sprintf('%s %s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f gen',...
        datasource_string,vn,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
        if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell);
    if (if_norand==1)
        tstring=cat(2,tstring,' NO rand');
    end
    if (if_norand==-1)
        tstring=cat(2,tstring,' Frozen rand');
    end
    disp(sprintf(' PCA of average btc statistic for %s',vn));
    if_fig=0;
    for iNR=1:length(N_list)
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        switch pca_views_howcomb(pca_var)
            case 0
                pca_avg_vs_data=squeeze(btc_avg_vs_reorg(iNR,:,iviews,:));
            case 1
                pca_avg_vs_data=(squeeze(btc_avg_vs_reorg(iNR,:,iviews(1),:))+squeeze(btc_avg_vs_reorg(iNR,:,iviews(2),:)))/2;
            case -1
                pca_avg_vs_data=squeeze(btc_avg_vs_reorg(iNR,:,iviews(1),:))-squeeze(btc_avg_vs_reorg(iNR,:,iviews(2),:));
        end
        havedata_pca{iNR,pca_var}=find(all(~isnan(pca_avg_vs_data),1));
        disp(sprintf(' %s N=%1.0f R=%1.0f S=%1.0f%s %5.0f sets available',vn,N,R,S,Splus,length(havedata_pca{iNR,pca_var})));
        if length(havedata_pca{iNR,pca_var})>=2
            if (if_fig==0) %open the figure if not opened yet
                hf_pca{pca_var}=figure;
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',tstring);
                set(gcf,'Position',[50 50 1200 750]);
                if_fig=1;
            end
            subplot(nr_pca,nc_pca,iNR);
            imagesc(pca_avg_vs_data',[-1 1]*max(abs(btc_avg_vs_reorg(:))));
            title(sprintf('%s N=%1.0f R=%1.0f S=%1.0f%s',vn,N,R,S,Splus));
            ylabel('subj');
            set(gca,'XTick',[1:btc_n],'FontSize',7);
            set(gca,'XTickLabel',codel');    
            %
            %do the PCA
            pca_avg_vs_data_h=pca_avg_vs_data(:,havedata_pca{iNR,pca_var});
            [u,s,v]=svd(pca_avg_vs_data_h); % pca_avg_vs_data=u*s*v', where s is diagonal, columns of u are orthonormal, columns of v are orthonormal
            %force u *[0 1 1 1 1 0 0 0 0 0] to be positive
            signs=sign(pca_upos*u);
            u=u*diag(signs);
            s=diag(signs)*s;
            %
            s_diag=diag(s);
            disp(cat(2,sprintf(' N=%2.0f R=%3.0f NR=%4.0f: ',N,R,NR),' eigenvalues',sprintf('%7.4f',abs(s_diag(:)'))));
            pca_nkeep_mins=min(pca_nkeep,length(s_diag));
            pca_wts(havedata_pca{iNR,pca_var},[1:pca_nkeep_mins],pca_var,iNR)=v(:,[1:pca_nkeep_mins])*s(1:pca_nkeep_mins,1:pca_nkeep_mins);
            disp(cat(2,'                    frac var expl',sprintf('%7.4f',(s_diag(:)'.^2)/sum(s_diag.^2))));
            %scree plot
            subplot(nr_pca,nc_pca,iNR+nc_pca);
            semilogy([1:min(btc_n,length(s_diag))],s_diag.^2/sum(s_diag.^2),'k.');
            hold on;
            for ipc=1:min(pca_nplot,length(s_diag))
                hp=semilogy(ipc,s_diag(ipc).^2/sum(s_diag.^2),'k.');
                set(hp,'MarkerFaceColor',color_pcs(mod(ipc-1,length(color_pcs))+1));
                set(hp,'MarkerEdgeColor',color_pcs(mod(ipc-1,length(color_pcs))+1));
            end
            set(gca,'XLim',[0 btc_n]);
            set(gca,'XTick',[1:btc_n]);
            set(gca,'XTickLabel',[1:btc_n],'FontSize',7);
            set(gca,'YLim',[0.001 1]);
            ylabel('frac var');
            %pc plot
            subplot(nr_pca,nc_pca,iNR+2*nc_pca);
            for ipc=1:min(pca_nplot,length(s_diag))
                hp=plot([1:btc_n],u(:,ipc),'k.-');
                set(hp,'Color',color_pcs(mod(ipc-1,length(color_pcs))+1));
                set(hp,'LineWidth',2);
                hold on;
            end
            plot([0 btc_n],[0 0],'k');
            set(gca,'XLim',[0 btc_n]);
            set(gca,'XTick',[1:btc_n],'FontSize',7);
            set(gca,'XTickLabel',codel');
            set(gca,'YLim',[-1 1]);
            ylabel('wt');
            %projections onto first two pcs, colored by BIRADS or colored by VGF
            for icolor_var=1:n_clin_colorings
                subplot(nr_pca,nc_pca,iNR+(2+icolor_var)*nc_pca);
                if (icolor_var==1)
                    clin_col=1; %BIRADS
                else
                    clin_col=clin_cols(pca_var);
                end
                for isubj=1:count_subjs
                    wt=(clin_data_corr(isubj,clin_col)-min(clin_data_corr(:,clin_col)))/(max(clin_data_corr(:,clin_col))-min(clin_data_corr(:,clin_col)));
                    wt=max(min(wt,1),0);
                    hp=plot3(pca_wts(isubj,1,pca_var,iNR),pca_wts(isubj,2,pca_var,iNR),wt,'k.');
                    %determine a color by expanding the clinical-correlation column into the available range
                    set(hp,'Color',[1-wt wt]*clin_colors);
                    hold on;
                end
                view(2);
                htt=title(cat(2,'x: pc1, y: pc2, z:',clin_data_corr_labels{clin_col}),'FontSize',7);
                set(htt,'Interpreter','none');
                xymax=max(max(abs(pca_wts(:,[1:2],pca_var,iNR))));
                set(gca,'XLim',[-1 1]*xymax);
                set(gca,'YLim',[-1 1]*xymax);
                plot([-1 1]*xymax,[0 0],'k');
                plot([0 0],[-1 1]*xymax,'k');
                axis square;
            end
        end %havedata_pca>=2?
    end %iNR
    if (if_fig)
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end
    %
    %now make 3-d scattergrams
    %
    disp(sprintf(' PCA scat of average btc statistic for %4s',vn));
    if_fig=0;
    for iNR=1:length(N_list)
        if length(havedata_pca{iNR,pca_var})>=2
            if (if_fig==0)
                hf_pca_scat{pca_var}=figure;
                set(gcf,'NumberTitle','off');
                set(gcf,'Name',tstring);
                set(gcf,'Position',[50 50 1200 750]);
                if_fig=1;
            end
            N=N_list(iNR);
            NR=NR_list(iNR);
            R=NR/N;
            S=S_list(iNR);
            Splus=Sstep_char{Sstep_list(iNR)+1};
            %projections onto first two pcs, colored by BIRADS or colored by VGF
            for icolor_var=1:n_clin_colorings
                subplot(nr_pca_scat,nc_pca,iNR+(icolor_var-1)*nc_pca);
                if (icolor_var==1)
                    clin_col=1; %BIRADS
                else
                    clin_col=clin_cols(pca_var);
                end
                for isubj=1:count_subjs
                    wt=(clin_data_corr(isubj,clin_col)-min(clin_data_corr(:,clin_col)))/(max(clin_data_corr(:,clin_col))-min(clin_data_corr(:,clin_col)));
                    wt=max(min(wt,1),0);
                    hp=plot3(pca_wts(isubj,1,pca_var,iNR),pca_wts(isubj,2,pca_var,iNR),clin_data_corr(isubj,clin_col),'k.');
                    %determine a color by expanding the clinical-correlation column into the available range
                    set(hp,'Color',[1-wt wt]*clin_colors);
                    hold on;
                end
                view(-25,18);
                xlabel('pc1','FontSize',7);
                ylabel('pc2','FontSize',7);
                if (icolor_var==1)
                    title(sprintf('%s N=%1.0f R=%1.0f S=%1.f%s',vn,N,R,S,Splus));
                end
                zlabel(clin_data_corr_labels{clin_col},'FontSize',7,'Interpreter','none');
                xymax=max(max(abs(pca_wts(:,[1:2],pca_var,iNR))));
                set(gca,'XLim',[-1 1]*xymax);
                set(gca,'YLim',[-1 1]*xymax);
                zp=[min(clin_data_corr(:,clin_col)),max(clin_data_corr(:,clin_col))];
                if min(zp)==max(zp) %in case they are all the same
                    zp=zp+[-0.1 0.1];
                end
                set(gca,'ZLim',zp);
                for imm=1:2
                    plot3([-1 1]*xymax,[0 0],[zp(imm) zp(imm)],'k');
                    plot3([0 0],[-1 1]*xymax,[zp(imm) zp(imm)],'k');
                end
                plot3([0 0],[0 0],zp,'k');
                axis vis3d
            end
        end %length(havedata_pca>=2)
        %
    end %iNR
    if (if_fig)
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end
end %pca_vars
disp(' ');
disp('correlation of principal components with clinical variables');
disp(' ');
if (have_clin)
    for iNR=1:length(N_list);
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        S=S_list(iNR);
        Splus=Sstep_char{Sstep_list(iNR)+1};
        disp(sprintf(' N=%3.0f R=%3.0f S=%1.0f%s',N,R,S,Splus));
        for ipca=1:pca_nkeep
            disp(sprintf(' principal component %2.0f',ipca));
            disp(pca_corr_hdfmt);
            [r,p]=corrcoef([clin_data_corr,squeeze(pca_wts(:,ipca,:,iNR))]);
            for icorr=1:size(clin_data_corr,2)
                disp(sprintf(pca_corr_dafmt,clin_data_corr_labels{icorr},r(icorr,size(clin_data_corr,2)+[1:pca_vars])));
                disp(sprintf(pca_corr_dafmt,'p',p(icorr,size(clin_data_corr,2)+[1:pca_vars])));
            end
        end
    end
    disp('Note that the colored scatter-plots of clinical data vs pc1 and pc2 can be rotated into 3D.');
else
    disp('[no clinical data]');
end
