%ffdm_btc_plot_gen: basic plots for binary texture coordinate analysis of ffdm datasets, and related
%
%  plots and analyzes local image statistics of patches and their covariances
%  -- calculated in ffdm_btc_calc_gen (also whitened there)
%
% run this after ffdm_btc_calc_gen, or on saved workspace.
% analysis parameters inherited from ffdm_btc_calc_gen.
%
% for natural image stats, ffdm_btc_plot_nistats will plot in similar format
%
% data for mean and stdv plots are in btc_[avg|std]_v and btc_[avg|std]_v_jsem
%    dim 1 is iNR
%    dim 2 is btc coordinate, assumed in JV standard order
%    dim 3 is view
%    dim 4 is average type (1: per sample, 2: per subject)
%  btc_[avg|std]_allv(_jsem) has corresponding quantities averaged across
%  all views, and btc_[avg|std]_allv(_jsem)_reorg does this for L<->R flips of views, if ffdm_btc_calc_gen_reorg has been run
%
% data to be plotted for covariance plots are in patches_btc{iNR}(ibtc,ipatch) or patches_btc_reorg{iNR}(ibtc,ipatch)
%   covariances are computed here; for individual views, second coord of patches_btc is
%   selected via find(patches_view{iNR}==iview)
%    
%  See also:  FFDM_BTC_SCAT_GEN, FFDM_EXTRACT_DEMO, PREWHITEN_MRI_DEMO, BTC_DEFINE,
%  GLIDER_MAPUBI, GETCORRS_P2X2, BTC_CORRS2VEC, BTC_VEC2LETCODE, FFDM_BTC_CALC_GEN, BTC_VFLIP, 
%  NISTATS_PLOT_DEMO, FFDM_BTC_REORG, FFDM_BTC_CALC_GEN_REORG, FFDM_NRPLOT, FFDM_BTC_PLOT_NISTATS,
%  GET_NRLIST_AVAIL, GET_NRVIEW, FFDM_UPMC_NVMS_CALC, FFDM_BTC_S2NR.
%
% much code and formatting from nistats_plot_demo
%
% 21Jan20:  begin to add plotting via ffdm_nrplot
% 04Feb20:  used as back end of ffdm_btc_plot_nistats
% 05Feb20:  allow option to not plot "all views" if there is only one view
% 11Feb20:  use get_nrlist_avail to get NR_ptrs and get_nrview for useful 3D views
% 14Feb20:  fix bug so that all NR plots are labelled
% 26Oct21:  misc documentation fixes, add tstring_have for compatibility with ffdm_upmc_nvms_calc
% 27Oct21:  eliminate un-needed dimension on btc_[avg|std]_allv
% 28Oct21:  add btc_order_JV and btc_order_HX, preprocess_downsample if not already defined
% 18Aug22:  option to call ffdm_btc_s2nr to take s and n into account
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
btc_dict=btc_define;
if ~exist('dp') dp=0.15; end %position offset for horizontal axis
btc_posits_all=[1 2-dp 2+dp 3-dp 3+dp 4.5-3*dp 4.5-dp 4.5+dp 4.5+3*dp 6];
if ~exist('lw_var') lw_var=1.5; end %line width for std dev plots
if ~exist('lw_cov') lw_cov=1.5; end; %line width for covariance plots
if ~exist('np_cov') np_cov=48; end; %number of points in a covariance ellipse
if ~exist('cov_scale') cov_scale=0.9; end %size of ellipses
if ~exist('btc_marker_size') btc_marker_size=8; end
if ~exist('eb_width') eb_width=dp; end %size of error bar cap
%
if ~exist('data_source')
    data_source=1;
end
if ~exist('tstring_have')
    tstring_have=0;
end
%
if ~exist('avg_scale') avg_scale=[0 0]; end
avg_scale=getinp('scale for plotting average btc stats (0 to autoscale, [-0.05 0.05] suggested for fixed scale)','f',[-1 1],avg_scale);
if ~exist('std_scale') std_scale=[0 0]; end
std_scale=getinp('scale for plotting standard dev of btc stats (0 to autoscale, [0 0.2] suggested for fixed scale)','f',[0 1],std_scale);
if ~exist('eb_mult') eb_mult=2; end
eb_mult=getinp('error bar, as multiple of sem (0 for none)','f',[0 4],eb_mult);
%
%reorganize to choose the vert-mirror-flipped btc coordinates for the R view
%since these are single-vie statistics, they can be calculated by just permuting btc indices
%
%btc_avg_vs_reorg=ffdm_btc_reorg(btc_avg_vs,view_names,'btc_avg_vs'); %removed 04Feb20, not used below
btc_avg_v_reorg=ffdm_btc_reorg(btc_avg_v,view_names,'btc_avg_v');
btc_std_v_reorg=ffdm_btc_reorg(btc_std_v,view_names,'btc_std_v');
btc_avg_v_jsem_reorg=ffdm_btc_reorg(btc_avg_v_jsem,view_names,'btc_avg_v_jsem');
btc_std_v_jsem_reorg=ffdm_btc_reorg(btc_std_v_jsem,view_names,'btc_std_v_jsem');
%
% if ffdm_btc_calc_gen_reorg has been run, it will set this to 1
if ~exist('if_allv_reorg')
    if_allv_reorg=0;
end
if if_allv_reorg==0
    disp('Only one view, or no L<->R reorganization of views.  Averages across views will only be calculated based on native coordinates.');
else
    disp('Averages across views will be calculated without and with coordinate reorganization for L<->R flip');
end
if_plot_allv=getinp('1 to plot ''all views'' (redundant if there is only one view)','d',[0 1],double(count_views>1));
btc_orders=struct();
%set up coordinate order
if ~exist('btc_eLife_order')
    btc_eLife_order='gcbdevwtua'; %from NIstats_plot_demo, should be same as btc_order_HX with g in front
end
if ~exist('btc_order_HX')
    btc_order_HX='cbdevwtua'; %parameter order used in Hermundstad et al. and Xu et al.
end
if ~exist('btc_order_JV')
    btc_order_JV=btc_dict.codel; %gbcdetuvwa
end
btc_orders.eLife.codel=btc_eLife_order;
btc_orders.eLife.label='Hermundstad and Xu data files';
btc_orders.HX.codel=btc_order_HX;
btc_orders.HX.label='Hermundstad and Xu plot order';
btc_orders.JV.codel=btc_dict.codel;
btc_orders.JV.label='Victor et al. standard';
btc_orders_fields=fieldnames(btc_orders);
%
%groups for plotting as function of N and R
btc_plot_groups_list={{'b','c','d','e'},{'bc','de','tuvw','a'}};
btc_plot_colors_list={{[1 0 0],[0 1 0],[0.35 0 1],[0 0.35 1]},{[1 0 0],[0 1 0],[0.35 0.35 1],[0 0 0]}};
btc_plot_strings={'bcde','bc de tuvw a'};
%
if ~strcmp(btc_eLife_order,cat(2,'g',btc_orders.HX.codel)) %check that the eLife ordering in NIstats_plot_demo is same as "Hermundstad-Xu" ordering in ffm_btc_calc_gen
    warning('coordinate order mismatch.');
end
%
for iord=1:length(btc_orders_fields)
    btc_orders_field=btc_orders_fields{iord};
    disp(sprintf(' order %2.0f->%13s:  %30s  [%10s]',iord,btc_orders_field,btc_orders.(btc_orders_field).label,btc_orders.(btc_orders_field).codel));
end
iord=getinp('choice','d',[1 length(btc_orders_fields)]);
%
btc_order_field=btc_orders_fields{iord};
btc_order_codel=btc_orders.(btc_order_field).codel;
btc_order_label=btc_orders.(btc_order_field).label;
btc_n_plot=length(btc_order_codel);
btc_order_source=zeros(1,btc_n_plot);
%btc_order_source is the position in btc_[avg|std]_... on dim 3 of each coordinate to plot
for ilet=1:btc_n_plot
    btc_order_source(ilet)=find(btc_dict.codel==btc_order_codel(ilet));%this works because reordering does not change the grouping of btc coords {b,c}, {d,e},{t,u,v,w}
end
btc_posits=sort(btc_posits_all(btc_order_source));
%
for iNR=1:length(N_list)
    disp(sprintf(' %2.0f-> N(downsampling)=%2.0f R=%3.0f NR=%4.0f',iNR,N_list(iNR),R_list(iNR),NR_list(iNR)));
end
%select values of N and R, and set up color code
%NR_ptrs=getinp('selection for plotting','d',[1 length(N_list)],[1:length(N_list)]);
if ~exist('preprocess_downsample') %added 28Oct21
    preprocess_downsample=1;
end
[nrlist_avail,NR_ptrs]=get_nrlist_avail(N_list,R_list,preprocess_downsample);
nr_view=get_nrview('ask');
%N colors range from green to blue, depending on N
%R colors get darker with increasing R
N_min=min(N_list(NR_ptrs));
N_max=max(N_list(NR_ptrs));
R_min=min(R_list(NR_ptrs));
R_max=max(R_list(NR_ptrs));
if ~exist('color_loN') color_loN=[0.5 1.0 0.5]; end %light green
if ~exist('color_hiN') color_hiN=[0.5 0.5 1.0]; end %light blue
NR_colors=zeros(length(NR_ptrs),3);
for iNR_ptr=1:length(NR_ptrs)
    iNR=NR_ptrs(iNR_ptr);
    N=N_list(iNR);
    R=R_list(iNR);
    if (N_min==N_max)
        N_mix=1;
    else
        N_mix=1-(log(N/N_min)/log(N_max/N_min));
    end
    color_N=N_mix*color_loN+(1-N_mix)*color_hiN;
    if (R_min==R_max)
        R_mix=1;
    else
        R_mix=1-(log(R/R_min)/log(R_max/R_min));
    end
    color_NRmax=color_N;
    color_NRmax(find(color_NRmax==min(color_NRmax)))=0; %set this color to be a dark version of color_N
    NR_colors(iNR_ptr,:)=R_mix*color_N+(1-R_mix)*color_NRmax; %color_N if R=R_min, color_NRmax if R=R_max
    disp(sprintf(' trace %2.0f: iNR=%3.0f, (N,R)=(%2.0f,%3.0f), color(rgb) is =[%5.3f %5.3f %5.3f]',iNR_ptr,iNR,N,R,NR_colors(iNR_ptr,:)));
end
%
mean_std_labels={'mean','std dev'};
ylims=zeros(length(mean_std_labels),2);
ylims(1,:)=[min(btc_avg_v(:)) max(btc_avg_v(:))];
ylims(2,:)=[min(btc_std_v(:)) max(btc_std_v(:))];
if any(avg_scale~=0) ylims(1,:)=avg_scale; end
if any(std_scale~=0) ylims(2,:)=std_scale; end
for iview=1:count_views+if_plot_allv*(1+if_allv_reorg)
    if (iview<=count_views)
        view_string=view_names{iview};
    elseif (iview==count_views+1)
        view_string='all views';
    else
        view_string='all views, reorg';
    end
    hf_mean_std=figure;
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',cat(2,'Single-Axis Analyses (btc mean and std): ',view_string));
    set(gcf,'Position',[50 50 1200 700]);
    %
    if tstring_have==0 %0 if from ffdm_btc_calc_gen, 1 if from ffdm_upc_nvms_calc
        tstring=sprintf('%s %s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f gen',...
            datasource_string,view_string,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
            if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell);
        if (if_norand==1)
            tstring=cat(2,tstring,' NO rand');
        end
        if (if_norand==-1)
            tstring=cat(2,tstring,' Frozen rand');
        end
    end
    vplot=cell(1,length(mean_std_labels));
    vplot_jsem=cell(1,length(mean_std_labels));
    %extract quantity to plot and its standard error of mean, and implement coordinate order choice
    %vplot{1} and vplot_jsem{1} are coordinate averages
    %vplot{2} and vplot_jsem{2} are coordinate standard devs
    %
    %patches_thisview{iNR} is an array of btc stats for each patch (dim1: JV order, dim2: patch)
    patches_thisview=cell(1,length(N_list));
    if (iview<=count_views)
        vplot{1}=reshape(btc_avg_v(:,btc_order_source,iview,:),[length(N_list),btc_n_plot,btc_avg_types]);%[iNR, btc, view, avtype]-> iNR, btc, avg type
        vplot{2}=reshape(btc_std_v(:,btc_order_source,iview,:),[length(N_list),btc_n_plot,btc_avg_types]);%[iNR, btc, view, avtype]-> iNR, btc, avg type
        vplot_jsem{1}=reshape(btc_avg_v_jsem(:,btc_order_source,iview,:),[length(N_list),btc_n_plot,btc_avg_types]);%[iNR, btc, view, avtype]-> iNR, btc, avg type
        vplot_jsem{2}=reshape(btc_std_v_jsem(:,btc_order_source,iview,:),[length(N_list),btc_n_plot,btc_avg_types]);%[iNR, btc, view, avtype]-> iNR, btc, avg type
        for iNR=1:length(N_list)
            patches_thisview{iNR}=patches_btc{iNR}(btc_order_source,find(patches_view{iNR}==iview));
        end
    elseif iview==count_views+1
        vplot{1}=btc_avg_allv(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot{2}=btc_std_allv(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot_jsem{1}=btc_avg_allv_jsem(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot_jsem{2}=btc_std_allv_jsem(:,btc_order_source,:); %was :,:), changed 27Oct21
        for iNR=1:length(N_list)
            patches_thisview{iNR}=patches_btc{iNR}(btc_order_source,:);
        end
    else
        vplot{1}=btc_avg_allv_reorg(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot{2}=btc_std_allv_reorg(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot_jsem{1}=btc_avg_allv_jsem_reorg(:,btc_order_source,:); %was :,:), changed 27Oct21
        vplot_jsem{2}=btc_std_allv_jsem_reorg(:,btc_order_source,:); %was :,:), changed 27Oct21
        for iNR=1:length(N_list)
            patches_thisview{iNR}=patches_btc_reorg{iNR}(btc_order_source,:);
        end
    end
    for ims=1:length(mean_std_labels)
        for btc_avg_type=1:btc_avg_types
            subplot(btc_avg_types,length(mean_std_labels),btc_avg_type+(ims-1)*length(mean_std_labels));
            hl=[];
            ht=[];
            for iNR_ptr=1:length(NR_ptrs)
                iNR=NR_ptrs(iNR_ptr);
                N=N_list(iNR);
                NR=NR_list(iNR);
                R=NR/N;
                hp=plot(btc_posits,vplot{ims}(iNR,:,btc_avg_type),'o','MarkerSize',btc_marker_size);
                set(hp,'Color',NR_colors(iNR_ptr,:));
                set(hp,'LineWidth',lw_var);
                hold on;
                hl=[hl;hp];
                ht=strvcat(ht,sprintf('N=%2.0f, R=%3.0f',N,R));
                %
                if eb_mult>0
                    for ibtc=1:length(btc_posits)
                        heb=plot(repmat(btc_posits(ibtc),1,2),vplot{ims}(iNR,ibtc,btc_avg_type)+vplot_jsem{ims}(iNR,ibtc,btc_avg_type)*[-1 1]*eb_mult);
                        set(heb,'Color',NR_colors(iNR_ptr,:));
                        set(heb,'LineWidth',lw_var);
                        for ieb=1:2
                            hebcap=plot(btc_posits(ibtc)+eb_width*[-1 1]/2,repmat(vplot{ims}(iNR,ibtc,btc_avg_type)+vplot_jsem{ims}(iNR,ibtc,btc_avg_type)*(3-2*ieb)*eb_mult,1,2));
                            set(hebcap,'Color',NR_colors(iNR_ptr,:));
                            set(hebcap,'LineWidth',lw_var);                           
                        end
                    end
                end
            end
            set(gca,'YLim',ylims(ims,:));
            set(gca,'XTickLabel',btc_order_codel');
            set(gca,'XTick',btc_posits);
            set(gca,'XLim',[btc_posits(1)-0.5 btc_posits(end)+.5]);
            xlabel('btc coord');
            if (ims==1)
                plot(get(gca,'XLim'),[0 0],'k--');
            end
            %
            ylabel(mean_std_labels{ims});
            title(btc_avg_types_labels{btc_avg_type});
            %
            hleg=legend(hl,ht,'Location','Best');
            set(hleg,'FontSize',7);
        end %avg type
    end %ims mean and std
    axes('Position',[0.02,0.05,0.01,0.01]); %for text
    text(0,0,cat(2,view_string,', coordinate order: ',btc_order_label),'Interpreter','none');
    axis off;
    axes('Position',[0.02,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
    axis off;
    %
    % plot as function of N and R, this view
    %
    for ibtc_plot_groups=1:length(btc_plot_groups_list)
        btc_plot_groups=btc_plot_groups_list{ibtc_plot_groups};
        %
        opts_nrplot_def=[];
        opts_nrplot_def.view=nr_view;
        opts_nrplot_def.N_label='N_c_o_r_r';
        opts_nrplot_def.N_range=[1 32];
        opts_nrplot_def.R_range=[16 512];
        %
        opts_nrplot_def.colors=btc_plot_colors_list{ibtc_plot_groups};
        opts_nrplot_def.labels=btc_plot_groups;
        %
        hf_nrplot=figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'Single-Axis NRplot ',btc_plot_strings{ibtc_plot_groups},' (btc mean and std): ',view_string));
        set(gcf,'Position',[50 50 1200 700]);
        for ims=1:length(mean_std_labels)
            opts_nrplot_def.v_range=ylims(ims,:);
            opts_nrplot_def.v_label=mean_std_labels{ims};
            for btc_avg_type=1:btc_avg_types
                vals=NaN(length(NR_ptrs),length(btc_plot_groups),3);
                %vplot{1} and vplot_jsem{1} are coordinate averages
                %vplot{2} and vplot_jsem{2} are coordinate standard devs
                for iNR_ptr=1:length(NR_ptrs)
                    iNR=NR_ptrs(iNR_ptr);
                    for igroup=1:length(btc_plot_groups) %compute average within each group, and rms avg for sem's
                        btc_string=btc_plot_groups{igroup};
                        ngroup=length(btc_string);
                        btcval=zeros(ngroup,1);
                        btcval_jsem=zeros(ngroup,1);
                        for ivar=1:ngroup %get statistic (mean or std) and jsem for each btc coord
                            ibtc=find(btc_string(ivar)==btc_dict.codel); %which coordinate in standard order 
                            ivp=find(btc_order_source==ibtc); %position in vplot
                            btcval(ivar)=vplot{ims}(iNR,ivp,btc_avg_type);
                            btcval_jsem(ivar)=vplot_jsem{ims}(iNR,ivp,btc_avg_type);
                        end
                        vals(iNR_ptr,igroup,1)=mean(btcval); %mean
                        vals(iNR_ptr,igroup,[2 3])=reshape(mean(btcval)+eb_mult*[-1 1]*sqrt(mean(btcval_jsem.^2)),[1 1 2]); %combine jsem's as RMS
                    end %igroup
                end %iNR_ptr
                opts_nrplot=opts_nrplot_def;
                opts_nrplot.ha=subplot(btc_avg_types,length(mean_std_labels),btc_avg_type+(ims-1)*length(mean_std_labels));
                ffdm_nrplot(vals,N_list(NR_ptrs)*preprocess_downsample,R_list(NR_ptrs),opts_nrplot);
                title(cat(2,mean_std_labels{ims},':',btc_avg_types_labels{btc_avg_type}));
            end %btc_avg
        end %ims
        axes('Position',[0.02,0.01,0.01,0.01]); %for text
        text(0,0,tstring,'Interpreter','none');
        axis off;
    end %ibtc_plot_groups_list
    %
    %
    % plot covariances
    %
    circ_mtx=[cos(2*pi*[0:np_cov]'/np_cov),sin(2*pi*[0:np_cov]'/np_cov)];
    figure;
	set(gcf,'NumberTitle','off')
    set(gcf,'Position',[100 100 800 750]);
    ptitle=cat(2,'Normalized Covariances: ',view_string);
    set(gcf,'Name',ptitle);
	hlcov=[];
    htcov=[];
    for iNR_ptr=1:length(NR_ptrs)
        iNR=NR_ptrs(iNR_ptr);
        N=N_list(iNR);
        NR=NR_list(iNR);
        R=NR/N;
        %
        %covmtx=cov(NIstats.set(iset).stats);
        covmtx=cov(patches_thisview{iNR}');
        for icoord=1:btc_n_plot-1
            for jcoord=icoord+1:btc_n_plot
                covs=covmtx([icoord jcoord],[icoord jcoord]);
                xypts=circ_mtx*inv(sqrtm(covs));
                xypts=cov_scale*xypts/max(abs(xypts(:)))/2;
                hp=plot(icoord+xypts(:,1),btc_n_plot+1-jcoord+xypts(:,2),'k-','LineWidth',lw_cov);
                hold on;
                set(hp,'Color',NR_colors(iNR_ptr,:));
            end
        end
        hlcov=[hlcov;hp];
        htcov=strvcat(htcov,sprintf('N=%2.0f, R=%3.0f',N,R));
    end  
    plot([0 btc_n_plot+1],[btc_n_plot+1 0],'k');
    title(ptitle);
	hleg=legend(hlcov,htcov,'Location','NorthEast');
    set(hleg,'FontSize',7);
    set(gca,'XLim',[0 btc_n_plot+1]);
	set(gca,'XTick',[1:btc_n_plot]);
	set(gca,'XTickLabel',btc_order_codel');
    set(gca,'YLim',[0 btc_n_plot+1]);
    set(gca,'YTick',[1:btc_n_plot]);
    set(gca,'YTickLabel',flipud(btc_order_codel'));
    xlabel('btc coord');
    ylabel('btc coord');
    axis square;
    %
    axes('Position',[0.02,0.05,0.01,0.01]); %for text
    text(0,0,cat(2,view_string,', coordinate order: ',btc_order_label),'Interpreter','none');
    axis off;
    axes('Position',[0.02,0.01,0.01,0.01]); %for text
    text(0,0,tstring,'Interpreter','none');
    axis off;
end %iview
