%ffdm_btc_spec_gen: plots and analysis for spectra used in binary texture coordinate calculation
%
%  plots spectra of patches
%
% run this after ffdm_btc_calc_gen, or on saved workspace.
% analysis parameters inherited from ffdm_btc_calc_gen
%
% key variable is  ps_avg_v{iNR}(:,:,iview,ps_avg_type), already normalized for N and R by ffdm_btc_calc_gen
%   plots are in terms of cycles per original pixel, prior to downsampling
%   by N (but neglecting any preprocessing prior to ffdm_btc_calc_gen, e.g., in mlis_alglib_pilot)
% As in ffdm_btc_calc_gen, power spectra are not normalized by gray-level
% range, which is in roi_range (e.g., [0 255] for raw bmp, [0 4096] for raw tiff)
%
% Determination of slope on log-log coordinates:
%  limits are from flo_cycles_per_patch cycles per patch to fhi_Nyq fraction of the Nyquist,
%  patch size is N*R, lowest nonzero spatial frqeuency in FFT is 1/(N*R*fft2_padfactor).
%  Nyquist, in cycles per pixel, is (1/(2*N)).  
%  defaults are flo_cycles_per_patch=2 and fhi_Nyq=0.9, as in Xu et al., Frontiers 2019
%
%  21Jan20:  add plotting of spectral params
%  11Feb20:  use get_nrlist_avail to get NR_ptrs and get_nrview for useful 3D views
%  29Apr22:  add option to collapse across views
%  18Aug22:  documentation note:  analyses that differ only by S are identical, since
%    the blocking by S in ffdm_btc_calc_gen happens after power spectra are calculated.
%    so ffdm_btc_s2nr is irrelevant.
%
%  See also:  FFDM_BTC_CALC_GEN, FFDM_BTC_SCAT_GEN, FFDM_BTC_PLOT_GEN, REGRESS, FFDM_NRPLOT,
%    GET_NRLIST_AVAIL, GET_NRVIEW.

%
% ps_avg_types_labels={'per patch','per subj'}; %ways to average the power spectra across patches
% ps_avg_types=length(ps_avg_types_labels);
% ps_avg_types_symbs={'-',':'};

%
if_collapse=0; %set to collapse across views
if (count_views>1)
    if_collapse=getinp('1 to collapse across views','d',[0 1],0);
end
%
tstring=sprintf('%s pad=%1.0f whiten=%s marg=%1.0f pxlmean %1.0f roimean %1.0f cosbell %1.0f gen',...
    datasource_string,fft2_padfactor,whiten_avg_types_labels{whiten_avg_type_use},margin_width,...
    if_spec_subpxlmean,if_spec_subroimean,if_spec_cosbell);
ncols_pss=3; % heatmap, scattergram, ratios
nrows_pss=ps_avg_types;
if ~exist('flo_cycles_per_patch') flo_cycles_per_patch=2.0; end
if ~exist('fhi_Nyq') fhi_Nyq=0.9; end
flo_cycles_per_patch=getinp('lower spatial frequency limit, in cycles per patch, for determination of spectral slope','f',[0 8],flo_cycles_per_patch);
fhi_Nyq=getinp('upper spatial frequency limit, as fraction of Nyquist, for determination of spectral slope','f',[0 1],fhi_Nyq);
if ~exist('regr_alpha') regr_alpha=0.05; end
regr_alpha=getinp('confidence limit alpha-value','f',[0 1],regr_alpha);
if ~exist('if_plot_spec') if_plot_spec=1; end
if_plot_spec=getinp('1 to plot spectra','d',[0 1],if_plot_spec);
if ~exist('if_plot_params') if_plot_params=1; end
if_plot_params=getinp('1 to plot spectral params','d',[0 1],if_plot_params);
if ~exist('slope_plot_range') slope_plot_range=[1 6]; end
%
%definitions for one-dimensional slices
dsdefs=struct();
dsdefs.labels={'horiz','vert','diag samesign','diag oppsign','all'};
dsdefs.sfmult=[1 1 sqrt(2) sqrt(2) NaN];
dsdefs.colors={[1 0 0],[0 1 0],[0.35 0 1],[0 0.35 1],[0.5 0.5 0.5]};
dsdefs.markers={'.','.','x','o','.'};
regr=cell(length(N_list),count_views,ps_avg_types,length(dsdefs.labels));
sfreqs_per_pixel_plot=[1/NR_max 1/2]; %spatial frqeuency range to plot, cycles per pixel
%
[nrlist_avail,NR_ptrs]=get_nrlist_avail(N_list,R_list,preprocess_downsample);
nr_view=get_nrview('ask');
%
if (if_collapse)
    count_views=1;
    view_names={'avg'};
    view_nos=1;
    for iNR=1:length(NR_list)
        ps_avg_v{iNR}=mean(ps_avg_v{iNR},3); %average across views
    end
    disp('Power spectra from individual views collapsed to their average.')
end
%
for iNR_ptr=1:length(NR_ptrs)
    iNR=NR_ptrs(iNR_ptr);
    N=N_list(iNR);
    NR=NR_list(iNR);
    R=NR/N;
    fft_length=fft2_padfactor*R;
    fft_half=fft_length/2;
    sfreqs_per_pixel=min([0:fft_length-1],fliplr([1:fft_length]))/fft_length/N; %spatial frequencies of each bin in FFT
    sfreqs_per_pixel_lims=[flo_cycles_per_patch/NR,fhi_Nyq*(1./(2*N))]; %limits for fitting power law
    %
    dsf=cell(1,length(dsdefs.labels)); %frequencies, in cycles per pixel, for spectral estimates
    dsf_mask=zeros(fft_length,fft_length);
    dsf_mask((fft_half+1):end,:)=1; %remove duplicates beginning at Nyquist on first coord
    dsf_mask(:,fft_half+1)=1; %remove Nyquist on second coord
    dsf_mask(1,1)=1; %remove zero frequency
    %
    sf2=repmat(sfreqs_per_pixel.^2,fft_length,1);
    dsf_mat=sqrt(sf2+sf2'); %spatial frequency array
    %
    fit_mask_ok=cell(1,length(dsdefs.labels)); %this will have 1's at the frequencies that are within the range used for regression fitting
    for ids=1:length(dsdefs.labels)
        if strcmp(dsdefs.labels{ids},'all')
            dsf{ids}=dsf_mat(:);
            dsf{ids}(dsf_mask==1)=NaN;
        else
            dsf{ids}=sfreqs_per_pixel(2:fft_half)*dsdefs.sfmult(ids);
        end
        fit_mask_ok{ids}=dsf{ids}>=(sfreqs_per_pixel_lims(1)) & (dsf{ids}<=sfreqs_per_pixel_lims(2)) & (~isnan(dsf{ids}));
    end
    disp(sprintf(' N=%2.0f R=%3.0f: %9.0f values in FFT, %6.0f unique vals at nonzero non-Nyquist freqs, power law fit lims: %6.4f to %6.4f cycles/pixel',...
        N,R,length(dsf{ids}),sum(~isnan(dsf{ids})),sfreqs_per_pixel_lims));
        for ids=1:length(dsdefs.labels)
            disp(sprintf('   subset %15s:  %9.0f of %9.0f values to be used for fitting',...
                dsdefs.labels{ids},sum(fit_mask_ok{ids}),sum(~isnan(dsf{ids}))));
        end
    for iview=1:count_views
        %
        dsy=cell(1,length(dsdefs.labels),ps_avg_types); %samples from the spectrum
        for ps_avg_type=1:ps_avg_types
            ps=fftshift(ps_avg_v{iNR}(:,:,iview,ps_avg_type));
            ps_v=ps(fft_half+1+[0:fft_half-1],fft_half+1)'; %begin at  zero spatial frequency
            ps_h=ps(fft_half+1,fft_half+1+[0:fft_half-1]);
            ps_ssdiag=diag(ps);
            ps_osdiag=diag(flipud(ps));
            ps_ss=ps_ssdiag(1:fft_half);
            ps_os=ps_osdiag(1:fft_half);
            %
            dsy{1,ps_avg_type}=ps_h(2:end);
            dsy{2,ps_avg_type}=ps_v(2:end);
            dsy{3,ps_avg_type}=flipud(ps_ss(2:end))';
            dsy{4,ps_avg_type}=flipud(ps_os(2:end))';
            dsy{5,ps_avg_type}=ps_avg_v{iNR}(:,:,iview,ps_avg_type);
            dsy{5,ps_avg_type}=dsy{5,ps_avg_type}(:);
        end
        %
        tstring2=sprintf('%s (N=%1.0f R=%1.0f)',view_names{iview},N,R); 
        if (if_plot_spec)
            hs=figure;
            set(gcf,'NumberTitle','off');
            set(gcf,'Name',cat(2,'spectra: ',tstring2));
            set(gcf,'Position',[50 50 1200 750]);
        end
        for ps_avg_type=1:ps_avg_types
            tstring3=sprintf('%s (avg %s)',view_names{view_nos(iview)},ps_avg_types_labels{ps_avg_type});
            for ids=1:length(dsdefs.labels)
                regr_y=log(dsy{ids,ps_avg_type}(fit_mask_ok{ids}));
                regr_y=regr_y(:);
                regr_f=log(dsf{ids}(fit_mask_ok{ids}));
                regr_x=[ones(length(regr_y),1),regr_f(:)];
                [reg.b,reg.b_int,reg.r,reg.r_int,reg.stats]=regress(regr_y,regr_x,regr_alpha);
                %stats contains, in the following order:
                % the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
                regr{iNR,iview,ps_avg_type,ids}=reg;
            end
            %disp(tstring3);
            if (if_plot_spec)
                irow=ps_avg_type;
                icol=1; %heatmap
                subplot(nrows_pss,ncols_pss,icol+(irow-1)*ncols_pss)
                imagesc(max(log10(ps_max)-ps_log10range,log10(fftshift(ps_avg_v{iNR}(:,:,iview,ps_avg_type)))),log10(ps_max)+[-ps_log10range 0]);
                title(tstring3);
                set(gca,'XTick',0.5+[0 fft_half fft_length]);
                set(gca,'XTickLabel',[-1/(2*N) 0 1/(2*N)]);
                set(gca,'XLim',0.5+[0 fft_length]);
                set(gca,'YTick',0.5+[0 fft_half fft_length])
                set(gca,'YTickLabel',[-1/(2*N) 0 1/(2*N)]);
                set(gca,'YLim',0.5+[0 fft_length]);
                set(gca,'YDir','normal') %so that positive Y is plotted at the top, in contrast to Matlab'e standard imagesc convention
                axis square;
                %
                icol=2; %scattergram
                % horizontal axis are shown in red, vertical in green, oblique (same sign) in magenta dots, oblique (opposite sign) in magenta x
                subplot(nrows_pss,ncols_pss,icol+(irow-1)*ncols_pss)
                hl=[];
                ht=[];
                for idsrev=1:length(dsdefs.labels)
                    ids=length(dsdefs.labels)+1-idsrev; %plot in reverse order to not occlude 1d slices
                    dsy{ids,ps_avg_type}(dsf{ids}==NaN)=NaN; %don't plot zero-frequency values
                    hds=loglog(dsf{ids},dsy{ids,ps_avg_type},dsdefs.markers{ids});
                    hl=[hl;hds];
                    ht=strvcat(ht,dsdefs.labels{ids});
                    set(hds,'Color',dsdefs.colors{ids});
                    hold on;
                    %plot regression line
                    b=regr{iNR,iview,ps_avg_type,ids}.b;
                    ps_fit=exp(b(1)+b(2)*log(sfreqs_per_pixel_plot));
                    hds_reg=loglog(sfreqs_per_pixel_plot,ps_fit,'k');
                    set(hds_reg,'Color',dsdefs.colors{ids});
                end
                for ilim=1:2
                    loglog(sfreqs_per_pixel_lims(ilim)*[1 1],get(gca,'YLim'),'b');
                end
                legend(hl,ht,'Location','NorthEast');
                xlabel('cy/pixel');
                ylabel(cat(2,'power',' (avg ',ps_avg_types_labels{ps_avg_type},')'));
                set(gca,'XLim',sfreqs_per_pixel_plot);
                set(gca,'YLim',ps_max*[ps_minfac,1]);
                %
                icol=3; %compare H with V, and oblique axes
                subplot(nrows_pss,ncols_pss,icol+(irow-1)*ncols_pss)
                loglog(sfreqs_per_pixel(2:fft_half),dsy{2,ps_avg_type}./dsy{1,ps_avg_type},'k');
                hold on;
                loglog(sfreqs_per_pixel(2:fft_half)*sqrt(2),dsy{3,ps_avg_type}./dsy{4,ps_avg_type},'k:');
                loglog([1/NR_max 0.5],[1 1],'k--');
                legend({'vert/horiz','same-sign/opp-sign','unity'},'Location','SouthWest','FontSize',8);
                xlabel('cy/pixel');
                ylabel(cat(2,'power ratio',' (avg ',ps_avg_types_labels{ps_avg_type},')'));
                set(gca,'XLim',sfreqs_per_pixel_plot);
                set(gca,'YLim',10.^[-2 2]);
            end %if_plot_spec
       end %avg type
       if (if_plot_spec)
           axes('Position',[0.02,0.01,0.01,0.01]); %for text
           text(0,0,tstring,'Interpreter','none');
           axis off;
           axes('Position',[0.02,0.05,0.01,0.01]); %for text
           text(0,0,tstring2,'Interpreter','none');
           axis off;
       end %if_plot_spec     
    end %iview
end %iNR
%show fitted params in more sensible order
for ps_avg_type=1:ps_avg_types
    disp(' ');
    disp(sprintf('power spectra as average %s',ps_avg_types_labels{ps_avg_type}));
    for iview=1:count_views
        disp(sprintf('view %s: best fit of power spec = A*(sf)^h',view_names{view_nos(iview)}));
        for ids=1:length(dsdefs.labels)
            disp(dsdefs.labels{ids})
            for iNR_ptr=1:length(NR_ptrs)
                iNR=NR_ptrs(iNR_ptr);
                N=N_list(iNR);
                NR=NR_list(iNR);
                R=NR/N;
                reg=regr{iNR,iview,ps_avg_type,ids};
                disp(sprintf(' N=%2.0f R=%3.0f: A %0.7g  [%0.7g %0.7g] h %6.3f [%6.3f %6.3f]',...
                    N,R,exp(reg.b(1)),exp(reg.b_int(1,:)),reg.b(2),reg.b_int(2,:)));
            end %iNR
        end %ids
    end %iview
end %ps_avg_type
%composite plot, using ffdm_nrplot
if (if_plot_params)
    opts_nrplot_def=[];
    opts_nrplot_def.view=nr_view;
    opts_nrplot_def.N_label='N_c_o_r_r';
    opts_nrplot_def.N_range=[1 32];
    opts_nrplot_def.R_range=[16 512];
    opts_nrplot_def.v_range=slope_plot_range;
    opts_nrplot_def.v_label='h (spectral slope)';
    opts_nrplot_def.zplane_draw=2;
    for ps_avg_type=1:ps_avg_types
        %each page: one average type, a subplot for each view
        tstring_params=ps_avg_types_labels{ps_avg_type};
        figure;
        %on each subplot: spectral slopes for each view, along different axes
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'slopes by view: average ',tstring_params));
        set(gcf,'Position',[50 50 1000 800]);
        %
        [nr,nc]=nicesubp(count_views);
        for iview=1:count_views
            vals=NaN(length(NR_ptrs),length(dsdefs.labels),3);
            for ids=1:length(dsdefs.labels)
                for iNR_ptr=1:length(NR_ptrs)
                    iNR=NR_ptrs(iNR_ptr);
                    reg=regr{iNR,iview,ps_avg_type,ids};
                    vals(iNR_ptr,ids,1)=-reg.b(2);
                    vals(iNR_ptr,ids,[2:3])=-reshape(reg.b_int(2,:),[1 1 2]);
                end %iNR
            end %ids
            opts_nrplot=opts_nrplot_def;
            opts_nrplot.colors=dsdefs.colors;
            opts_nrplot.labels=dsdefs.labels;
            opts_nrplot.ha=subplot(nr,nc,iview);
            ffdm_nrplot(vals,N_list(NR_ptrs)*preprocess_downsample,R_list(NR_ptrs),opts_nrplot);
            title(view_names{view_nos(iview)});
        end
        axes('Position',[0.02,0.01,0.01,0.01]); %for text
        text(0,0,cat(2,tstring,'; average ',tstring_params),'Interpreter','none');
        axis off;
        %
        figure;
        %on each subplot: spectral slopes for axis, different views
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,'slopes by axis: average ',tstring_params));
        set(gcf,'Position',[50 50 1000 800]);
        %
        [nr,nc]=nicesubp(length(dsdefs.labels));
        for ids=1:length(dsdefs.labels)
            vals=NaN(length(NR_ptrs),count_views,3);
            for iview=1:count_views
                for iNR_ptr=1:length(NR_ptrs)
                    iNR=NR_ptrs(iNR_ptr);
                    reg=regr{iNR,iview,ps_avg_type,ids};
                    vals(iNR_ptr,iview,1)=-reg.b(2);
                    vals(iNR_ptr,iview,[2:3])=-reshape(reg.b_int(2,:),[1 1 2]);
                end %iNR
            end %ids
            opts_nrplot=opts_nrplot_def;
            opts_nrplot.colors=view_colors;
            opts_nrplot.labels=view_names;
            opts_nrplot.ha=subplot(nr,nc,ids);
            opts_nrplot_used=ffdm_nrplot(vals,N_list(NR_ptrs)*preprocess_downsample,R_list(NR_ptrs),opts_nrplot);
            title(dsdefs.labels{ids});
        end
        axes('Position',[0.02,0.01,0.01,0.01]); %for text
        text(0,0,cat(2,tstring,'; average ',tstring_params),'Interpreter','none');
        axis off;
    end %ps_avg_type
end %if_plot_params

