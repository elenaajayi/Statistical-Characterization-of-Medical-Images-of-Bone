function [results,opts_used]=bone_pspec_plot(pwrspec,fft2_padfactor,name_label,opts)
% [results,opts_used]=bone_pspec_plot(pwrspec,fft2_padfactor,name_label,opts)
% plots the power spectrum from quantitites computed in bone_pspec_demo and
% computes the power-law slope.  Uses similar logic to ffdm_btc_spec_gen.
%
%  Input:
%   pwrspec: 2-dimensional power spectrum, after fftshift has been applied
%   fft2_padfactor: padding factor, 1 if not supplied
%   name_label: a text name_label to appear on plots, empty if not supplied
%   opts: plotting and analysis options
%      if_plot: 1 to plot (default)
%      if_regress: 1 to do regression (default)
%      if_showregr:  1 to show regression analysis (default)
%      if_norm_patch_size: 1 to normalize so that spectra are approx independent of patch size
%      tick_fracs: positions of ticks, in cycles/pixel, relative to center (1=Nyquist)
%         defaults to [-.75:.25:.75], which will display as [-.375:.125:.375]
%      slope_colors: colors to use for slopes on scatterplot, defaults to {'b','m','r'}; 
%      slope_vals: values to use on scatterplot, defaults to [-2 -3 -4]
%      slope_thickness: thickness of slope lines, defaults to 2
%      flo_cycles_per_patch: lower spatial frequency limit, in cycles per patch, for determination of spectral slope, defaults to 2
%      fhi_Nyq: upper spatial frequency limit, as fraction of Nyquist, for determination of spectral slope, defaults to 0.9
%      regr_alpha: p-value for regression confidence limits for determination of spectral slope, defaults to 0.05
%
% Output:
%   results: results of analysis
%     results.regression: regression analysis (power-law analysis) for various selections of power spectrum
%   opts_used: options used, after defaults filled in
%
%  See also:
%
%   BONE_PSPEC_DEMO, FILLDEFAULT, REGRESS, FFDM_BTC_SPEC_GEN.
if (nargin<=1)
    fft2_padfactor=1;
end
if (nargin<=2)
    name_label='';
end
if (nargin<=3)
    opts=struct();
end
%
opts=filldefault(opts,'if_plot',1);
opts=filldefault(opts,'if_regress',1); 
opts=filldefault(opts,'if_showregr',1);
opts=filldefault(opts,'if_norm_patch_size',1);
opts=filldefault(opts,'tick_fracs',[-0.75:.25:.75]);
opts=filldefault(opts,'slope_colors',{'b','m','r'});
opts=filldefault(opts,'slope_vals',[-2 -3 -4]);
opts=filldefault(opts,'slope_thickness',2);
opts=filldefault(opts,'flo_cycles_per_patch',2.0); % lower spatial frequency limit, in cycles per patch, for determination of spectral slope
opts=filldefault(opts,'fhi_Nyq',0.9); % upper spatial frequency limit, as fraction of Nyquist, for determination of spectral slope
opts=filldefault(opts,'regr_alpha',0.05);
%
results=[];
opts_used=opts;
%
nfft=length(pwrspec);
patch_size=nfft/fft2_padfactor;
if (opts.if_norm_patch_size)
    pwrspec=pwrspec/(patch_size).^4; %patch size^2 is area, but square again for pwrspec
end
if (opts.if_plot)
    %
    % plot power spectrum as heatmap
    %
    figure;
    set(gcf,'NumberTitle','off'); %turn off numbered figure
    set(gcf,'Name',name_label);
    set(gcf,'Position',[100 100 1200 700]);
    zerofreq=nfft/2+1;
    imagesc(log10(pwrspec));
    set(gca,'YDir','normal'); %17Jul21 added so that y-axis goes from low to high
    title('power spectrum')
    tick_locs=zerofreq+(fft2_padfactor)*patch_size/2*opts.tick_fracs;
    cycles_per_patch=opts.tick_fracs*(patch_size/2); %patch_size/2 is the Nyquist freq
    cycles_per_pixel=cycles_per_patch/patch_size;
    tick_vals=cycles_per_pixel;
    set(gca,'XTick',tick_locs);
    set(gca,'XTickLabel',tick_vals);
    xlabel('cycles per pixel');
    set(gca,'YTick',tick_locs);
    set(gca,'YTickLabel',tick_vals);
    ylabel('cycles per pixel');
    axis equal;
    axis tight;
    colorbar;
    % add a name_label
    axes('Position',[0.02,0.02,0.01,0.01]); %for text
    text(0,0,name_label,'Interpreter','none');
    axis off;
    drawnow;
end %if_plot
%
%show power spectrum as function of frequency magnitude
%
%create an array of frequency magnitudes, of size equal to that of pwrspec
%
freqs=([0:nfft-1]-nfft/2)/nfft; %cycles per pixel
fsq_1d=(repmat(freqs,nfft,1)).^2; %freq squared
fsq=fsq_1d+fsq_1d'; %fx2+fy2
fmag=sqrt(fsq);
%
%for scattergrams, select based on frequencies that are less than Nyquist, i.e., 0.5 cycle/pixel, and also >0
fselect=(fmag<0.5) & (fmag>0);  
%for regression, select based on flo_cycles_per_patch and fhi_Nyq
fselect_regr=(fmag<0.5*opts.fhi_Nyq) & (fmag>=opts.flo_cycles_per_patch/patch_size);
%define the slices
dsdefs=struct();
dsdefs.labels={'horiz','vert','diag samesign','diag oppsign','all'};
dsdefs.colors={[1 0 0],[0 1 0],[0.35 0 1],[0 0.35 1],[0.5 0.5 0.5]};
dsdefs.markers={'.','.','x','o','.'};
ndsdefs=length(dsdefs.labels);
% compute filters for each kind of slice
filt=zeros(nfft,nfft,ndsdefs); %define arrays with 1's selecting the approriate slice
filt(nfft/2+1,:,1)=1; %horizontal
filt(:,nfft/2+1,2)=1; %vertical
filt(:,:,3)=eye(nfft); %same-sign diagonal
filt(:,:,4)=fliplr(eye(nfft)); %opposite-sign diagonal
filt(:,:,5)=ones(nfft); %all
%
if (opts.if_plot)
    figure; %scattergram of power spectrum as function of frequency
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',name_label);
    set(gcf,'Position',[100 100 1200 700]);
    %
    subplot(1,2,1);
    loglog(fmag(fselect(:)),pwrspec(fselect(:)),'k.');
    xlabel('frequency, cycles/pixel');
    ylabel('power spectrum');
    hold on;
    legend_text={'power spectrum'};
    %
    %add lines of specific slopes
    %
    if length(opts.slope_vals)>0
        %ymid=geomean(get(gca,'YLim'))/10.^[mean(opts.slope_vals)/2]; %get a typical y-value and center
        ymid=exp(mean(log(get(gca,'YLim'))))/10.^[mean(opts.slope_vals)/2]; %get a typical y-value and center
        xslope=[.01 .1]; %x-values plotted for standard slopes
        for islope=1:length(opts.slope_vals)
            hp=plot(xslope,[ymid ymid*10.^(opts.slope_vals(islope))]);
            set(hp,'Color',opts.slope_colors{1+mod(islope-1,length(opts.slope_colors))});
            set(hp,'LineWidth',opts.slope_thickness);
            legend_text{end+1}=sprintf('slope = %3.1f',opts.slope_vals(islope));
        end
    end
    xlims=get(gca,'XLim'); %get limits of panel so same limits are used for other plot
    ylims=get(gca,'YLim');
    legend(legend_text,'Location','NorthEast','Interpreter','none');
    %
    % scattergram of specific slices
    %
    subplot(1,2,2); 
    %
    lh=[];
    for idsrev=1:ndsdefs %plot the slices in reverse order so that "all" does not cover everything
        ids=ndsdefs+1-idsrev;
        slice_select=filt(:,:,ids) & fselect;
        hds=loglog(fmag(slice_select(:)),pwrspec(slice_select(:)),dsdefs.markers{ids});
        lh=[lh;hds];
        set(hds,'Color',dsdefs.colors{ids});
        hold on;
    end
    xlabel('frequency, cycles/pixel');
    ylabel('power spectrum');
    set(gca,'XLim',xlims);
    set(gca,'YLim',ylims);
    % 
end
%
%do regression
%
results.regression=cell(ndsdefs,1);
if (opts.if_regress)
    for ids=1:ndsdefs %plot the slices in reverse order so that "all" does not cover everything
        slice_select=filt(:,:,ids) & fselect_regr;
        regr_y=log(pwrspec(slice_select(:)));
        regr_f=log(fmag(slice_select(:)));
        regr_x=[ones(length(regr_y),1),regr_f(:)];
        reg=struct();
        reg.label=dsdefs.labels{ids};
        [reg.b,reg.b_int,reg.r,reg.r_int,reg.stats]=regress(regr_y,regr_x,opts.regr_alpha);
        %stats contains, in the following order:
        % the R-square statistic, the F statistic and p value for the full model, and an estimate of the error variance.
        results.regression{ids}=reg; %save regression 
        %show key regression results
        if (opts.if_showregr)
            disp(sprintf('slope for %20s: %7.3f confidence limit: [%7.3f %7.3f]',dsdefs.labels{ids},reg.b(2),reg.b_int(2,:)));
        end
    end
end
if (opts.if_plot) %show regression limits and add label
    plot(repmat(opts.flo_cycles_per_patch/patch_size,1,2),ylims,'k--');
    plot(repmat(opts.fhi_Nyq*0.5,1,2),ylims,'k--');
    legend(lh,fliplr(dsdefs.labels));
    axes('Position',[0.02,0.02,0.01,0.01]); %for text
    text(0,0,name_label,'Interpreter','none');
    axis off;
    drawnow;
end %if_plot
return
end

