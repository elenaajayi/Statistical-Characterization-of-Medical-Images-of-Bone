function opts_used=ffdm_nrplot(vals,N_list,R_list,opts)
% opts_used=ffdm_nrplot(vals,N_list,R_list,opts) plots one or more parameters as a function of 
% N (downsampling) and R (region size, in downsampled units)
%   N and R are plotted on log scales, ticks at octave intervals
%   values are plotted on either linear or log scale, as specified by
%
% vals: 3d array
%  dim 1: each (N_list(k),R_list(k)) point
%  dim 2: which parameter
%  dim 3: 1 for value, 2 for low clim, 2 for high clim
% N_list: list of values of N
% R_list: list of values of R
% opts: options
%   opts.ha: handle to axis; created if none given
%   opts.N_range: range of N-values, defaults to [1 32], should be powers of 2
%   opts.N_label: label for N-axis, defaults to 'N (downsampling)'
%   opts.N_connect: connect points that have same values of N, defaults to 1
%   opts.R_range: range of R-valyes, defaults to [16 512], should be powers of 2
%   opts.R_label: label for R-axis, defaults to 'R (region size)'
%   opts.R_connect: connect points that have same values of R, defaults to 1
%   opts.NR_connect: connect points that have same values of N*R, defaults to 1
%   opts.v_range: value range, defaults to [0 1];
%   opts.v_label: value label, defaults to 'value'
%   opts.v_scale: 'linear' or 'log', how values are plotted, defaults to 'linear'
%   opts.eb_size: error bar size (fraction of octave), defaults to 0.05
%   opts.eb_shape: [x,y] coords for error bar shape, defaults to [1 -1 -1 1 1;1 1 -1 -1 1], which yields a square
%       opts.eb_shape=[1 -1;0 0] yields a line segment parallel to N_list axis; opts.eb_shape=[1 -1 0 0 0;0 0 0 -1 1] yields cross-hairs
%   opts.colors:  colors for lines and points, if empty, all are 'k';
%       if non-empty, should be a cell array of values such as {[0 1 .4],'k'}, or a single value, 
%       or an array whose second dimension is 2, or a string of color specifiers
%       These are cycled through if fewer than size(vals,2)
%   opts.labels: cell array of strings, length =size(vals,2), as labels for legends.  If empty, then no legend
%   opts.linewidth: line width, defaults to 2
%   opts.linewidth_eb: error bar line widths, defaults to 1
%   opts.legend_location: legend location, defaults to 'Best', ignored if opts.labels is empty
%   opts.zplane_list: list of planes of z-axis values to draw
%   opts.Box, XGrid, YGrid, ZGrid, GridLineStyle are corresponding matlab
%      grid options (with N_grid corresponding to X_grid, YGrid to N_grid, ZGrid to v_grid) default to 'on',
%      Box defaults to 'on'
%      GridLineStyle defaults to ':'
%      view: 3D view params, defaults to [azimuth, elev]=[-30 15] (view(3) is [-37.5 30])
%      FontSize: defaults to 10, also is applied to labels
%
% opts_used: options used
%
% 07Feb20:  add zplane_list
%
%   See also:  FFDM_BTC_SCAT_GEN, FFDM_BTC_PLOT_GEN, FFDM_BTC_SPEC_GEN, FFDM_NRPLOT_TEST, FILLDEFAULT, 
%      GET_NRLIST_AVAIL, GET_NRVIEW.
if (nargin<=3)
    opts=[];
end
%define options that correspond to standard Matlab properties
opts_std.Box.val='on';
opts_std.N_grid.MLname='XGrid';
opts_std.N_grid.val='on';
opts_std.R_grid.MLname='YGrid';
opts_std.R_grid.val='on';
opts_std.v_grid.MLname='ZGrid';
opts_std.v_grid.val='on';
opts_std.GridLineStyle.val=':';
opts_std.view.val=[-30 15];
opts_std.FontSize.val=10;
%
opts=filldefault(opts,'ha',[]);
opts=filldefault(opts,'N_range',[1 32]);
opts=filldefault(opts,'N_label','N (downsampling)');
opts=filldefault(opts,'N_connect',1);
opts=filldefault(opts,'R_label','R (region size)');
opts=filldefault(opts,'R_range',[16 512]);
opts=filldefault(opts,'R_connect',1);
opts=filldefault(opts,'NR_connect',1);
opts=filldefault(opts,'v_range',[0 1]);
opts=filldefault(opts,'v_label','value');
opts=filldefault(opts,'v_scale','linear');
opts=filldefault(opts,'eb_size',0.05);
opts=filldefault(opts,'eb_shape',[1 -1 -1 1 1;1 1 -1 -1 1]);
opts=filldefault(opts,'colors','k');
opts=filldefault(opts,'labels',cell(0));
opts=filldefault(opts,'linewidth',2);
opts=filldefault(opts,'linewidth_eb',1);
opts=filldefault(opts,'zplane_draw',[]);
opts=filldefault(opts,'linewidth_zplane',1);
opts=filldefault(opts,'legend_location','Best');
%
std_fields=fieldnames(opts_std);
for ifield=1:length(std_fields)
    fn=std_fields{ifield};
    opts=filldefault(opts,fn,opts_std.(fn).val);
end
%
npts=size(vals,1);
nparams=size(vals,2);
if ~isempty(opts.ha)
    ha=opts.ha;
else
    ha=gca;
    opts.ha=ha;
end
%
%determine color specifiers for each parameter
%
colorspecs=cell(0);
colors_supplied=opts.colors;
if iscell(colors_supplied)
    ncolors=length(colors_supplied);
    for ic=1:ncolors
        colorspecs{ic}=colors_supplied{ic};
    end
elseif ischar(colors_supplied)
    ncolors=length(colors_supplied);
    for ic=1:ncolors
        colorspecs{ic}=colors_supplied(ic);
    end
elseif isnumeric(colors_supplied) & (size(colors_supplied,2)==3)
    ncolors=size(colors_supplied,1);
    for ic=1:ncolors
        colorspecs{ic}=colors_supplied(ic,:);
    end
else
    ncolors=1;
    colorspecs{1}='k';
end
opts.colorspecs=colorspecs;        
%
leg_labels=[];
leg_handles=[];
for iparam=1:nparams
    colorspec=colorspecs{1+mod(iparam-1,ncolors)};
    have_handle=0;
    for ipt=1:npts
        style_point='k.';
        style_line='k';
        if ~isnan(vals(ipt,iparam,1))
            N_val=log(N_list(ipt))/log(2);
            R_val=log(R_list(ipt))/log(2);
            hp=plot3(N_val,R_val,vals(ipt,iparam,1),style_point);
            set(hp,'Color',colorspec);
            hold on;
            if ~isempty(opts.labels) & have_handle==0
                leg_labels=strvcat(leg_labels,opts.labels{iparam});
                leg_handles=[leg_handles;hp];
                have_handle=1;
            end
            for icl=1:2
                if ~isnan(vals(ipt,iparam,1+icl))
                    hp=plot3(repmat(N_val,2,1),repmat(R_val,2,1),squeeze(vals(ipt,iparam,[1 1+icl])),style_line);
                    set(hp,'Color',colorspec);
                    set(hp,'LineWidth',opts.linewidth_eb);
                    hold on;
                    if opts.eb_size>0 %draw an error bar cap
                        hp=plot3(N_val+opts.eb_size*opts.eb_shape(1,:),R_val+opts.eb_size*opts.eb_shape(2,:),repmat(vals(ipt,iparam,1+icl),1,size(opts.eb_shape,2)),style_line);
                        set(hp,'Color',colorspec);
                        set(hp,'LineWidth',opts.linewidth_eb);
                    end
                end
            end %confidence limit
            if ~isempty(opts.zplane_draw)
                for iz=1:length(opts.zplane_draw)
                    hp=plot3(repmat(N_val,2,1),repmat(R_val,2,1),[vals(ipt,iparam,1) opts.zplane_draw(iz)],'k:');
                    set(hp,'Color',colorspec);
                    set(hp,'LineWidth',opts.linewidth_zplane);                   
                end
            end
        end %value present?
    end
    %optionally connect points with same N and adjacent R, same R and adjacent N, or same product NR and adjacent R/N
    for ictype=1:3
        switch ictype
            case 1
                if_connect=opts.N_connect;
                match_list=N_list;
                adj_list=R_list;
            case 2
                if_connect=opts.R_connect;
                match_list=R_list;
                adj_list=N_list;
            case 3
                if_connect=opts.NR_connect;
                match_list=N_list.*R_list;
                adj_list=R_list./N_list;
        end
        style_line='k';
        if if_connect
            match_unique=unique(match_list);
            for ic=1:length(match_unique)
                pts_select=find(match_list==match_unique(ic));
                if length(pts_select)>1
                    adjvals=adj_list(pts_select);
                    [adjvals_sorted,adjvals_indices]=sort(adjvals);
                    %connect sequential values indexed by pts_select(adjvals_indices)
                    for cline=1:length(pts_select)-1
                        inds=pts_select(adjvals_indices(cline+[0 1]));
                        hp=plot3(log(N_list(inds))/log(2),log(R_list(inds))/log(2),vals(inds,iparam,1),style_line);
                        set(hp,'Color',colorspec);
                        set(hp,'LineWidth',opts.linewidth);
                    end
                end %enough points?
            end %ic: each match
        end %if_connect
    end %ictype
end %next param
%
XLim_noeb=round(log(opts.N_range)/log(2));
set(gca,'XLim',XLim_noeb+[-1 1]*max(abs(opts.eb_shape(1,:)))*opts.eb_size);
set(gca,'XTick',XLim_noeb(1):XLim_noeb(2));
set(gca,'XTickLabel',2.^get(gca,'XTick'));
xlabel(opts.N_label,'FontSize',opts.FontSize);
%
YLim_noeb=round(log(opts.R_range)/log(2));
set(gca,'YLim',YLim_noeb+[-1 1]*max(abs(opts.eb_shape(2,:)))*opts.eb_size);
set(gca,'YTick',YLim_noeb(1):YLim_noeb(2));
set(gca,'YTickLabel',2.^get(gca,'YTick'));
ylabel(opts.R_label,'FontSize',opts.FontSize);
%
if ~isempty(opts.zplane_draw)
    xlims=get(gca,'XLim');
    ylims=get(gca,'YLim');
    for iz=1:length(opts.zplane_draw)
        hp=plot3(xlims([1 1 2 2 1]),ylims([1 2 2 1 1]),repmat(opts.zplane_draw(iz),1,5),'k');
        set(hp,'LineWidth',opts.linewidth_zplane);
    end
end
%
%legend
if ~isempty(opts.labels) & ~isempty(leg_labels)
    if ~iscell(leg_labels)
        leg_labels=cellstr(leg_labels);
    end
    legend(leg_handles,leg_labels,'Location',opts.legend_location,'FontSize',opts.FontSize);
end
%
set(gca,'ZLim',opts.v_range);
set(gca,'ZScale',opts.v_scale);
zlabel(opts.v_label,'FontSize',opts.FontSize);
%
%set matlab standard fields
for ifield=1:length(std_fields)
    fn=std_fields{ifield};
    if isfield(opts_std.(fn),'MLname')
        MLname=opts_std.(fn).MLname;
    else
        MLname=fn;
    end
    set(gca,MLname,opts.(fn));
end
axis vis3d
%
opts_used=opts;
return

