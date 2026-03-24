%usic_read_all_demo:  demonstrate usic_read_all and check for duplicate images
%
%  See also: USIC_READ_ALL, MRIX_READ_ALL_DEMO.
%
if ~exist('usic_path') usic_path=[]; end
if ~exist('usic_file') usic_file=[]; end
if ~exist('preproc') preproc=struct(); end
if ~exist('colormap_name') colormap_name='gray'; end
if ~exist('npatches_to_show') npatches_to_show=6; end
[patches,subpatches,meta,preproc_used]=usic_read_all(usic_path,usic_file,preproc);
disp('Ultrasound data read.');
meta
preproc_used
preproc_used_string=sprintf('scale: [%7.2f %7.2f] XYOff: [%4.0f %4.0f] XYLen [%4.0f %4.0f] downsample %1.0f',...
    preproc_used.scale,preproc_used.XYOff,preproc_used.XYLen, preproc_used.downsample);
disp(preproc_used_string);
maxmax=0;
for isubj=1:max(unique(meta.PID_num))
    subpatch_nos=find(meta.PID_num==isubj)';
    matchstring=sprintf(' first patch is %4.0f',subpatch_nos(1));
    if (length(subpatch_nos))>1
    	for isp=2:length(subpatch_nos)
            maxdiff=max(abs(subpatches{subpatch_nos(isp)}(:)-subpatches{subpatch_nos(1)}(:)));
            matchstring=cat(2,matchstring,sprintf(' c/w %4.0f: diff %5.3f',subpatch_nos(isp),maxdiff));
        end
        maxmax=max(maxmax,maxdiff);
    end
    disp(sprintf('PID %2.0f: %s, number of subpatches: %3.0f (check: %3.0f),%s',...
        isubj,meta.PID{subpatch_nos(1)},length(subpatch_nos),meta.patch_count_table(isubj,1),matchstring));
end
%show patches
zrange_patch=[+Inf -Inf];
for ipatch=1:meta.npatches
    zrange_patch(1)=min(zrange_patch(1),min(patches{ipatch}(:)));
    zrange_patch(2)=max(zrange_patch(2),max(patches{ipatch}(:)));
end
if npatches_to_show>0
    if (npatches_to_show==1)
        patch_list=1;
    else
        patch_list=round(1+(meta.npatches-1)*[0:npatches_to_show-1]/(npatches_to_show-1));
    end
    [nr,nc]=nicesubp(npatches_to_show,0.7);
    tstring=meta.usic_file;
    for ipsp=1:2 %patch or subpatch
        if (ipsp==1)
            tstring2='patch';
            zrange=zrange_patch;
        else
            tstring2='subpatch';
            zrange=[0 1];
        end
        figure;
        set(gcf,'Position',[50 100 1200 750]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',cat(2,tstring2,'es from ',tstring));
        for ipatch_to_show=1:npatches_to_show
            ipatch=patch_list(ipatch_to_show);
            subplot(nr,nc,ipatch_to_show);
            if (ipsp==1)
                imagesc(patches{ipatch},zrange);
            else
                imagesc(subpatches{ipatch},zrange);
            end
            colormap(colormap_name);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            xlabel(sprintf('%s %2.0f PID %2.0f (%s) ex %1.0f',tstring2,ipatch,meta.PID_num(ipatch),meta.PID{ipatch},meta.Examp_num(ipatch)),'FontSize',8);
            axis equal;
            axis tight;
        end
    end
    axes('Position',[0.02,0.02,0.01,0.01]); %for text
    text(0,0,sprintf('%ses from %s: range [%7.2f %7.2f]',tstring2,tstring,zrange),'Interpreter','none');
    axis off;
    if (ipsp==2)
        axes('Position',[0.02,0.06,0.01,0.01]); %for text
        text(0,0,cat(2,'preprocessing used: ',preproc_used_string));
        axis off;
    end
end
