%mrix_read_all_demo:  demonstrate mrix_read_all 
%
%  See also: MRIX_READ_ALL, USIC_READ_ALL_DEMO.
%
if ~exist('mrix_path') mrix_path=[]; end
if ~exist('mrix_file') mrix_file=[]; end
if ~exist('colormap_name') colormap_name='gray'; end
if ~exist('nsubjs_to_show') nsubjs_to_show=6; end
if ~exist('nslices_to_show') nslices_to_show=8; end
%
%call in special mode to get subdirectory options
subdir_choices=mrix_read_all(mrix_path);
for ichoice=1:length(subdir_choices)
    disp(sprintf(' %2.0f->%s',ichoice,subdir_choices{ichoice}));
end
ds_list=getinp('choice(s)','d',[1 length(subdir_choices)],[1:length(subdir_choices)]);
subpatches=cell(1,length(ds_list));
meta=cell(1,length(ds_list));
lims=cell(1,length(ds_list));
for ds_ptr=1:length(ds_list)
    ids=ds_list(ds_ptr);;
    mrix_subdir=subdir_choices{ids};
    disp(sprintf(' processing %s',mrix_subdir));
    [subpatches{ids},meta{ids},lims{ids}]=mrix_read_all(mrix_path,mrix_subdir,mrix_file);
    disp(sprintf('found %5.0f subjects, %5.0f slices, limits: [%8.3f %8.3f]',size(meta{ids}.patch_count_table,1),meta{ids}.npatches,lims{ids}));
    nsubjs=size(meta{ids}.patch_count_table,1);
    if nsubjs_to_show>0
        if (min(nsubjs,nsubjs_to_show)==1)
            subj_list=1;
        else
            subj_list=round(1+(nsubjs-1)*[0:min(nsubjs,nsubjs_to_show)-1]/(min(nsubjs,nsubjs_to_show)-1));
        end
        nr=length(subj_list);
        nc=nslices_to_show;
        tstring=meta{ids}.mrix_subdir;
        zrange=lims{ids};
        figure;
        set(gcf,'Position',[50 100 1200 750]);
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',tstring);
        for isubj_ptr=1:length(subj_list)
            isubj=subj_list(isubj_ptr);
            nslices=meta{ids}.patch_count_table(isubj);
            if (nslices_to_show==1)
                slice_list=round((1+nslices)/2); %middle slice
            else
                slice_list=round(1+(nslices-1)*[0:nslices_to_show-1]/(nslices_to_show-1));
            end
            for islice_ptr=1:nslices_to_show
                islice=slice_list(islice_ptr);
                subplot(nr,nc,islice_ptr+(nslices_to_show)*(isubj_ptr-1));
                ipatch=intersect(find(meta{ids}.PID_num==isubj),find(meta{ids}.Examp_num==islice));
                if length(ipatch)~=1
                    warning(sprintf(' for subject %3.0f and slice %3.0f, %3.0f patches found',isubj,islice,length(ipatch)));
                else
                    imagesc(subpatches{ids}{ipatch},zrange);
                    colormap(colormap_name);
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
                    xlabel(sprintf('p%5.0f: s%3.0f ex%3.0f',ipatch,isubj,islice),'FontSize',8);
                    axis equal;
                    axis tight;
                end
            end %islice_ptr
        end %isubj_ptr
        %
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,sprintf('MRIs from %s: range [%7.2f %7.2f]',tstring,zrange),'Interpreter','none');
        axis off;
    end %npatches_to_show=0
end %ds_ptr
