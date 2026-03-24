%bone_dbase_demo:  demo of reading the bone database original images and rois
%
%   See also:  BONE_READ_XLS, BONE_PSPEC_DEMO.
%
file_name=getinp('bone database file name','s',[0 1],'bonedatabase.xlsx');
[im_data_all,roi_data_all,roi_sizes,roi_origimg,im_files,opts_used]=bone_read_xls(file_name,setfield([],'if_log',1));
n_entries=length(im_data_all);
%
show_choice=-1;
while (show_choice~=0)
    show_choice=getinp('entry number to show (0 to quit)','d',[0 n_entries]);
    if (show_choice>0)
        im_file=im_files{show_choice};
        im_data=im_data_all{show_choice};
        roi_list=find(roi_origimg==show_choice); %list of rois that are in this image
        n_roi=length(roi_list);
        im_range=[min(im_data(:)) max(im_data(:))]; %image range
        disp(sprintf(' image values range from %7.3f to %7.3f',im_range));
        figure;
        set(gcf,'NumberTitle','off');
        set(gcf,'Name',im_file);
        set(gcf,'Position',[100 100 1200 800]); %position and size of figure window: [lowX lowY widthX widthY]
        [n_rows,n_cols]=nicesubp(1+n_roi,.7); %nice arrangement of subplots
        subplot(n_rows,n_cols,1); %where to put next plot
        imagesc(im_data,[0 im_range(2)]);
        axis equal; %equal scales for pixels horiz and vert
        axis tight;
        title(im_file);
        colormap gray;
        %retrieve the ROIs
        for i_roi=1:n_roi
            roi_suffix=zpad(i_roi,2); %pad the roi number to two digits
            roi_data=roi_data_all{roi_list(i_roi)};
            disp(sprintf(' roi %2.0f:  (height,width) is %5.0f x %5.0f',i_roi,size(roi_data)));
            subplot(n_rows,n_cols,1+i_roi); %next subplot, first subplot reserved for full image
            imagesc(roi_data,[0 im_range(2)]);
            axis equal; %equal scales for pixels horiz and vert
            axis tight;
            title(sprintf('%s: roi %2.0f',im_file,i_roi));
            colormap gray;
        end %each ROI
    end
end
disp('done.');
