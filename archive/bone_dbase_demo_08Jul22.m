%bone_dbase_demo:  demo of reading the bone database and its images
%
file_name=getinp('bone database file name','s',[0 1],'bonedatabase.xlsx');
[numeric,text,raw]=xlsread(file_name);
n_entries=size(text,1)-1;
disp(sprintf('database has %3.0f entries and a header row',n_entries));
headers=text(1,:);
disp(headers);
%search in headers to get proper columns, in case columns are added at a later time
%
source_col=strmatch('Site',headers); %column for source/site
width_col=strmatch('Width',headers); %column for Width
height_col=strmatch('Height',headers); %column for Height
nroi_col=strmatch('how many ROI?',headers); %column for number of ROIs
orig_col=strmatch('Original Image',headers); %column for image name
type_col=strmatch('Bone type',headers); %column for bone type
show_choice=-1;
while (show_choice~=0)
    show_choice=getinp('entry number to show (0 to quit)','d',[0 n_entries]);
    if (show_choice>0)
        im_file=raw{1+show_choice,orig_col};
        im_size_db=[raw{1+show_choice,height_col} raw{1+show_choice,width_col}]; %size of image (height, width) according to database
        n_roi=raw{1+show_choice,nroi_col};
        disp(sprintf('Requested entry is from %10s, image %15s, size is %5.0f(h) x %5.0f(w); %3.0f ROIs, type: %s',...
            raw{1+show_choice,source_col},im_file,im_size_db,n_roi,raw{1+show_choice,type_col}));
        % read the image
        im_data=imread(im_file);
        if size(im_data,3)==3 %is the image RGB?
            im_data=rgb2gray(im_data); %if so, convert to gray
        end
        disp(sprintf(' full image (height, width) is %5.0f x %5.0f',size(im_data)));
        if all(size(im_data)==im_size_db)
            disp('size agrees with database entry');
        else
            disp('size *disagrees* with database entry');
        end
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
        %read the ROIs
        im_file_base=im_file(1:find(im_file=='.')-1); %remove extension from im_file
        im_file_ext=im_file(find(im_file=='.'):end); %just im_file extension
        for i_roi=1:n_roi
            roi_suffix=zpad(i_roi,2); %pad the roi number to two digits
            roi_file=cat(2,im_file_base,'_',roi_suffix,im_file_ext);
            roi_data=imread(roi_file);
            if size(roi_data,3)==3
                roi_data=rgb2gray(roi_data);
            end
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

