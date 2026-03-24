function [im_data_all,roi_data_all,roi_sizes,roi_origimg,im_files,opts_used]=bone_read_xls(file_name,opts)
%[im_data_all,roi_data_all, roi_sizes,roi_origimg,im_files,opts_used]=bone_read_xls(file_name,opts)
%
% reads a database of images (typically bone radiographs or scintigraphy)
%
% file_name: spreadsheet file name
% opts: options
%   opts.if_log: 1 to log progress
%   opts.nroi_col_tag: string that marks the column with number of ROI's
%   opts.roi_infix: infix for roi file names,
%       file names will be [base file name]_[roi_infix]01.[ext], etc, where
%       ext is the extension of the original image file, typically jpg
%       roi_infix='' for original roi files
%       roi_infix='i' for roi files restricted to interior of image
%   opts.file_path: path to spreadsheet, empty if omitted
%   opts.im_path: path to image files, empty if omitted
%   opts.roi_path: path to roi files, empty if omittted
%
% im_data_all: cell array of original images; number of entries read from database is length(im_data_all)
% roi_data_all: cell array of rois, number of rois read is length(roi_data_all)
% roi_sizes: array of size (nimgs,2), size of original images
% roi_orig_img: original image number for each roi
% im_files: cell array of original image file names
% opts_used: options used
%
% 12Sep22: add options to change file name of roi files (opts.roi_infix) and columm for n_roi (opts.nroi_col_tag)
% 22Sep22: fix bug when number of roi entry is blank
% 24Sep22: add opts.file_path, opts.im_path,opts.roi_path
%
% See also:  FILLDEFAULT, BONE_DBASE_DEMO, BONE_PSPEC_DEMO, BONE_BTC_DEMO.
%
if (nargin<=1)
    opts=[];
end
opts=filldefault(opts,'if_log',0);
opts=filldefault(opts,'nroi_col_tag','how many ROI?');
opts=filldefault(opts,'roi_infix','');
opts=filldefault(opts,'file_path','');
opts=filldefault(opts,'im_path','');
opts=filldefault(opts,'roi_path','');
opts_used=opts;
[numeric,text_xls,raw]=xlsread(cat(2,opts.file_path,file_name));
n_entries=size(text_xls,1)-1;
headers=text_xls(1,:);
if opts.if_log
    disp(sprintf('database has %3.0f entries and a header row',n_entries));
    disp(headers);
end
%search in headers to get proper columns, in case columns are added at a later time
%
source_col=strmatch('Site',headers); %column for source/site
width_col=strmatch('Width',headers); %column for Width
height_col=strmatch('Height',headers); %column for Height
%nroi_col=strmatch('how many ROI?',headers); %column for number of ROIs
nroi_col=strmatch(opts.nroi_col_tag,headers);
orig_col=strmatch('Original Image',headers); %column for image name
type_col=strmatch('Bone type',headers); %column for bone type
im_data_all=cell(0);
roi_data_all=cell(0);
roi_sizes=zeros(0,2); %roi_sizes(:,1) is height, roi_sizes(:,2) is width
roi_origimg=zeros(0,1); %roi_origimg is original image number
n_rois=0;
%
for i_entry=1:n_entries %inspect each entry
    im_file=raw{1+i_entry,orig_col};
    im_size_db=[raw{1+i_entry,height_col} raw{1+i_entry,width_col}]; %size of image (height, width) according to database
    n_roi=raw{1+i_entry,nroi_col};
    if opts.if_log
        disp(sprintf('Requested entry is from %10s, image %15s, size is %5.0f(h) x %5.0f(w); %3.0f ROIs, type: %s',...
            raw{1+i_entry,source_col},im_file,im_size_db,n_roi,raw{1+i_entry,type_col}));
    end
    % read the image
    im_data=imread(cat(2,opts.im_path,im_file));
    im_files{i_entry}=im_file;
    if size(im_data,3)==3 %is the image RGB?
        im_data=rgb2gray(im_data); %if so, convert to gray
    end
    im_data_all{i_entry}=im_data;
    im_range=[min(im_data(:)) max(im_data(:))]; %image range
    if (opts.if_log)
        disp(sprintf(' full image (height, width) is %5.0f x %5.0f',size(im_data)));
        if all(size(im_data)==im_size_db)
            disp('size agrees with database entry');
        else
            disp('size *disagrees* with database entry');
        end
        disp(sprintf(' image values range from %7.3f to %7.3f',im_range));
    end
    if isnan(n_roi)
        disp('no rois listed');
    else
        % read the individual rois
        im_file_base=im_file(1:find(im_file=='.')-1); %remove extension from im_file
        im_file_ext=im_file(find(im_file=='.'):end); %just im_file extension
        for i_roi=1:n_roi
            roi_suffix=zpad(i_roi,2); %pad the roi number to two digits
            roi_file=cat(2,im_file_base,'_',opts.roi_infix,roi_suffix,im_file_ext);
            roi_data=imread(cat(2,opts.roi_path,roi_file));
            if size(roi_data,3)==3
                roi_data=rgb2gray(roi_data);
            end
            n_rois=n_rois+1; %one more roi
            roi_data_all{n_rois}=roi_data;
            roi_sizes(n_rois,:)=size(roi_data);
            roi_origimg(n_rois)=i_entry;
            if (opts.if_log)
                disp(sprintf(' roi %2.0f:  (height,width) is %5.0f x %5.0f',i_roi,size(roi_data)));
            end
        end %each ROI
    end
end
return
