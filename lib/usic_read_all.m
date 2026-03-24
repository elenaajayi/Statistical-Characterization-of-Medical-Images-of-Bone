function [patches,subpatches,meta,preproc_used]=usic_read_all(usic_path,usic_file,preproc)
% [patches,subpatches,meta,preproc_used]=usic_read_all(usic_path,usic_file,preproc) reads all of the 
% carotid ultrasound images as collected by Amanda Simons
%
%  Note:  this reads the patches AND preprocesses them into subpatches by choosing a fixed window and downsampling.
%   Also(unless keepdup=1) removes apparent duplicates.
%
% usic_path: path to image file.  If current folder, specify as .\. If empty, defaults to ...\Documents\JV\teaching\ASusi
% usic_file: file name of image file.  If empty, defaults to 'US_struct_nocursor.mat' 
%     if patient name D data cannot be found, then patient ID data is read from
%     a file with _Patient_Name appended to file name 
% preproc: preprocessing, [] if not supplied, modeled after ffdm_read_subpatch
%    preproc.XYLen: length (x,y) of patch to sample, if not supplied, defaults to [426 460]
%    preproc.XYOff: offset into patch, defaults to [75 211] if not supplied
%      above defaults are chosen to recapitulate Simons' values, which remove the text labelling:
%      row_length=size(US_table_nocursor.Image{1}(75:500,211:670),1);
%      col_length=size(US_table_nocursor.Image{1}(75:500,211:670),2);
%    prepoc.downsample: downsampling factor, defaults to 1, must divide XYLen
%    preproc.scale: full grayscale range, typically 0 to 255 for 8-bit uint
%
% Here, X is first coord, Y is second, considered as viewing a matrix, which is imagesc conventions:
%    X is 1 at top and increases down
%    Y is 1 at left and increases right
%
% patches: cell array containing patch images, without preprocessing
% subpatches:  cell array of preprocessed patches
% meta: patch metadata structure, fields modeled after those returned by ffdm_readmetadata
%   usic_path: path to data file
%   usic_file: file name read
%   npatches: number of patches
%   PID: patient ID
%   PID_num: patient ID number
%   PID_orig: original patient ID, may have a suffix (1). Images with the same PID but different PID_num's are identical
%   View_num: ones(npatches,1)
%   View: view , cell array of length npatches (always 'std')
%   view_names: {'std'}
%   XYPos: position of patch in original image (always (1,1))
%   XYLen: length of original patch
%   patch_count_table: number of patches for each patient and view, array of size max(PID_orig_num), length(view_name)
%   Examp_num: which example for this  patient (size is npatches, 1)
% preproc_used: preprocessing used
%
%  See also: FILLDEFAULT, FFDM_READMETADATA, FFDM_READ_SUBPATCH, USIC_READ_ALL_DEMO, MRIX_READ_ALL.
%
usic_path_default='\Documents\JV\teaching\ASusi\';
usic_file_default='US_struct_nocursor.mat';
usic_XYLen_def=[426 460];
usic_XYOff_def=[75 211];
if (nargin<1)
    usic_path=[];
end
if (nargin<2)
    usic_file=[];
end
if (nargin<3)
    preproc=[];
end
patches=cell(0);
subpatches=cell(0);
meta=struct;
preproc_used=struct;
%
%file locations for stimulus database
if isempty(usic_path)
    [dos_err,user_profile]=dos('echo %USERPROFILE%');
    user_profile=deblank(user_profile);
    usic_path=cat(2,user_profile,usic_path_default);
end
if isempty(usic_file)
    usic_file=usic_file_default;
end
meta.usic_path=usic_path;
meta.usic_file=usic_file;
%
%read the data file and Patient Name file if needed
%
fullname=cat(2,usic_path,filesep,usic_file);
if ~exist(fullname,'file')
    error(sprintf('cannot find %s',fullname));
end
s=getfield(load(fullname),'US_struct_nocursor');
disp(sprintf('Ultrasound Image data read from  %s',strrep(fullname,'\','/')));
Patient_Name=[];
if isfield(s,'Patient_Name')
    Patient_Name=getfield(s,'Patient_Name');
    disp(sprintf('Ultrasound Patient ID data read from  %s',strrep(fullname,'\','/')));
end
if isempty(Patient_Name)
    disp(sprintf('Ultrasound Patient ID data not found in %s',strrep(fullname,'\','/')));
    %form a file name by appending _Patient_Name to filename, whether or not .mat was specified
    fullname_PN=strrep(strrep(cat(2,fullname,'.mat'),'.mat.mat','.mat'),'.mat','_Patient_Name.mat');
    if ~exist(fullname_PN,'file')
        warning(sprintf('cannot find %s',fullname_PN));
        return
    end
    s_PN=load(fullname_PN);
    if isfield(s_PN,'Patient_Name')
        Patient_Name=getfield(s_PN,'Patient_Name');
        disp(sprintf('Ultrasound Patient ID data read from  %s',strrep(fullname_PN,'\','/')));
    end
    if isempty(Patient_Name)
        disp(sprintf('Ultrasound Patient ID data not found in %s',strrep(fullname_PN,'\','/')));
    end
end
%
npatches=length(s.Image);
meta.npatches=npatches;
meta.PID_orig=Patient_Name;
meta.PID=strrep(Patient_Name,'(1)','');
%
[PID_unique,ix,jx]=unique(meta.PID);
meta.PID_num=jx;
meta.XYLen=s.Image_Dimensions;
meta.XYPos=ones(npatches,2);
meta.views={'std'};
meta.View_num=ones(npatches,1);
meta.View=cell(npatches,1);
meta.patch_count_table=zeros(length(PID_unique),length(meta.views));
meta.Examp_num=zeros(npatches,1);
for ipatch=1:npatches
    meta.View{ipatch}=meta.views{1};
    meta.patch_count_table(meta.PID_num(ipatch),meta.View_num(ipatch))=1+meta.patch_count_table(meta.PID_num(ipatch),meta.View_num(ipatch));
    meta.Examp_num(ipatch)=sum(meta.PID_num(1:ipatch)==meta.PID_num(ipatch));
end
% 
preproc=filldefault(preproc,'scale',[0 255]);
preproc=filldefault(preproc,'XYOff',usic_XYOff_def); %could also be [0 0]
preproc=filldefault(preproc,'XYLen',usic_XYLen_def); % could also be min(meta.XYLen,[],1);
preproc=filldefault(preproc,'downsample',1);
preproc_used=preproc;
%
XY=preproc.XYLen;
d=preproc.downsample;
%
if any(mod(XY,d)>0)
    error(sprintf('requested downsampling (%3.0f %3.0f region size, downsampling by %3.0f) requires non-integer number of pixels per patch',XY,d));
end
if any(preproc.XYOff+XY>min(meta.XYLen,[],1)) | any(preproc.XYOff<0)
    disp(preproc)
    error(sprintf('region requested from patch exceeds patch boundaries, smallest patch is of size %3.0f %3.0f',min(meta.XYLen,[],1)));
end
%
for ipatch=1:npatches
    patches{ipatch,1}=double(s.Image{ipatch});
    region=patches{ipatch}(preproc.XYOff(1)+[1:XY(1)],preproc.XYOff(2)+[1:XY(2)]); %cut out requested region
    region_scaled=(double(region)-preproc.scale(1))/diff(preproc.scale); %rescale
    %downsample by block averaging
    region_scaled=reshape(region_scaled,[d,XY(1)/d,d,XY(2)/d]);
    subpatches{ipatch,1}=reshape(sum(sum(region_scaled,3),1)/d/d,XY/d);    
end
return
