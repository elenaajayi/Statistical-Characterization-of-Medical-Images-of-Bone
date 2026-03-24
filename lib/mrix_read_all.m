function [subpatches,meta,lims]=mrix_read_all(mrix_path,mrix_subdir,mrix_file)
% [subpatches,meta,lims]=mrix_read_all(mrix_path,mrix_subdir,mrix_file) reads MRI patches as collected by Yueyang Xu
%   and used in Frontiers (Xu et al., 2019)
%
% orients images so that top of head is on top, front of head is on right
%  for ADNI and OASIS, this means flipping on Y
%  for TNS databases (HV and Patient), this means tranposing and then flipping on Y
% reorientation is determined by string supplied in mrix_subdir
%
%  Note:  MRI data have already been preprocessed into ROIs of identical size, 64 x 64, and these are
%   the subpatches -- so there is no option to retrieve full patches, or for alternative preprocessing.
%
%  mrix_path: path to MRI files, if empty, defaults to ...\Dropbox\From_SamXu\data\
%      if ONLY mrix_path is supplied, then other variables are ignored and
%      subpatches is returned as a cell array of possible  values for mrix_subdir
%      If supplied (or default) path is not found, also looks for
%      ...\Documents\Dropbox\...
%  mrix_subdir: subdirectory in mrix path, typically 'ADNI Statistics','HV Statistics', Oasis Statistics', or 'Patient Statistics'
%  mrix_file; file name contining data, typically 'sampleList.mat', defaults to 'sampleList.mat'.  Data expected in 'sampleList'
%
% subpatches:  cell array of subpatches.  Note that MRI data have already been processed by Xu prior to
%   creation of sampleList.  So there is no option to retrieve full patches, or for alternative preprocessing
% meta: (sub)patch metadata structure, fields modeled after those returned by ffdm_readmetadata
%   mrix_path: path to data file
%   mrix_subdir: subdirectory containing data file
%   mrix_file: data file
%   npatches: number of patches, equals number of subpatches
%   nsubjs: number of subjects (patients); this is 1 for simulations as there is only one brain
%   PID: patient (subject) ID
%   PID_num: patient ID number
%   image_file: original image file for each patient
%   View_num: ones(npatches,1)
%   View: view , cell array of length npatches (always 'std')
%   view_names: {'std'}
%   XYPos: position of patch in original image (always (1,1))
%   XYLen: length of original patch
%   patch_count_table: number of patches for each patient and view, array of size max(PID_orig_num), length(view_name)
%   Examp_num: which example for this  patient (size is npatches, 1)
% lims: high and low values across all subpatches
%
%  See also: FILLDEFAULT, FFDM_READMETADATA, FFDM_READ_SUBPATCH, USIC_READ_ALL, MRIX_READ_ALL_DEMO, FFDM_BTC_CALC_GEN, NIMG_READMETADATA.
%
mrix_path_default='\Dropbox\From_SamXu\data\';
dropbox_path_alts={'\Dropbox','Documents\Dropbox'}; %alternate paths to dropbox from %USERPROFILE%
mrix_subdir_default='Statistics'; %but a non-default value should be supplied, such as 'ADNI Statistics'
mrix_file_default='sampleList.mat';
if (nargin<1)
    mrix_path=[];
end
if (nargin<2)
    mrix_subdir=[];
end
if (nargin<3)
    mrix_file=[];
end
subpatches=cell(0);
lims=[];
meta=struct;
%file locations for stimulus database
if isempty(mrix_path)
    [dos_err,user_profile]=dos('echo %USERPROFILE%');
    user_profile=deblank(user_profile);
    mrix_path=cat(2,user_profile,mrix_path_default);
end
if ~exist(mrix_path,'file')
    ifound=0;
    ialt=0;
    while (ialt<length(dropbox_path_alts)) & ifound==0
        ialt=ialt+1;
        mrix_path_try=strrep(mrix_path,'Dropbox',dropbox_path_alts{ialt});
        ifound=exist(mrix_path_try,'file');
        if (ifound==0)
            disp(sprintf('cannot find path %s',mrix_path_try))
        else
            disp(sprintf('found path       %s',mrix_path_try))
            mrix_path=mrix_path_try;
        end
    end
    if ~exist(mrix_path,'file')
        error(sprintf('cannot find %s or alternates',mrix_path));
    end
end
meta.mrix_path=mrix_path;
ascii_linesep=10;
if (nargin<=1) %special behavior for only one input
    subdirs=cell(0);
    [doserr,dosout]=dos(cat(2,'dir ',mrix_path,' /b'));
    while ~isempty(dosout)
        [subdir,dosout]=strtok(dosout,char(ascii_linesep));
        if ~isempty(deblank(subdir))
            subdirs{end+1}=deblank(subdir);
        end
    end
    subpatches=subdirs;
    return
end
meta.mrix_subdir=mrix_subdir;
if isempty(mrix_file) mrix_file=mrix_file_default; end
meta.mrix_file=mrix_file;
%
%read the data file
%
fullname=cat(2,mrix_path,filesep,mrix_subdir,filesep,mrix_file);
if ~exist(fullname,'file')
    error(sprintf('cannot find %s',fullname));
end
s=getfield(load(fullname),'sampleList');
disp(sprintf('MRI data read from  %s',strrep(fullname,'\','/')));
nsubjs=length(s);
meta.views={'std'};
meta.patch_count_table=zeros(nsubjs,length(meta.views));
%first determine total number of slices
for isubj=1:nsubjs
    meta.patch_count_table(isubj,1)=size(s{isubj}.imagepixel,3);
end
npatches=sum(meta.patch_count_table(:));
lims=[+Inf,-Inf];
meta.npatches=npatches;
meta.XYLen=zeros(npatches,2);
meta.XYPos=ones(npatches,2);
meta.View_num=ones(npatches,1);
meta.View=cell(npatches,1);
meta.Examp_num=zeros(npatches,1);
meta.PID_num=zeros(npatches,1);
meta.PID=cell(npatches,1);
meta.image_file=cell(nsubjs,1);
meta.convRegion=zeros(npatches,4);
meta.slicelist=zeros(npatches,1);
subpatches=cell(npatches,1);
index=0;
for isubj=1:nsubjs
    nslices=meta.patch_count_table(isubj,1);
    meta.XYLen(index+[1:nslices],:)=repmat(size(s{isubj}.imagepixel(:,:,1)),nslices,1);
    meta.convRegion(index+[1:nslices],:)=s{isubj}.convRegion;
    if isfield(s{isubj},'slicelist')
        meta.slicelist(index+[1:nslices])=s{isubj}.slicelist;
    else
        meta.slicelist(index+[1:nslices])=[1:nslices]';
    end
    meta.image_file{isubj}=strrep(s{isubj}.filename,'F:\SharedFolder\output','');
    meta.PID_num(index+[1:nslices])=isubj;
    meta.Examp_num(index+[1:nslices])=[1:nslices];
    %retrieve slices
    PIDstring=cat(2,'P',zpad(isubj,4));
    lims(1)=min(lims(1),min(s{isubj}.imagepixel(:)));
    lims(2)=max(lims(2),max(s{isubj}.imagepixel(:)));
    %orient images so that in imagesc display, top of head is at top and front of head is at right
    imagepixel_raw=double(s{isubj}.imagepixel);
    switch lower(mrix_subdir)
        case {'adni statistics','oasis statistics','adni','oasis'}
            images_oriented=flipdim(imagepixel_raw,2);
            if (isubj==1)
                disp('images flipped on dimension 2');
            end
        case {'hv statistics','patient statistics','tns'}
            images_oriented=flipdim(permute(imagepixel_raw,[2 1 3]),2);
            if (isubj==1)
                disp('images transposed and flipped on dimension 2');
            end
        otherwise        
            images_oriented=imagepixel_raw;
            if (isubj==1)
                disp('no transormation appied to images');
            end
    end
    for islice=1:nslices
        subpatches{index+islice}=images_oriented(:,:,islice);
        meta.PID{index+islice}=PIDstring;
        meta.View{index+islice}=meta.views{1};
    end
    index=index+nslices;
end
disp(sprintf('found %5.0f subjects, %5.0f slices',nsubjs,npatches));
return
