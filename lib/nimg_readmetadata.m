function meta=nimg_readmetadata(nimg_meta_path,nimg_meta_file,NR_orig)
% meta=nimg_readmetadata(nimg_meta_path,nimg_meta_file,NR_orig) reads the metadata for a collection of
%  natural images from Hermundstad et al., eLife.
%
%  this includes pointers to patches that are in focus, where focus is determined at multiple different scales as defined by NR_orig
%
%  see C:\Users\jdvicto\Dropbox\Textures\RawNIstats\ImagePatchStats\code\extractingPatches\README.rtf
%
%  individual images play the role of "patients", but also, grouping of images into galleries is available
%
%  nimg_meta_path: path to metadata file, if empty, then 'C:\Users\jdvicto\Dropbox\Textures\RawNIstats\ImagePatchStats\data\'
%  nimg_meta_file: metadata filename, if empty, then 'PPandNIstats_allAnalyses.mat'
%  NR_orig: [Norig Rorig], the downsampling and region size value used to determine which patches are in focus.  Both must match.
%     If empty or omitted, then a pair is requested
%
% meta: patch metadata structure, fields modeled after those returned by ffdm_readmetadata
%   nimg_meta_path: metadata path
%   nimg_meta_file: metadata file name
%   npatches: number of patches
%   NR_orig: downsampling and patch size (in downsampled checks) used for selecting in-focus patches
%   nsubjs: number of subjects (patients); here interpreted as number of galleries
%   PID: patient (subject) ID
%   PID_num: patient ID number, ASCII equiv of PID
%   image_file: original image file for each patient, typically contains a directory (gallery), e.g., 'cd01A/DSC_0002_LUM.mat'
%   View_num: ones(npatches,1)
%   View: view , cell array of length npatches (always 'std')
%   views: {'std'}
%   XYPos: position of patch in original image, in original pixels.  
%      XYPos reverses the order in dataNI.indA(iNR_select).focus.ic.[x y], so that
%      XYPos(:,1) refers to focus.ic.y, i.e., the *first* coordinate of the LUM.mat files, which is top to bottom.
%      XYPos(:,2) refers to focus.ic.x, i.e., the *second* coordinate of the LUM.mat files, which is top to bottom.
%   XYLen: length of original patch, this is Norig*Rorig
%   patch_count_table: number of patches for each patient and view, array of size max(PID_orig_num), length(view_name)
%   Examp_num: which example for this  patient (size is npatches, 1)
%      fields related to how the images are organized into galleries (specific to Penn database, not related to fields in ffdm_readmetadata)
%      These fields are also indepenent of the choice of N and R
%   gal_names_sorted: cell(ngals,1), gallery (diretory) names, in original Penn index order
%   gal_img_list:  cell(n_gals,1), list of image numbers in each gallery
%   gal_num: size(n_imgs,1) gallery number for each image
%
%  26Oct21: minor typos fixed
%
%  See also: FILLDEFAULT, USIC_READ_ALL, MRIX_READ_ALL, FFDM_BTC_CALC_GEN, NIPATCHES_GETGALS, NIMG_READ_PATCHES, NIMG_SELECT_PATCHES.
%
meta=struct;
if isempty(nimg_meta_path) nimg_meta_path='C:\Users\jdvicto\Dropbox\Textures\RawNIstats\ImagePatchStats\data\'; end
if isempty(nimg_meta_file) nimg_meta_file='PPandNIstats_allAnalyses.mat'; end
if (nargin<=2)
    NR_orig=[];
end
if length(NR_orig)~=2
    NR_orig=[];
end
if isempty(NR_orig)
    NR_orig=[0 0]; %so match won't fail
end
meta.nimg_meta_path=nimg_meta_path;
meta.nimg_meta_file=nimg_meta_file;
dataNI=getfield(load(cat(2,nimg_meta_path,nimg_meta_file)),'dataNI');
n_NR=length(dataNI.indA);
disp(sprintf('image patch information loaded with %2.0f (N,R) sets',n_NR));
%
N_list=zeros(n_NR,1);
R_list=zeros(n_NR,1);
NR_list=zeros(n_NR,1);
for iNR=1:n_NR
    indA=dataNI.indA(iNR);
    N_list(iNR)=indA.N;
    R_list(iNR)=indA.R;
    NR_list(iNR)=N_list(iNR)*R_list(iNR);
end
iNR_select=find(all(repmat(NR_orig,n_NR,1)==[N_list R_list],2));
if length(iNR_select)~=1
    disp('(N,R) value pair either not recognized or not supplied.');
    for iNR=1:n_NR
        indA=dataNI.indA(iNR);
        ic_all=indA.ic;
        ic_focus=indA.focus.ic;
        disp(sprintf(' set %3.0f: N=%3.0f, R=%3.0f, NR=%4.0f, total patches: %8.0f (tot area: %15.0f), in-focus patches: %8.0f (tot area: %15.0f)',...
            iNR,N_list(iNR),R_list(iNR),NR_list(iNR),length(ic_all.image),length(ic_all.image)*NR_list(iNR)^2,...
            length(ic_focus.image),length(ic_focus.image)*NR_list(iNR)^2));
    end
    iNR_select=getinp('choice','d',[1 n_NR]);
end
indA=dataNI.indA(iNR_select);
ic_focus=dataNI.indA(iNR_select).focus.ic;
npatches=length(ic_focus.image);
disp(sprintf(' This is %7.0f in-focus patches, all of size %3.0f x %3.0f pixels, considered as %3.0f x %3.0f checks of %3.0f x %3.0f pixels each',...
    npatches,indA.N*indA.R,indA.N*indA.R,indA.R,indA.R,indA.N,indA.N));
NR=indA.N*indA.R;
meta.npatches=npatches;
meta.NR_orig=[indA.N,indA.R];
meta.nsubjs=length(ic_focus.name);
meta.PID=cellstr(num2str(ic_focus.image'));
meta.PID_num=ic_focus.image';
meta.image_file=cell(npatches,1);
meta.View_num=ones(npatches,1);
meta.View=cellstr(repmat('std',npatches,1));
meta.views={'std'};
%   exchange the coordinates to take into account AH conventions  that .x refers to second (long) coordinate
%   also position is in units of original pixels but ic_focus.[x,y] are in units of NR
meta.XYPos=1+([ic_focus.y' ic_focus.x']-1)*NR;
%
meta.XYLen=repmat(NR,npatches,2);
meta.patch_count_table=zeros(length(ic_focus.name),1);
meta.Examp_num=zeros(npatches,1);
img_no_prev=0;
for ipatch=1:npatches
    img_no=ic_focus.image(ipatch);
    %determine which patch within the image, by seeing if image number has changed
    if (img_no==img_no_prev)
        examp_num=examp_num+1;
    else
        examp_num=1;
    end
    meta.Examp_num(ipatch)=examp_num;
    img_no_prev=img_no;
    meta.image_file{ipatch}=ic_focus.name{img_no};
    meta.patch_count_table(img_no)=meta.patch_count_table(img_no)+1;
end
%fields related to organization of files into gallery
%determine a list of unique gallery names, in the original order 
file_names=ic_focus.name;
n_imgs=length(file_names);
gal_names=cell(n_imgs,1);
for img_index=1:n_imgs
    im_fullname=deblank(file_names{img_index});
    slash_pos=min(strfind(im_fullname,'/'));
    gal_names{img_index}=im_fullname(1:slash_pos-1);
end
[gal_names_unique,gal_first]=unique(gal_names);
[gal_first_sorted,gal_first_index]=sort(gal_first); %want the galleries in original order
n_gals=length(gal_names_unique);
meta.gal_names_sorted=cell(n_gals,1);
meta.gal_img_list=cell(n_gals,1);
meta.Gal_num=zeros(n_imgs,1);
for igal=1:n_gals
    meta.gal_names_sorted{igal}=gal_names_unique{gal_first_index(igal)};
    meta.gal_img_list{igal}=strmatch(meta.gal_names_sorted{igal},gal_names);
    meta.Gal_num(meta.gal_img_list{igal})=igal;
end
return
