%bone_btc_idemo:  compute local image statistics (binary texture coordinates) from database and its images
%
% modified from bone_pspec_demo
%
%  uses logic of ffdm_btc_calc_gen for subdivision of bigpatches into patches
%   and then subsamplng of patches by N, so that patches are of size N*R
%   Note that N*R must divide the size of bigpatch
%
%  Note also that the downsampling by N precedes the FFT, so the size of the unpadded region 
%    for FFT is R x R, where patch_size=N*R.
%  Note that "stepping" option (S) of ffdm_btc_calc_gen is not implemented.
%
%  03Aug22
%   * use subroutine for reading
%   * exclude a bigpatch if any subpatch is constant
%   * add heatmap of btc coords
%   * fix some prompts
%
%    See also:  FFDM_BTC_CALC_GEN, BTC_DEFINE, BONE_PSPEC_DEMO, GLIDER_MAPUBI.
%
jpeg_max=256; %maximum value in any jpeg file
%
btc_dict=btc_define; %create a structure that defines the binary texture coordinates
btc_checkdef=[0 0;0 1;1 0;1 1]; %define the offsets for a 2 x 2 block
btc_ng=2; %number of gray levels (just black and white)
btc_n=length(btc_dict.codel); %number of coords, typically 10
btc_reshape=repmat(btc_ng,1,size(btc_checkdef,1)); %typically, [2 2 2 2];
btc_nconfigs=btc_ng.^size(btc_checkdef,1); %number of block configurations, typically, 16=2^4 
%
file_name=getinp('bone database file name','s',[0 1],'bonedatabase.xlsx');
opts=struct;opts.if_log=1;opts.roi_infix='i';% This invokes bone_read.xlsx for sub orig img
opts.nroi_col_tag='how many iroi'; %String for the header roi that calls for iroi
[im_data_all,roi_data_all,roi_sizes,roi_origimg,im_files,opts_used]=bone_read_xls(file_name,opts);
n_entries=length(im_data_all);
n_rois=length(roi_data_all);
minsize=min(roi_sizes(:));
minsize_pwr2=2^floor(log2(minsize));
disp(sprintf('minimum roi size, minimum power-of-2: %5.0f %5.0f',minsize,minsize_pwr2));
%
%changes from here down w.r.t. bone_pspec_demo
%
%  patches in bone_pspec_demo become bigpatches here
%  bigpatches are cut into patches of eddge length N*R, 
%  where N is the downsampling of pixels R is the number of downsampled values in one edge of the patch
%
%ask for bigpatch size (suggest minsize_pwr2, but need not be)
%
ifok=0;
while (ifok==0)
    bigpatch_size=getinp('bigpatch size, should be power of 2','d',[32 1024],minsize_pwr2);
    ifok=(bigpatch_size==2.^floor(log2(bigpatch_size)));
end
%
% set up a range of N and R values, for which N and R are both powers of 2,
% N*R is a factor of bigpatch_size, and R is never less than 16.
%
if ~exist('R_min_log2') R_min_log2=4; end
R_min=2.^R_min_log2;
% from ffdm_btc_calc_gen
NR_max_poss_log2=floor(log(min(bigpatch_size))/log(2)); %was roi_minsize ffdm_btc_calc_gen
NR_max_poss=2.^NR_max_poss_log2;
N_list_log2_def=[1:NR_max_poss_log2-R_min_log2];
if_NRS_ok=0;
Sstep_char={' ','+'}; %+ will indicate that phases are stepped
while (if_NRS_ok==0)
    N_list_log2=getinp('values of blocking (log 2(N))','d',[0 max(N_list_log2_def)],[0:max(N_list_log2_def)]);
    N_list=2.^N_list_log2;
    NR_log2_list=getinp('values of patch size in pixels (log 2(NR))','d',[N_list_log2(1)+2 NR_max_poss_log2],repmat(NR_max_poss_log2,1,length(N_list)));
    NR_list=2.^NR_log2_list;
    NR_max_log2=NR_max_poss_log2; 
    NR_max=2.^NR_max_log2;
    R_list=NR_list./N_list;
    disp('analysis parameters so far:')
    for iNR=1:length(N_list)
        fprintf(' %2.0f-> blocking (N) =%4.0f    blocks per patch edge (R) =%4.0f    patch size (NR)=%4.0f\n',...
            iNR,N_list(iNR),R_list(iNR),N_list(iNR)*R_list(iNR));
    end
%     if_S_ok=0;
%     while if_S_ok==0
%         S_list_log2=getinp('scales for calculating image statstics (log2(S), must be less than correpsponding values of R/4)',...
%             'd',[0 round(log(max(R_list))/log(2))],zeros(1,length(N_list)));
%         S_list=2.^S_list_log2;
%         if_S_ok=all(S_list<=(R_list/4));
%     end
%     Sstep_list=getinp('1 to step the starting phase for scaling','d',[0 1],zeros(1,length(N_list)));
%
%     disp(sprintf('analysis parameters including scaling, maximum total patch size (N*R)=%5.0f',NR_max));
%     for iNR=1:length(N_list)
%         disp(sprintf(' %2.0f-> blocking (N) =%4.0f    blocks per patch edge (R) =%4.0f    patch size (NR)=%4.0f   scale for btc statistics, in blocks(S)=%4.0f step phase: %s',...
%             iNR,N_list(iNR),R_list(iNR),N_list(iNR)*R_list(iNR),S_list(iNR),Sstep_char{Sstep_list(iNR)+1}));
%     end
     if_NRS_ok=getinp('1 if ok','d',[0 1]);
end %if_NRS_ok
%
%divide each roi into bigpatches
%
bigpatches=zeros(bigpatch_size,bigpatch_size,0); % each slice is a square bigpatch of size bigpatch_size
bigpatch_origroi=zeros(0,1); %indicates which roi is the source of each bigpatch
bigpatch_origloc=zeros(0,2); %indicates where the bigpatch is located in the roi
%
NR_min=min(N_list.*R_list);
disp(sprintf('minimum patch size: %4.0f',NR_min));
NR_exclude=NR_min;
if_exclude=getinp('1 to exclude constant patches','d',[0 1],1);
if (if_exclude==1)
    NR_exclude=getinp('patch size to use for excluding','d',[4 NR_max],NR_exclude);
end
n_bigpatches=0;
n_bigpatches_excluded=0;
for i_roi=1:n_rois %look at each bigpatch   
    nbigpatches_xy=floor(roi_sizes(i_roi,:)/bigpatch_size);
    disp(sprintf(' working on roi %3.0f [%4.0f x %4.0f], from original image %3.0f',...
        i_roi,roi_sizes(i_roi,:),roi_origimg(i_roi)));
    row_starts=1+floor(roi_sizes(i_roi,1)*[0:nbigpatches_xy(1)-1]/nbigpatches_xy(1)); 
    row_ends=[row_starts(2:end)-1,roi_sizes(i_roi,1)];
    col_starts=1+floor(roi_sizes(i_roi,2)*[0:nbigpatches_xy(2)-1]/nbigpatches_xy(2));
    col_ends=[col_starts(2:end)-1,roi_sizes(i_roi,2)];
    %choose a random bigpatch position for each row and column of subdivision
    for irow=1:nbigpatches_xy(1)
        for icol=1:nbigpatches_xy(2)
            row_start_range=[row_starts(irow):row_ends(irow)-bigpatch_size+1]; %possible row starting position
            col_start_range=[col_starts(icol):col_ends(icol)-bigpatch_size+1]; %possible col starting position
            row_start=row_start_range(ceil(length(row_start_range)*rand(1)));
            col_start=col_start_range(ceil(length(col_start_range)*rand(1)));
            %disp([irow icol row_start col_start])
            %
            %create a new bigpatch
            bigpatch_data=roi_data_all{i_roi}(row_start:row_start+bigpatch_size-1,col_start:col_start+bigpatch_size-1); %finally cut out a bigpatch
            %test to see whether any patch would be constant
            %subdivide bigpatch_data into a stack of patches of size NR_exclude
            patch_data=permute(reshape(bigpatch_data,[NR_exclude,bigpatch_size/NR_exclude,NR_exclude,bigpatch_size/NR_exclude]),[1 3 2 4]);
            patch_data_min=min(min(patch_data,[],2),[],1);% minimum within each slice
            patch_data_max=max(max(patch_data,[],2),[],1);% maximum within each slice
            if all(patch_data_min(:)~=patch_data_max(:)) | (if_exclude==0)
                n_bigpatches=n_bigpatches+1;
                bigpatches(:,:,n_bigpatches)=bigpatch_data; %finally cut out a bigpatch
                bigpatch_origroi(n_bigpatches)=i_roi; %which roi is the origin of this bigpatch
                bigpatch_origloc(n_bigpatches,:)=[row_start,col_start]; %location of this bigpatch in the original roi
            else
                n_bigpatches_excluded=n_bigpatches_excluded+1;
            end
        end
    end
    disp(sprintf(' inspected %4.0f bigpatches from roi %3.0f',nbigpatches_xy(1)*nbigpatches_xy(2),i_roi));
end
disp(sprintf(' made a total of %5.0f bigpatches, %5.0f bigpatches excluded',n_bigpatches,n_bigpatches_excluded));
%
%ask whether to use rng('default') or rng('shuffle')
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time','d',[0 1]);
if (if_frozen)
    rng('default');
else
    rng('shuffle');
end

%options for power spectrum calculation
wintype_desc={'none','cb1'};
wintype=getinp('window type (0 for none, 1 for cosbell)','d',[0 1],1);
fft2_padfactor=getinp('padding factor for Fourier Transform (power of 2, 1 for no padding)','d',[1 8],2); 
fft2_padfactor=2.^floor(log2(fft2_padfactor));
%options for image statistics calculation
margin_width=getinp('margin size for calculating binary texture statistics','d',[0 floor(min(R_list)/2)-1],1);
%
bigpatches_show=getinp('list of bigpatches to show pipeline (0 for none, can be array)','d',[0 n_bigpatches]);
if max(bigpatches_show)>0
    iNR_show=getinp('setups to show','d',[1 length(N_list)]);
else
    iNR_show=0;
end
%
patches_blockcounts=cell(1,length(N_list));
patches_btc=cell(1,length(N_list));
patch_origroi=cell(1,length(N_list)); %which roi is the origin of this patch
%
for iNR=1:length(N_list)
    N=N_list(iNR);
    R=R_list(iNR);
    disp(sprintf('doing analysis %2.0f-> blocking (N) =%4.0f    blocks per patch edge (R) =%4.0f    patch size (NR)=%4.0f',...
        iNR,N,R,N*R));
    patch_size=N*R;
    patches_per_bigpatch_side=bigpatch_size/patch_size;
    patches_per_bigpatch=patches_per_bigpatch_side.^2;
    n_patches=n_bigpatches*patches_per_bigpatch;
    disp(sprintf(' creating %5.0f patches of size %3.0f from %5.0f  bigpatches of size %3.0f',n_patches,patch_size,n_bigpatches,bigpatch_size));
    %now subdivide bigpatches into patches
    patches=reshape(permute(reshape(bigpatches,[patch_size,patches_per_bigpatch_side,patch_size,patches_per_bigpatch_side,n_bigpatches]),[1 3 2 4 5]),[patch_size patch_size n_patches]);
    %do downsampling
    patches_downsampled=reshape(mean(mean(reshape(patches,[N R N R n_patches]),3),1),[R R n_patches]);
    %determine spectrum  (sub mean, window, pad, FFT, abs square, average)
    win=mlis_window_setup(R,wintype); %window type depends on R
    nfft=R*(fft2_padfactor); %fft size depends on R
    patch_fft=zeros(nfft,nfft,n_patches); %keep transforms of all patches
    windowed=zeros(R,R,n_patches);
    patches_whitened=zeros(R,R,n_patches);
    patches_binarized=zeros(R,R,n_patches);
    meansub_range=zeros(n_patches,2);
    for i_patch=1:n_patches
        patch_origroi{iNR}(i_patch)=bigpatch_origroi(1+floor((i_patch-1)/patches_per_bigpatch));
        meansub=patches_downsampled(:,:,i_patch)-mean(mean(patches_downsampled(:,:,i_patch))); %mean of everything in the patch number i_patch
        meansub_range(i_patch,:)=[min(meansub(:)) max(meansub(:))]; %needed for plotting
        if meansub_range(i_patch,2)<=meansub_range(i_patch,1)
            meansub_range(i_patch,2)=meansub_range(i_patch,1)+1;
        end
        windowed(:,:,i_patch)=meansub.*win;
        patch_fft(:,:,i_patch)=fft2(windowed(:,:,i_patch),nfft,nfft); %NOT shifted
    end
    %average
    pwrspec_unshifted=mean(abs(patch_fft).^2,3); %mean along dimension 3, i.e., across all patches, of magnitude-squared of patch Fourier transforms
    %whiten
    for i_patch=1:n_patches
        patch_whitened_padded=ifft2(patch_fft(:,:,i_patch)./sqrt(pwrspec_unshifted));
        patch_whitened_raw=patch_whitened_padded(1:R,1:R);
        %adjust mean and standard dev to match original patch
        windowed_single=windowed(:,:,i_patch);
        windowed_mean=mean(windowed_single(:));
        windowed_std=std(windowed_single(:));
        %
        patch_whitened_adj=patch_whitened_raw*windowed_std/std(patch_whitened_raw(:)); %adjust the standard dev
        patch_whitened_adj=patch_whitened_adj+windowed_mean-mean(patch_whitened_adj(:)); % adjust the mean
        patches_whitened(:,:,i_patch)=patch_whitened_adj;
    end
    %
    patches_blockcounts{iNR}=zeros(btc_nconfigs,n_patches);
    patches_btc{iNR}=zeros(btc_n,n_patches);
    for i_patch=1:n_patches
        %this section of code adapted from ffdm_btc_calc_gen
        %binarize
        patches_binarized(:,:,i_patch)=patches_whitened(:,:,i_patch)>median(reshape(patches_whitened(:,:,i_patch),[R*R,1]));
        %calculate block counts
        patch_marg=patches_binarized((1+margin_width):(end-margin_width),(1+margin_width):(end-margin_width),i_patch); %remove the margin
        blockcounts=glider_mapubi(patch_marg,btc_checkdef,btc_ng,setfield([],'mapubi_bc',0)); %NON-periodic boundary conditions
        p2x2=reshape(blockcounts,btc_reshape)/sum(blockcounts(:));
        patches_blockcounts{iNR}(:,i_patch)=blockcounts;
        %calculate binary texture coordinates (local image statistics)
        patches_btc{iNR}(:,i_patch)=btc_corrs2vec(getcorrs_p2x2(p2x2,0,1),btc_dict); %compute only btc stats (07Feb20)
    end
    %
    %show pipeline if requested
    nrows=3; %first row for bigpatch, row 2 for first patch in bigpatch, row 3 for last patch in bigpatch
    ncols=5; %stages in pipeline
    if ismember(iNR,iNR_show) 
        for i_bigpatch=bigpatches_show
            i_roi=bigpatch_origroi(i_bigpatch);
            im_file=im_files{roi_origimg(i_roi)};
            name_pipeline=sprintf(' bigpatch %3.0f from roi %3.0f from image %s N=%2.0f R=%4.0f',i_bigpatch,i_roi,im_file,N,R);       
            figure;
            set(gcf,'NumberTitle','off'); %turn off numbered figure
            set(gcf,'Name',name_pipeline);
            set(gcf,'Position',[100 100 1200 700]);
            %
            %show original bigpatch
            %
            subplot(nrows,ncols,1); %original
            imagesc(bigpatches(:,:,i_bigpatch),[0 jpeg_max]);
            colormap gray;
            title(sprintf('bigpatch %3.0f',i_bigpatch));
            axis equal;
            axis tight;
            colorbar;
            %
            %show first and last patches
            %
            for i_firstlast=1:2
                switch i_firstlast
                    case 1
                        patch_offset=1;
                        patch_label='first';
                    case 2
                        patch_offset=patches_per_bigpatch;
                        patch_label='last';
                end
                i_patch=patch_offset+(i_bigpatch-1)*patches_per_bigpatch;
                subplot(nrows,ncols,1+i_firstlast*ncols); %original
                imagesc(patches(:,:,i_patch),[0 jpeg_max]);
                colormap gray;
                title(sprintf('%s patch (%3.0f)',patch_label,i_patch));
                axis equal;
                axis tight;
                colorbar;
                %
                subplot(nrows,ncols,2+i_firstlast*ncols); %downsampled
                imagesc(patches_downsampled(:,:,i_patch),[0 jpeg_max]);
                colormap gray;
                title(sprintf('downsampled by %2.0f',N));
                axis equal;
                axis tight;
                colorbar;
                %
                subplot(nrows,ncols,3+i_firstlast*ncols); %windowed
                imagesc(windowed(:,:,i_patch),meansub_range(i_patch,:));
                colormap gray;
                title(sprintf('windowed (%s)',wintype_desc{wintype}));
                axis equal;
                axis tight;
                colorbar;
                %
                subplot(nrows,ncols,4+i_firstlast*ncols); %whitened
                imagesc(patches_whitened(:,:,i_patch),meansub_range(i_patch,:));
                colormap gray;
                title(sprintf('whitened (pad %2.0f)',fft2_padfactor));
                axis equal;
                axis tight;
                colorbar;
                %
                subplot(nrows,ncols,5+i_firstlast*ncols); %binarized
                imagesc(patches_binarized(:,:,i_patch),[0 1]);
                colormap gray;
                title('binarized');
                axis equal;
                axis tight;
                colorbar;
            end %i_firstlast
            %
            %add a label
            %
            axes('Position',[0.02,0.02,0.01,0.01]); %for text
            text(0,0,name_pipeline,'Interpreter','none');
            axis off;
            drawnow;
        end 
    end %if show
end %iNR
%
%plot heatmaps of binary texture coords
%
name_btc=sprintf('%s: %3.0f database entries, %3.0f rois, %3.0f bigpatches of size %3.0f x %3.0f, pad factor %2.0f, window type %1.0f (%s)',...
    file_name,n_entries,n_rois,n_bigpatches,bigpatch_size,bigpatch_size,fft2_padfactor,wintype,wintype_desc{wintype+1});figure;
set(gcf,'NumberTitle','off'); %turn off numbered figure
set(gcf,'Name',cat(2,'btc ',name_btc));
set(gcf,'Position',[100 100 1200 700]);
for iNR=1:length(N_list)
    subplot(1,length(N_list),iNR);
    imagesc(patches_btc{iNR}',[-1 1]);
    xlabel('coord');
    set(gca,'XTick',[1:btc_n]);
    set(gca,'XTickLabel',btc_dict.codel');
    ylabel('patch number');
    title(sprintf('N=%2.0f R=%3.0f',N_list(iNR),R_list(iNR)));
end
axes('Position',[0.02,0.02,0.01,0.01]); %for text
text(0,0,name_btc,'Interpreter','none');
axis off;
drawnow;
%
disp('done.');
