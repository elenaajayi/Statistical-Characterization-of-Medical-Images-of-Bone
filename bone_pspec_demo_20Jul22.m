%bone_pspec_demo:  compute power spectra from database and its images
%
% 17Jul22
%   * add labels (name_pipeline, name_pspec)
%   * padding: change specification from pad_factor where 0=no padding, 1 = padding  to
%         fft2_padfactor, where 1=no padding, and allow padding of 2, 4, 8, for consistency with ffdm series
%   * fix plot of power spectrum so that y-axis is increases from bottom to top
%
% to do:
%   * more images, 
%    * weighted average so each image counts the same amount?
%
jpeg_max=256; %maximum value in any jpeg file
%
file_name=getinp('bone database file name','s',[0 1],'bonedatabase.xlsx');
[numeric,text_xls,raw]=xlsread(file_name);
n_entries=size(text_xls,1)-1;
disp(sprintf('database has %3.0f entries and a header row',n_entries));
headers=text_xls(1,:);
disp(headers);
%search in headers to get proper columns, in case columns are added at a later time
%
source_col=strmatch('Site',headers); %column for source/site
width_col=strmatch('Width',headers); %column for Width
height_col=strmatch('Height',headers); %column for Height
nroi_col=strmatch('how many ROI?',headers); %column for number of ROIs
orig_col=strmatch('Original Image',headers); %column for image name
type_col=strmatch('Bone type',headers); %column for bone type
%
%go through all entries, reading the ROIs and getting their sizes
%
n_rois=0;
roi_data_all=cell(0);
roi_sizes=zeros(0,2); %roi_sizes(:,1) is height, roi_sizes(:,2) is width
roi_origimg=zeros(0,1); %roi_origimg is original image number
for i_entry=1:n_entries %inspect each entry
    im_file=raw{1+i_entry,orig_col};
    im_size_db=[raw{1+i_entry,height_col} raw{1+i_entry,width_col}]; %size of image (height, width) according to database
    n_roi=raw{1+i_entry,nroi_col};
    disp(sprintf('Requested entry is from %10s, image %15s, size is %5.0f(h) x %5.0f(w); %3.0f ROIs, type: %s',...
        raw{1+i_entry,source_col},im_file,im_size_db,n_roi,raw{1+i_entry,type_col}));
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
    % read the individual rois
    im_file_base=im_file(1:find(im_file=='.')-1); %remove extension from im_file
    im_file_ext=im_file(find(im_file=='.'):end); %just im_file extension
    for i_roi=1:n_roi
        roi_suffix=zpad(i_roi,2); %pad the roi number to two digits
        roi_file=cat(2,im_file_base,'_',roi_suffix,im_file_ext);
        roi_data=imread(roi_file);
        if size(roi_data,3)==3
            roi_data=rgb2gray(roi_data);
        end
        n_rois=n_rois+1; %one more roi
        roi_data_all{n_rois}=roi_data;
        roi_sizes(n_rois,:)=size(roi_data);
        roi_origimg(n_rois)=i_entry;
        disp(sprintf(' roi %2.0f:  (height,width) is %5.0f x %5.0f',i_roi,size(roi_data)));
    end %each ROI
end
minsize=min(roi_sizes(:));
minsize_pwr2=2^floor(log2(minsize));
disp(sprintf('minimum roi size, minimum power-of-2: %5.0f %5.0f',minsize,minsize_pwr2));
%
%ask for patch size (suggest minsize_pwr2, but need not be)
%
ifok=0;
while (ifok==0)
    patch_size=getinp('patch size, should be power of 2','d',[32 1024],minsize_pwr2);
    ifok=(patch_size==2.^floor(log2(patch_size)));
end
%
%now divide each roi into patches
patches=zeros(patch_size,patch_size,0); % each slice is a square patch of size patch_size
patch_origroi=zeros(0,1); %indicates which roi is the source of each patch
patch_origloc=zeros(0,2); %indicates where the patch is located in the roi
%
%ask whether to use rng('default') or rng('shuffle')
if_frozen=getinp('1 for frozen random numbers, 0 for new random numbers each time','d',[0 1]);
if (if_frozen)
    rng('default');
else
    rng('shuffle');
end
%
n_patches=0;
for i_roi=1:n_rois %look at each patch   
    npatches_xy=floor(roi_sizes(i_roi,:)/patch_size);
    disp(sprintf(' working on roi %3.0f [%4.0f x %4.0f], from original image %3.0f',...
        i_roi,roi_sizes(i_roi,:),roi_origimg(i_roi)));
    row_starts=1+floor(roi_sizes(i_roi,1)*[0:npatches_xy(1)-1]/npatches_xy(1)); 
    row_ends=[row_starts(2:end)-1,roi_sizes(i_roi,1)];
    col_starts=1+floor(roi_sizes(i_roi,2)*[0:npatches_xy(2)-1]/npatches_xy(2));
    col_ends=[col_starts(2:end)-1,roi_sizes(i_roi,2)];
    %choose a random patch position for each row and column of subdivision
    for irow=1:npatches_xy(1)
        for icol=1:npatches_xy(2)
            row_start_range=[row_starts(irow):row_ends(irow)-patch_size+1]; %possible row starting position
            col_start_range=[col_starts(icol):col_ends(icol)-patch_size+1]; %possible col starting position
            row_start=row_start_range(ceil(length(row_start_range)*rand(1)));
            col_start=col_start_range(ceil(length(col_start_range)*rand(1)));
            %disp([irow icol row_start col_start])
            %
            %create a new patch
            n_patches=n_patches+1;
            patches(:,:,n_patches)=roi_data_all{i_roi}(row_start:row_start+patch_size-1,col_start:col_start+patch_size-1); %finally cut out a patch
            patch_origroi(n_patches)=i_roi; %which roi is the origin of this patch
            patch_origloc(n_patches,:)=[row_start,col_start]; %location of this patch in the original roi
        end
    end
    disp(sprintf(' made %4.0f patches from roi %3.0f',npatches_xy(1)*npatches_xy(2),i_roi));
end
disp(sprintf(' made a total of %5.0f patches',n_patches));
%do mean subtraction
meansub=zeros(patch_size,patch_size,n_patches); 
for i_patch=1:n_patches
    meansub(:,:,i_patch)=patches(:,:,i_patch)-mean(mean(patches(:,:,i_patch))); %mean of everything in the patch number i_patch
end
%this is another way to do the same thing, without a for loop
%   meansub=patches-repmat(mean(mean(patches,1),2),[patch_size patch_size 1]);
meansub_range=[min(meansub(:)) max(meansub(:))];
%
wintype_desc={'none','cb1'};
wintype=getinp('window type (0 for none, 1 for cosbell)','d',[0 1],1);
win=mlis_window_setup(patch_size,wintype); %make a window
%multiply by window
windowed=zeros(patch_size,patch_size,n_patches); 
for i_patch=1:n_patches
    windowed(:,:,i_patch)=meansub(:,:,i_patch).*win;
%    meansub(:,:,i_patch)=patches(:,:,i_patch)-mean(mean(patches(:,:,i_patch))); %mean of everything in the patch number i_patch
end
%
i_patch=-1;
while (i_patch~=0)
    i_patch=getinp('patch to show pipeline (0 to end)','d',[0 n_patches]);
    if (i_patch>0)
        i_roi=patch_origroi(i_patch);
        im_file=raw{1+roi_origimg(i_roi),orig_col};
        name_pipeline=sprintf(' patch %3.0f from roi %3.0f from image %s',i_patch,i_roi,im_file);       
        figure;
        set(gcf,'NumberTitle','off'); %turn off numbered figure
        set(gcf,'Name',name_pipeline);
        set(gcf,'Position',[100 100 1200 700]);
        %show original patch
        subplot(1,3,1); %original
        imagesc(patches(:,:,i_patch),[0 jpeg_max]);
        colormap gray;
        title(name_pipeline);
        axis equal;
        axis tight;
        colorbar;
        %show mean-subtracted patch
        subplot(1,3,2); %mean-sub
        imagesc(meansub(:,:,i_patch),meansub_range);
        colormap gray;
        title('mean subtracted');
        axis equal;
        axis tight;
        colorbar;
        %show windowed patch
        subplot(1,3,3); %windowed
        imagesc(windowed(:,:,i_patch),meansub_range);
        colormap gray;
        title('windowed');
        axis equal;
        axis tight;
        colorbar;
        %add a label
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,name_pipeline,'Interpreter','none');
        axis off;
        drawnow;
    end
end
%
% do Fourier transform and shift
%
fft2_padfactor=getinp('padding factor for Fourier Transform (power of 2, 1 for no padding)','d',[1 8],2); 
fft2_padfactor=2.^floor(log2(fft2_padfactor));
%
nfft=patch_size*(fft2_padfactor);
patch_fftsq=zeros(nfft,nfft,n_patches); 
for i_patch=1:n_patches
    patch_fftsq(:,:,i_patch)=fftshift(abs(fft2(windowed(:,:,i_patch),nfft,nfft)).^2);
end
%average
pwrspec=mean(patch_fftsq,3); %mean along dimension 3, i.e., across all patches
%plot
name_pspec=sprintf(' %s: %2.0f database entries, %3.0f rois, %3.0f patches of size %3.0f x %3.0f, pad factor %2.0f, window type %1.0f (%s)',...
    file_name,n_entries,n_rois,n_patches,patch_size,patch_size,fft2_padfactor,wintype,wintype_desc{wintype+1});
figure;
set(gcf,'NumberTitle','off'); %turn off numbered figure
set(gcf,'Name',name_pspec);
set(gcf,'Position',[100 100 1200 700]);
zerofreq=nfft/2+1;
imagesc(log10(pwrspec));
set(gca,'YDir','normal'); %17Jul21 added so that y-axis goes from low to high
title('power spectrum')
tick_fracs=[-0.75:.25:.75]; %tick locations relative to center
tick_locs=zerofreq+(fft2_padfactor)*patch_size/2*tick_fracs;
cycles_per_patch=tick_fracs*(patch_size/2); %patch_size/2 is the Nyquist freq
cycles_per_pixel=cycles_per_patch/patch_size;
tick_vals=cycles_per_pixel;
set(gca,'XTick',tick_locs);
set(gca,'XTickLabel',tick_vals);
xlabel('cycles per pixel');
set(gca,'YTick',tick_locs);
set(gca,'YTickLabel',tick_vals);
ylabel('cycles per pixel');
axis equal;
axis tight;
colorbar;
% add a label
axes('Position',[0.02,0.02,0.01,0.01]); %for text
text(0,0,name_pspec,'Interpreter','none');
axis off;
drawnow;
%
%show power spectrum as function of frequency magnitude as scattergram
%
%create an array of frequency magnitudes, of size equal to that of pwrspec
%
freqs=([0:nfft-1]-nfft/2)/nfft; %cycles per pixel
fsq_1d=(repmat(freqs,nfft,1)).^2; %freq squared
fsq=fsq_1d+fsq_1d'; %fx2+fy2
fmag=sqrt(fsq);
%
%select based on frequencies that are less than Nyquiat, i.e., 0.5 cycle/pixel
%
fselect=(fmag<0.5); 
%
figure; %scattergram of power spectrum as function of frequency
set(gcf,'Position',[100 100 1200 700]);
%
subplot(1,2,1);
loglog(fmag(fselect(:)),pwrspec(fselect(:)),'k.');
xlabel('frequency, cycles/pixel');
ylabel('power spectrum');
%add lines of specific slopes
%
hold on;
%ymid=geomean(get(gca,'YLim')); %get a typical y-value
ymid=exp(mean(log(get(gca,'YLim')))); %get a typical y-value
xslope=[.01 .1]; %x-values plotted for standard slopes
plot(xslope,[ymid ymid/10^2],'b');
plot(xslope,[ymid ymid/10^3],'m');
plot(xslope,[ymid ymid/10^4],'r');
legend({'power spectrum','slope -2','slope -3','slope -4'},'Location','NorthEast');
%
% plot specific slices as scattergram
%    
subplot(1,2,2);
loglog(fmag(nfft/2+1,:),pwrspec(nfft/2+1,:),'r.');
hold on;
loglog(fmag(:,nfft/2+1),pwrspec(:,nfft/2+1),'g.');
loglog(diag(fmag),diag(pwrspec),'bo');
loglog(diag(fliplr(fmag)),diag(fliplr(pwrspec)),'bx');
legend({'0 deg','90 deg','oblique','oblique'},'Location','NorthEast');
% add a label
axes('Position',[0.02,0.02,0.01,0.01]); %for text
text(0,0,name_pspec,'Interpreter','none');
axis off;
drawnow;
%
%now use power spectrum to "whiten" the original image patches -- 
%
pwrspec_unshifted=fftshift(pwrspec);
patch_whitened=zeros(patch_size,patch_size,n_patches); 
for i_patch=1:n_patches
    patch_fft=fft2(windowed(:,:,i_patch),nfft,nfft);% Fourier transform of the windowed patch
    patch_whitened_padded=ifft2(patch_fft./sqrt(pwrspec_unshifted));
    patch_whitened_raw=patch_whitened_padded(1:patch_size,1:patch_size);
    %adjust mean and standard dev to match original patch
    windowed_single=windowed(:,:,i_patch);
    windowed_mean=mean(windowed_single(:));
    windowed_std=std(windowed_single(:));
    %
    patch_whitened_adj=patch_whitened_raw*windowed_std/std(patch_whitened_raw(:)); %adjust the standard dev
    patch_whitened_adj=patch_whitened_adj+windowed_mean-mean(patch_whitened_adj(:)); % adjust the mean
    patch_whitened(:,:,i_patch)=patch_whitened_adj;
end
i_patch=-1;
while (i_patch~=0)
    i_patch=getinp('patch to show whitened version (0 to end)','d',[0 n_patches]);
    if (i_patch>0)
        i_roi=patch_origroi(i_patch);
        im_file=raw{1+roi_origimg(i_roi),orig_col};
        name_pipeline=sprintf(' patch %3.0f from roi %3.0f from image %s',i_patch,i_roi,im_file);       
        figure;
        set(gcf,'NumberTitle','off'); %turn off numbered figure
        set(gcf,'Name',name_pipeline);
        set(gcf,'Position',[100 100 1200 700]);
        %show original patch
        subplot(1,3,1); %original
        imagesc(patches(:,:,i_patch),[0 jpeg_max]);
        colormap gray;
        title(name_pipeline);
        axis equal;
        axis tight;
        colorbar;
        %show windowed patch
        subplot(1,3,2); %windowed
        imagesc(windowed(:,:,i_patch),meansub_range);
        colormap gray;
        title('windowed');
        axis equal;
        axis tight;
        colorbar;
        %add a label
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,name_pipeline,'Interpreter','none');
        axis off;
        drawnow;
        %show whitened patch
        subplot(1,3,3); %windowed
        imagesc(patch_whitened(:,:,i_patch),meansub_range);
        colormap gray;
        title('windowed and whitened');
        axis equal;
        axis tight;
        colorbar;
        %add a label
        axes('Position',[0.02,0.02,0.01,0.01]); %for text
        text(0,0,name_pipeline,'Interpreter','none');
        axis off;
        drawnow;

    end
end

disp('done.');

