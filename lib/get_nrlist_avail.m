function [nrlist_avail,NR_ptrs]=get_nrlist_avail(N_list,R_list,corr_fac,NR_nominal)
% [nrlist_avail,NR_ptrs]=function [nrlist_avail,NR_ptrs]=get_nrlist_avail(N_list,R_list,corr_fac,NR_nominal)
% creates a structure nrlist_avail whose fields indicate the values of N (downsamplings) and R (region size, in downsampled checks)
% and, if called with at least 2 arguments, selects from a set of available values of N and R. In this case, user has option of free choice or
% selection of options that intersect with (N,R) values from field of nrlist_avail
%
%   input:
%       * if none are supplied, then only nrlist_avail is returned
%       * if N_list is supplied but empty, then a nice display of nrlist_avail is also produced
%       * if N_list and R_list are both supplied (must be same length), then NR_ptrs is determined from keyboard inputs 
%     N_list: available values of N, downsampling, not including preprocessing
%     R_list: available values of R, region size, in downsampled checks
%       Note that (N_list,R_list) cannot contain duplicate pairs.
%     corr_fac: correction factor to take into account downsampling at preprocessing stage,
%        typically equal to preprocess_downsample from ffdm_btc_calc_gen; 1 if not supplied or empty
%     NR_nominal: nominal NR value (number of pixels across patch) for UPMC data (patches are irregular), 512 if not supplied
%   output: 
%     nrlist_avail: values of N and R availaable in full analyses via ffdm_btc_calc_gen and nistats.mat, with fields:
%       ffdm_unmod: unmodified ffdm, src=2 
%       ffdm_mlis: ffdm with modified local image statistics, src=3
%       usic: carotid ultrasound, src=4
%       mrix_[adni|oasis|hv|patient]: mri: adni, oasis, tns healthy volunteers, tns patients), src=5
%       nistats: natural images, src=6
%       upmc_imgs: Univ of Pittsburgh Med Ctr mammograms, src=7
%       upmc_nvms: Univ of Pittsburgh Med Ctr multivariate Gaussian surrogates, src=8
%       nibs: natural image database with parsed boundary segments, src=9
%    each of these has fields N_list, R_list, NR_list, N_list_corr NR_list_corr (corrected for preprocessing)
%        src (source, integer, corresponding to ffdm_btc_plot_[gen|nistats])
%        src_label (string, e.g., 'ffdm_unmod')
%     NR_ptrs:  selected pointers into N_list, R_list (only if those args are present)
%
% 28Oct21: add upmc (source 7) and nvms (source 8) and NR_nominal
% 29Jun22: add nibs (source 9), fix documentation
%
%    See also:  FFDM_BTC_CALC_GEN, FFDM_BTC_SPEC_GEN, FFDM_BTC_PLOT_GEN,
%       FFDM_BTC_PLOT_NISTATS, FFDM_NRPLOT, FFDM_UPMC_NVMS_CALC.
%
if (nargin<=2)
    corr_fac=[];
end
if isempty(corr_fac)
    corr_fac=1;
end
if (nargin<=3)
    NR_nominal=512;
end
%range for full analysis, must be powers of 2
N_range=[1,32];
R_range=[16,512];
%max region sizes
usic_size=256;
mrix_size=64;
%
ffdm_unmod_size=512;
ffdm_mlis_size=256;
%
nibs_size=256; %nominal patch size for boundary database
%
% downsamping at preprocessing stage
preprocess_downsamples=[1 1 2 1 1 1 1 1 1]; %2 for src=3 (ffdm_mlis), otherwise 1
nsrcs=length(preprocess_downsamples);
% source labels
src_labels={'ffdm_pilot','ffdm_unmod','ffdm_mlis','usic','mrix','nistats','upmc_imgs','upmc_nvms','nibs'};
%
nrlist_avail=struct();
NR_ptrs=[];
%
nrlist_avail.ffdm_unmod=allnr(ffdm_unmod_size,N_range,R_range);
nrlist_avail.ffdm_unmod.src=2;
%
nrlist_avail.ffdm_mlis=allnr(ffdm_mlis_size,N_range,R_range);
nrlist_avail.ffdm_mlis.src=3;
%
nrlist_avail.usic=allnr(usic_size,N_range,R_range);
nrlist_avail.usic.src=4;
%
mrix=allnr(mrix_size,N_range,R_range);
mrix.src=5;
nrlist_avail.mrix_adni=mrix;
nrlist_avail.mrix_oasis=mrix;
nrlist_avail.mrix_hv=mrix;
nrlist_avail.mrix_patient=mrix;
%
nrlist_avail.nistats.N_list=[1    1    1    1    1    2    2    2    2    4    4    4    4    8    8   12   16   20]; % from Hermundstad NIstats.mat
nrlist_avail.nistats.R_list=[128 32   48   64   80   32   48   64   80   32   48   64   80   32   48   32   32   32]; % from Hermundstad NIstats.mat
nrlist_avail.nistats.src=6;
%
upmc_nscales=6;
upmc.N_list=2.^[0:upmc_nscales-1];
upmc.R_list=NR_nominal./upmc.N_list;
nrlist_avail.upmc_imgs=upmc;
nrlist_avail.upmc_imgs.src=7;
nrlist_avail.upmc_nvms=upmc;
nrlist_avail.upmc_nvms.src=8;
%
nrlist_avail.nibs=allnr(nibs_size,N_range,R_range);
nrlist_avail.nibs.src=9;
%
%add other fields N_corr,  and NR_list fields
fns=fieldnames(nrlist_avail);
n_r_corrs_list=cell(1,nsrcs); %list of (n,r)-corrected pairs
for ifn=1:length(fns)
    fn=fns{ifn};
    src=nrlist_avail.(fn).src;
    preprocess_downsample=preprocess_downsamples(src);
    nrlist_avail.(fn).src_label=src_labels{src};
    nrlist_avail.(fn).NR_list=nrlist_avail.(fn).N_list.*nrlist_avail.(fn).R_list;
    nrlist_avail.(fn).N_list_corr=nrlist_avail.(fn).N_list*preprocess_downsample;
    nrlist_avail.(fn).NR_list_corr=nrlist_avail.(fn).NR_list*preprocess_downsample;
    if isempty(n_r_corrs_list{src})
        n_r_corrs_list{src}=[nrlist_avail.(fn).N_list_corr'  nrlist_avail.(fn).R_list'];
    end
end
%
if (nargin==0)
    return
else %determine pointers from input list
    col_header='N(uncorr)         R         NR(uncorr)       N(corr)        NR(corr)';
    col_fmt=repmat('   %4.0f to %4.0f',1,5);
    srcs_have=[];
    disp('analysis parameters in previously-analyzed datasets:');
    disp(cat(2,'  src      label        dataset  nsets       ',col_header));
    for isrc=1:length(preprocess_downsamples)
        for ifn=1:length(fns)
            fn=fns{ifn};
            if (nrlist_avail.(fn).src==isrc)
                srcs_have=[srcs_have,isrc];
                disp(sprintf(cat(2,' %2.0f    %-10s      %-12s %2.0f',col_fmt),...
                    isrc,nrlist_avail.(fn).src_label,fn,length(nrlist_avail.(fn).N_list),...
                    min(nrlist_avail.(fn).N_list),max(nrlist_avail.(fn).N_list),...
                    min(nrlist_avail.(fn).R_list),max(nrlist_avail.(fn).R_list),...
                    min(nrlist_avail.(fn).NR_list),max(nrlist_avail.(fn).NR_list),...
                    min(nrlist_avail.(fn).N_list_corr),max(nrlist_avail.(fn).N_list_corr),...
                    min(nrlist_avail.(fn).NR_list_corr),max(nrlist_avail.(fn).NR_list_corr)));
            end
        end %ifn
    end %isrc
    srcs_have=unique(srcs_have);
end
choices=length(N_list);
%if isempty(choices)
if isempty(N_list)
    return
end
N_list=N_list(:);
R_list=R_list(:);
NR_list=N_list.*R_list;
N_list_corr=N_list*corr_fac;
NR_list_corr=NR_list*corr_fac;
ifok=0;
while (ifok==0)
    disp('available analysis parameters in this dataset')
    disp(cat(2,'                                  set        ',col_header));
    choice_fmt=repmat('       %4.0f    ',1,5);
    for ic=1:choices
        disp(sprintf(cat(2,'                                 %3.0f->',choice_fmt),...
            ic,N_list(ic),R_list(ic),NR_list(ic),N_list_corr(ic),NR_list_corr(ic)));
    end
    how_choose=getinp('1 to enter set choices directly, 2 for intersection of analyzed datasets, 3 for union of analyzed datsets','d',[1 3]);
    NR_ptrs=[];
    n_r_corrs_have=[N_list_corr,R_list];
    switch how_choose
        case 1
            NR_ptrs=getinp('one or more choices of sets','d',[1 choices]);
            NR_ptrs=NR_ptrs(:);
        case 2
            logic_string='intersect';
        case 3
            logic_string='union';
    end
    if isempty(NR_ptrs)
        srcs_logic=getinp(sprintf('list of sources to %s',logic_string),'d',[min(srcs_have),max(srcs_have)]);
        srcs_logic=intersect(srcs_logic,srcs_have); %can only use sources that we actually have
        %operate on (N,R) pairs of corrected values
        for is=1:length(srcs_logic)
            src=srcs_logic(is);
            if (is==1)
                logic_list=n_r_corrs_list{src};
            else
                switch how_choose
                    case 2 %intersect
                        logic_list=intersect(logic_list,n_r_corrs_list{src},'rows');
                    case 3 %union
                        logic_list=union(logic_list,n_r_corrs_list{src},'rows');
                end
            end
            disp(sprintf(' in source %2.0f (%s), N_corr and R are:',src,src_labels{src}));
            disp(n_r_corrs_list{src}');
        end
        disp(sprintf(' after %s, N_corr and R are:',logic_string));
        disp(logic_list');
        %now find pointers to what is in common with N_list_corr and R_list
        n_r_corrs_selected=intersect(n_r_corrs_have,logic_list,'rows');
        disp(sprintf(' values in common with supplied N_corr and R are:'));
        disp(n_r_corrs_selected');
        for is=1:size(n_r_corrs_selected,1)
            NR_ptrs(is)=find(all((n_r_corrs_have==repmat(n_r_corrs_selected(is,:),[size(n_r_corrs_have,1) 1])),2));
        end
    end %NR_ptrs filled by logic
    NR_ptrs=sort(NR_ptrs(:));
    disp('chosen analysis parameters from this dataset')
    disp(cat(2,'                                  set        ',col_header));
    for ic=1:length(NR_ptrs)
        iptr=NR_ptrs(ic);
        disp(sprintf(cat(2,'                                 %3.0f->',choice_fmt),...
            iptr,N_list(iptr),R_list(iptr),NR_list(iptr),N_list_corr(iptr),NR_list_corr(iptr)));
    end
    %disp([NR_ptrs,N_list(NR_ptrs),R_list(NR_ptrs),n_r_corrs_have(NR_ptrs,:)]);
    ifok=getinp('1 if ok','d',[0 1],~isempty(NR_ptrs));
end %ifok
return

function s=allnr(npxls,N_range,R_range)
%fill s.N_list and s.R_list with all power-of-two values that multiply to
%less than npxls, and in which N<=N_max and R>=R_min
N_min=N_range(1);
N_max=N_range(2);
R_min=R_range(1);
R_max=R_range(2);
npxls_pwr2=2.^floor(log(npxls)/log(2));
pxls=min(R_max*N_max,npxls_pwr2);
N_list=[];
R_list=[];
while (pxls>=R_min*N_min)
    N=N_min;
    while (N<=N_max) & (pxls/N)>=R_min
        if (pxls/N<=R_max)
            N_list=[N_list,N];
            R_list=[R_list,pxls/N];
        end
        N=N*2;
    end
    pxls=pxls/2;
end
s.N_list=N_list;
s.R_list=R_list;
return
