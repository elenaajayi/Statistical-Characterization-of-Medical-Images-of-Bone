%img_gen_cartoons:  cartoon-ish diagrams of image generatation, 
% for projection (like natural images), transmission (for xrays), scatter (for scintigrams)
%
% overlaps and occlusions addressed "by hand" rather than with
% alpha-transparency, etc, so that full occlusion and x-ray-like transparency could be rendered
% and consequently, objects are all rectangular, so that overlaps can be easily computed.
%
% parallel projection only
%
% objects are defined in planes at progressively greater distance from the
% imaaging surface (z=0), with one object in each plane
%
% then the overlaps are calculated.
%
% Loosely inspired by generative_image_model_ideas.doc and medical image analyses
% in ../encl/Abbey, ../Dropbox/BoneStats, ../Dropbox/FromSamXu/Data, ../Dropbox/Textures/RawNIstats/ImagePatchStats
%
rng('default');
%
if ~exist('npxls') npxls=256; end
if ~exist('zoffs') zoffs=[1/40 1/10]*npxls; end %distances of objects from image plane
if ~exist('circ_qual') circ_qual=npxls; end %circule quqlity
if ~exist('zsep') zsep=0.01; end %separation in z to get occlusion rather than interleaving
if ~exist('viewang') viewang=[-45 12]; end %a view angle that shows the separation between the planes
midpt=(npxls+1)/2;
%
%params for scint rendering
%
if ~exist('decays_per_plane') decays_per_plane=7*10^3; end
%
if ~exist('obj_descs')
    obj_descs=cell(0);
    obj_descs{1}.type='rect';
    obj_descs{1}.lengths=[.15 .61]*npxls;
    obj_descs{1}.center=[.10 .18]*npxls;
    obj_descs{1}.color=repmat(.35,1,3);
    %circles don't export well to powerpoint and also harder to compute overlap
    %obj_descs{2}.type='circle';
    %obj_descs{2}.diam=.3*npxls;
    obj_descs{2}.type='rect';
    obj_descs{2}.lengths=[.46*npxls .23*npxls];
    obj_descs{2}.center=[-.13 .12]*npxls;
    obj_descs{2}.color=repmat(.20,1,3);
end
nplanes=length(obj_descs);
%
%compute overlaps
%
xy_lo=NaN(1,2);
xy_hi=NaN(1,2);
for iobj=2:nplanes
    for jobj=1:iobj-1
        ovlp=struct();
        ovlp.type='rect';
        ovlp.color=obj_descs{iobj}.color.*obj_descs{jobj}.color;
        for xy=1:2
            xy_lo(1,xy)=max(obj_descs{iobj}.center(xy)-obj_descs{iobj}.lengths(xy)/2,obj_descs{jobj}.center(xy)-obj_descs{jobj}.lengths(xy)/2);
            xy_hi(1,xy)=min(obj_descs{iobj}.center(xy)+obj_descs{iobj}.lengths(xy)/2,obj_descs{jobj}.center(xy)+obj_descs{jobj}.lengths(xy)/2);
        end
        if all (xy_hi>=xy_lo)
            ovlp.lengths=xy_hi-xy_lo;
            ovlp.center=(xy_hi+xy_lo)/2;
            obj_descs{end+1}=ovlp;
        else
            obj_descs{end+1}=[];
        end
    end
end
nobjs=length(obj_descs);
%
%set up thin objects as masks and as patch data
%
obj_masks=zeros(npxls,npxls,nplanes);
xy_offs=[1:npxls]-(1+npxls)/2;
xy_dists=zeros(npxls,2);
xy_range=NaN(2,2);
obj_patchdata=cell(1,nobjs);
obj_nsegs=zeros(1,nobjs);
obj_areafrac=zeros(1,nobjs); %fraction of the area of the imaging plane
for k=1:nobjs
    if ~isempty(obj_descs{k})
        for xy=1:2
            xy_dists(:,xy)=abs(xy_offs-obj_descs{k}.center(xy));
        end
        switch obj_descs{k}.type
            case 'rect'
                for xy=1:2
                    xy_in=find(xy_dists(:,xy)<=obj_descs{k}.lengths(xy)/2);
                    if ~isempty(xy_in)
                        xy_range(:,xy)=[min(xy_in) max(xy_in)];
                    end
                end
                obj_masks(xy_range(1,1):xy_range(2,1),xy_range(1,2):xy_range(2,2),k)=1;
                obj_patchdata{k}.X=midpt+obj_descs{k}.center(1)+obj_descs{k}.lengths(1)/2*[+1 +1 -1 -1];
                obj_patchdata{k}.Y=midpt+obj_descs{k}.center(2)+obj_descs{k}.lengths(2)/2*[+1 -1 -1 +1];
                obj_areafrac(k)=prod(obj_descs{k}.lengths)/npxls^2;
%             case 'circle'
%                 r2=repmat(xy_dists(:,1).^2,1,npxls)+repmat(xy_dists(:,2).^2,1,npxls)';
%                 obj_masks(:,:,k)=reshape(double(r2<=obj_descs{k}.diam^2/4),[npxls npxls 1]);
%                 obj_patchdata{k}.X=midpt+obj_descs{k}.center(1)+obj_descs{k}.diam*cos(2*pi*[1:circ_qual]/circ_qual)/2;
%                 obj_patchdata{k}.Y=midpt+obj_descs{k}.center(2)+obj_descs{k}.diam*sin(2*pi*[1:circ_qual]/circ_qual)/2;
%                 obj_areafrac{k}=pi*(obj_descs{k}.diam/2)^2/npxls^2;
            otherwise
                warning('unknown object type');
        end
        obj_nsegs(k)=length(obj_patchdata{k}.X);
    end %isempty
end
%
%draw 3D arrangements of patches with various imaging methods
%
nmeths=2;
meth_names={'overview','natural','xray','scint'};
nmeths=length(meth_names);
patch_layout=cell(nmeths,nplanes);
patch_bkgd=cell(nmeths,1);
for imeth=1:nmeths
    figure;
    set(gcf,'Position',[50 50 1200 850]);
    set(gcf,'NumberTitle','off');
    set(gcf,'Name',meth_names{imeth});
    patch_bkgd{imeth}=patch([0 0 npxls npxls],[0 npxls npxls 0],'m');
    hold on;
    set(patch_bkgd{imeth},'FaceColor',[1 1 1]);
    set(patch_bkgd{imeth},'EdgeColor',[0 0 0]);
    for k=1:nplanes
        patch_layout{imeth,k}=patch(obj_patchdata{k}.X,obj_patchdata{k}.Y,repmat(zoffs(k),1,obj_nsegs(k)),'m');
        set(patch_layout{imeth,k},'FaceColor',obj_descs{k}.color);
        set(patch_layout{imeth,k},'EdgeColor',obj_descs{k}.color);
    end
    switch meth_names{imeth}
        case {'natural'}
            for k=1:nplanes
                hp=patch(obj_patchdata{k}.X,obj_patchdata{k}.Y,repmat(zsep*k,1,obj_nsegs(k)),'m');
                set(hp,'FaceColor',obj_descs{k}.color);
                set(hp,'EdgeColor','none');
            end
        case 'xray'
            set(patch_bkgd{imeth},'FaceColor',[0 0 0]);
            for k=1:nobjs %include overlaps
                if obj_nsegs(k)>0
                    hp=patch(obj_patchdata{k}.X,obj_patchdata{k}.Y,repmat(zsep*k,1,obj_nsegs(k)),'m');
                    set(hp,'FaceColor',1-obj_descs{k}.color);
                    set(hp,'EdgeColor','none');
                end
            end
        case 'scint' %simulate decays in random directions
            for k=1:nplanes
                nscints=round(obj_areafrac(k)*decays_per_plane);
                dirs=randn(nscints,3);
                dirs=dirs./repmat(sum(dirs.^2,2),1,3); %normalize
                dirs=dirs(dirs(:,3)~=0,:);
                nscints=size(dirs,1);
                starts=midpt+repmat(obj_descs{k}.center,nscints,1)+repmat(obj_descs{k}.lengths,nscints,1).*(rand(nscints,2)-1/2);           
                dists=zoffs(k)./dirs(:,3); %distances traveled
                endpts=starts(:,[1:2])+repmat(dists,1,2).*dirs(:,1:2);
                hpt=plot3(endpts(:,1),endpts(:,2),zeros(nscints,1),'k.');
            end
    end
    set(gca,'XLim',[0 npxls]);
    set(gca,'YLim',[0 npxls]);
    set(gca,'ZLim',[0 npxls/2]);
    view(viewang);
    axis off;
end
