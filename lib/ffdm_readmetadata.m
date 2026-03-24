function [m,ma,fname_read]=ffdm_readmetadata(filename)
% [m,ma,fname_read]=ffdm_readmetadata(filename) reads the metadata from a csv file
%
% filename: file name (and path)
%  if extension is not given, looks for a *.mat file, if found, reads it and returns the variables m and ma
%  if *.mat is not found, then a csv file is sought, if found, read and parsed
%
% if extension is mat, then only a mat-file is read
% if the extension is csv, then only a csv-file is read
% if no file can be found, fname_read is empty and a warning is issued
%
% m: structure with fields corresponding to columns of filename
%   m.Examp_num is 1 for the first example of a (PID, view) combination, 2 for the second, etc.
% ma: arrays read
% fname_read: file name read
%
%  13Jan20:   Begin modifications to access fields for BIRADS and maker artifacts
%  15Jan20:   Change 'Junk_In_Patch' to 'Quality', with Quality=1 or 2 for Junk_In_Patch=0 or 1
%
%See also:  FFDM_SELECT_PATCHES, FFDM_READ_SUBPATCH, FFDM_SUBPATCH_DEMO.
%
m=[];
ma=[];
fname_read=[];
ext=filename(end-3:end);
%
fields_within_PID={'PID_num','BIRADS','UCD_NUM'}; %these fields, if prseent, should be constant within PID
if strmatch(ext,'.csv')
    if exist(filename,'file')
        [m,ma]=ffdm_readmetadata_csv(filename);
        fname_read=filename;
    else
        warning(sprintf('%s cannot be found.',filename));
    end
elseif strmatch(ext,'.mat') %fixed from 'mat' 13Jan20
    if exist(filename,'file')
        [m,ma]=ffdm_readmetadata_mat(filename);
        fname_read=filename;
    else
        warning(sprintf('%s cannot be found.',filename));
    end
else %extension is not csv or mat, so try adding those
    if exist(cat(2,filename,'.mat'),'file')
        filename_try=cat(2,filename,'.mat');
        [m,ma]=ffdm_readmetadata_mat(filename_try);
        fname_read=filename_try;
    elseif exist(cat(2,filename,'.csv'),'file')
        filename_try=cat(2,filename,'.csv');
        [m,ma]=ffdm_readmetadata_csv(filename_try);
        fname_read=filename_try;
    else
        warning(sprintf('%s, %s.mat and %s.csv cannot be found.',filename,filename,filename));
    end
end
if ~isempty(m)
    [unique_PID,ulist,uptrs]=unique(m.PID);
    uptrs_all=cell(length(unique_PID),1);
    for upid=1:length(unique_PID)
        uptrs_all{upid}=find(uptrs==upid); %points to all entries with each unique patient ID
    end
    for icheck=1:length(fields_within_PID);
        fieldname=fields_within_PID{icheck};
        if isfield(m,fieldname)
            values_inconsistent=[];
            for upid=1:length(unique_PID)
                values_this_pid=m.(fieldname)(uptrs_all{upid});
                if min(values_this_pid)~=max(values_this_pid)
                    values_inconsistent=[values_inconsistent,upid];
                end
            end
            if length(values_inconsistent)==0
                disp(sprintf(' consistency check for %10s      OK: %4.0f unique values across %4.0f patients',...
                    fieldname,length(unique(m.(fieldname))),length(unique_PID)));
            else
                warning(sprintf(' consistency check for %10s   FAILS: %4.0f unique values across %4.0f patients',...
                    fieldname,length(unique(m.(fieldname))),length(unique_PID)));
                for ifail=1:length(values_inconsistent)
                    fail_ptr=values_inconsistent(ifail);
                    disp(sprintf('   PID %3.0f->%s',fail_ptr,unique_PID{fail_ptr}));
                end
            end
        else
            disp(sprintf(' consistency check for %10s SKIPPED: field not found',fieldname));
        end %upid
    end %field to check
end
return

function [m,ma]=ffdm_readmetadata_mat(filename)
s=load(filename);
m=s.m;
ma=s.ma;
return

function [m,ma]=ffdm_readmetadata_csv(filename)
nheader_rows=1;
string_fields={'FName','PID','View'};
headers_rawnames={'Birads_Density','Junk_In_Patch'};
headers_stdnames={'BIRADS','Quality'};
num_fields={'XPos','YPos','XLen','YLen','Age','BIRADS','UCD_NUM','Quality'};
num_fields_temp={'XPos','YPos','XLen','YLen'}; %numeric fields that are can be deleted after reorganization
[num,txt,raw]=xlsread(filename);
ma.num=num;
ma.txt=txt;
ma.raw=raw;
%
nrows=size(raw,1);
ncols=size(raw,2);
nfiles=nrows-nheader_rows;
headers=cell(1,ncols);
for icol=1:ncols
    headers{icol}=strrep(raw{1,icol},' ','');
    %
    rename=strmatch(headers{icol},headers_rawnames,'exact'); %BIRADS became 'Birads_Density' in Annotated_10Jan19
    if ~isempty(rename)
        headers{icol}=headers_stdnames{rename};
    end   
end
%
m.headers=headers;
%string fields
for ifield=1:length(string_fields)
    fieldname=string_fields{ifield};
    col=strmatch(fieldname,headers,'exact');
    s=cell(nfiles,1);
    for ifile=1:nfiles
        s{ifile}=strrep(ma.txt{nheader_rows+ifile,col},' ','');
    end
    m.(fieldname)=s;
end
%numeric fields
for ifield=1:length(num_fields)
    fieldname=num_fields{ifield};
    col=strmatch(fieldname,headers,'exact');
    if ~isempty(col) %added since files prior to Annotated_10Jan19 don't have 'Junk_In_Patch', UCD_NUM; later files don't have 'Age'
        v=zeros(nfiles,1);
        for ifile=1:nfiles
            v(ifile)=ma.raw{nheader_rows+ifile,col};
        end
        m.(fieldname)=v;
    end
end
if isfield(m,'Quality')
    Quality_old=m.Quality;
    m.Quality=zeros(size(Quality_old,1),1);
    m.Quality(Quality_old==0)=1; %first quality is no junk in patch
    m.Quality(Quality_old==1)=2; %second quality is junk in patch
end
%special
PID_num=zeros(nfiles,1);
for ifile=1:nfiles
    PID_num(ifile)=str2num(strrep(strrep(m.PID{ifile},'P',''),' ',''));
end
m.PID_num=PID_num;
m.XYPos=[m.XPos m.YPos];
m.XYLen=[m.XLen m.YLen];
for ifield=1:length(num_fields_temp)
    m=rmfield(m,num_fields_temp{ifield});
end
m.view_names=unique(m.View)'; %this is consistent with view name order in ffdm_extract_demo
View_num=zeros(nfiles,1);
for ifile=1:nfiles
    View_num(ifile)=strmatch(m.View{ifile},m.view_names,'exact');
end
m.View_num=View_num;
m.npatches=nfiles; %total number of patches
%table of how many patches for each patient ID number (row) and view number (column)
m.patch_count_table=full(sparse(m.PID_num,m.View_num,1,max(m.PID_num),max(m.View_num)));
%determine example numbers: 
pid_view=[m.PID_num m.View_num];
[pv_unique,pv_i,pv_j]=unique(pid_view,'rows');
%
examp_num=zeros(size(pid_view,1),1);
for pv_u=1:length(pv_i)
    ptrs=find(pv_j==pv_u);
    examp_num(ptrs)=[1:length(ptrs)];
end
m.Examp_num=examp_num;
return
