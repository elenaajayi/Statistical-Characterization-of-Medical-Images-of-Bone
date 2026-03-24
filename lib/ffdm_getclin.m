function [s,read_info]=ffdm_getclin(filename,col_nums)
% [s,read_info]=ffdm_getclin(filename,col_nums) gets a structure of the clinical data for a mammogram database
%
% filename: xls file name, 'UCD_Mammo_Data.xlsx' if not supplied
% col_nums: a structure of the data to be extracted.  values from UCD_Mammo_Data used if not supplied.
%    col_nums.subj_num is the column number of my subject number, and must be supplied.
%
% s: a structure of the clinical data
% read_info: the informatoin used for reading
%    filename_used: filename used
%    xlsread: arguments returned by xlsread
%    colnums_used: column numbers of the fields
%
% see UCD_Mammo_Data_email_11May18.pdf
%
% See also:  FFDM_TEST, FFDM_EXTRACT_DEMO, FFDM_BTC_SCAT_GEN.
%
if nargin<1
    filename=[];
end
if nargin<2
    col_nums=[];
end
if isempty(filename)
    filename='UCD_Mammo_Data.xlsx';
end
if isempty(col_nums)
    col_nums.subj_num=1;
    col_nums.UCD_ID=2;
    col_nums.BIRADS=3;
    col_nums.VGF_L=4;
    col_nums.VGF_R=5;
end
read_info.filename_used=filename;
[num,txt,raw]=xlsread(filename);
read_info.xlsread.num=num;
read_info.xlsread.txt=txt;
read_info.xlsread.raw=raw;
read_info.col_nums=col_nums;
s=[];
subj_nums=num(:,col_nums.subj_num);
fields=fieldnames(col_nums);
for iv=1:length(fields)
    fn=fields{iv};
    if (strcmp(fn,'subj_num')==0)
        s.(fn)(subj_nums,1)=num(:,col_nums.(fn));
    end
end
subj_missing=setdiff(1:max(subj_nums),subj_nums);
read_info.subj_missing=subj_missing;
if ~isempty(subj_missing)
    warning(sprintf('At least one subject number is missing from %s.',filename))
end
return
