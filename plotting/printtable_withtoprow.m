function[] = printtable_withtoprow(M,SD,table_head,row_name,toprow,varargin)
% PRINTTABLE prints the array M in LaTeX table format
% M is an array of size m x n
% SD is an array of size m x n
% In our context, the values in M are averaged performance indicators
% If a standard deviation is available, we print the value in () also
% You must put SD as an array of m x n. If SD is not available just put []
% 
% table_head is a cell of strings containing the column header (size n)
% row_name is a cell of strings containing the row header (size m)
% 
% You can modify the format of printing (float, string, integer, etc in the for loop
% USAGE: 
% printtable(M,SD,table_head,row_name)
% printtable(M,[],table_head,row_name)
    

[m,n] = size(M);

if isempty(varargin)
    text_include = '';
else
n_top = n/length(toprow);
str_in = varargin{1};
text_include = '&';
for k=1:length(toprow)
    if k < length(toprow)
        text_include = [text_include,'\\multicolumn{',int2str(n_top),'}{c}{',str_in,'} & '];
    else
        text_include = [text_include,'\\multicolumn{',int2str(n_top),'}{c}{',str_in,'}'];
    end
end
   
text_include = [text_include,' \\\\ \n \\hline \n'];
end

str_table_head = ' & ';
for k=1:n
    if k < n
        str_table_head = [str_table_head table_head{k} ' & '];
    else
        str_table_head = [str_table_head table_head{k}];
    end
end
str_table_head = [str_table_head ' \\\\ \n'];

n_top = n/length(toprow);
str_top= ['&'];
for k=1:length(toprow)
    if k < length(toprow)
        str_top = [str_top,'\\multicolumn{',int2str(n_top),'}{c}{',toprow{k},'} & '];
    else
        str_top = [str_top,'\\multicolumn{',int2str(n_top),'}{c}{',toprow{k},'}'];
    end
end
str_top = strrep(str_top,'%','\\%%');
% \multicolumn{2}{|c|}{c1} & \multicolumn{2}{c|}{c2} \\ \hline

fprintf('\n');
fprintf(['\\begin{tabular}{' repmat('c',1,n+1) '} \\hline \n']);
fprintf([str_top,' \\\\ \n \\hline \n']);
fprintf(text_include);
fprintf([str_table_head,'\\hline \n']);
if isempty(SD)
        for ii = 1:m
        value_row = [row_name{ii} '& '];
        for jj=1:n
            if jj < n
                value_row = [value_row [num2str(M(ii,jj),'%2.1f'), ' & ']];
            else
                value_row = [value_row [num2str(M(ii,jj),'%2.1f')]];
            end            
        end
        value_row = [value_row ' \\\\ \n'];
        fprintf(value_row);
        end
    
else
        for ii = 1:m
        value_row = [row_name{ii} '& '];
        for jj=1:n
            if jj < n
                value_row = [value_row ['',num2str(M(ii,jj),'%2.1f'),'', ' (','',num2str(SD(ii,jj),'%2.1f'),'', ') & ']];
            else
                value_row = [value_row ['',num2str(M(ii,jj),'%2.1f'),'', ' (','',num2str(SD(ii,jj),'%2.1f'),'', ')']];
            end            
        end
        value_row = [value_row ' \\\\ \n'];
        fprintf(value_row);
        end
end
fprintf(' \\end{tabular} \n')
end
