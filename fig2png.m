function fig2png(ylims, xlims)
if nargin < 1 || isempty(ylims)
  ylims = [] ;
end
if nargin < 2 || isempty(xlims)
  xlims = [] ;
end
dirstruct=dir(pwd);
filenames=arrayfun(@(x) x.name,dirstruct,'UniformOutput',false);
TF=strfind(filenames,'.fig');
notfigcells=cellfun(@isempty,TF);
figfiles_start=(filenames(~notfigcells));
i=1;
figfiles=figfiles_start;
while i<(length(figfiles)+1)
  workingfiles=cellfun(@isempty,figfiles);
  figfiles=figfiles(~workingfiles);
  open(figfiles{i});
  if ~isempty(xlims)
    xlim(xlims)
  end
  if ~isempty(ylims)
    ylim(ylims)
  end
  
  tempfile=figfiles{i};

  filename=[tempfile(1:end-4) '.png'];
  pause(.1)
  print('-dpng', '-r300', filename)
  i=i+1;
  close all;
end
end