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
  
  %scrsz = get(0,'ScreenSize') ;
  %paperPosition = [0 0 scrsz(3) - 2 scrsz(4) - 2] ;
  set(gcf, 'PaperPositionMode', 'auto', 'InvertHardCopy', 'off', 'Position', [1 25 1920 950])

  
  tempfile=figfiles{i};

  filename=[tempfile(1:end-4) '.png'];
  pause(.1)
  print('-dpng', '-r72', filename)
  i=i+1;
  close all;
end
end