function [stimfile, subjectID, pen, site, sort, nreps,nspikes, toes, alltoes] = readtoe(infile)

fid = fopen(infile);
toeinfo = textscan(fid, '%s %s %s %s %s', 1, 'delimiter', '\n');
 stimfile = cell2mat(toeinfo{1});
 subjectID = cell2mat(toeinfo{2});
 pen = cell2mat(toeinfo{3});
 site = cell2mat(toeinfo{4});
 sort = cell2mat(toeinfo{5});
 
 %this is a really bad way of skipping over the "sort : ss001" label,...
 %which Jason does not seem to be keeping track of
 textscan(fid, '%s', 1, 'delimiter', '\n');
 
 nreps = cell2mat(textscan(fid, '%d', 1, 'delimiter', '\n') );
 nspikes = cell2mat(textscan(fid, '%d', nreps, 'delimiter', '\n'));
toestart=ftell(fid);

for i=1:nreps 
 toes{:,i} = textscan(fid, '%.6f32 %*n', nspikes(i) );
end

fseek(fid, toestart, 'bof');
alltoes = cell2mat(textscan(fid, '%.6f32 %*n'));

fclose(fid);
