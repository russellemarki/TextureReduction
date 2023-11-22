clearvars
close all
clc

maxL = 8;

if false
    load('GSH_Library')
else
    [almnsp, blmn, cLin] = calcLibrary(maxL, true);
end


cnt = almnsp_ind(maxL, maxL, maxL, maxL);

cnt - length(almnsp)

atest = cell(1, cnt + 2);

atest{1} = sprintf('static const double almnsp_dat[%i] = {', cnt);
atest{end} = sprintf('};\n');


for i = 1:cnt
    if ~any(almnsp(i)==[0, 1])
        atest{i+1} = sprintf('%.17g,', almnsp(i));
    else
        atest{i+1} = sprintf('%.1f,', almnsp(i));
    end
end

atest{end - 1} = atest{end - 1}(1:end-1);



atest = cat(2, atest{:});

cnt = numel(blmn);

btest = cell(1, cnt + 2);

btest{1} = sprintf('static const double blmn_dat[%i] = {', cnt);
btest{end} = sprintf('};\nreturn blmn_dat[m + %i*(mu-1) + %i*((l/2)-2)];\n', size(blmn, 1), size(blmn, 2)*size(blmn, 1));


for i = 1:cnt
    if ~any(blmn(i)==[0, 1])
        btest{i+1} = sprintf('%.17g,', blmn(i));
    else
        btest{i+1} = sprintf('%.1f,', blmn(i));
    end
end

btest{end - 1} = btest{end - 1}(1:end-1);
btest = cat(2, btest{:});


ctest = cell(1, maxL + 3);
ctest{1} = sprintf('static const size_t cubLin_dat[%i] = {', maxL+1);
ctest{end} = sprintf('};\n');
for i = 1:(maxL+1)
        ctest{i+1} = sprintf('%i,', cLin(i));
end
ctest{end - 1} = ctest{end - 1}(1:end-1);
ctest = cat(2, ctest{:});

clipboard('copy', cat(2, atest, newline, btest, newline, ctest))