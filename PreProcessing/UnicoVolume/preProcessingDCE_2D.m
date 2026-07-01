t = readtable('E:\CampusBiomedico\sintesiPazientiFinale.xlsx');
savePath = 'VolumiDCE_Sub_AllBox_2D';
name = t.name;
base = 'E:\CampusBiomedico\Estrazione_NEW';

for i=1:numel(name)
    row = t(i,:);
    disp(['Working on ' row.name{1}]);
    
    %----------------------------------------> determino la cartella di uscita
    if  row.cavo == 1
        specificPath = [savePath filesep '1_cavo'];
    else
        specificPath = [savePath filesep '0_cavo'];
    end
    
    if ~exist(specificPath, 'dir')
        mkdir(specificPath)
    end
    
    %----------------------------------------> load maschera
    niifile = dir([base filesep row.name{1} filesep '*.nii*']);
    maschera = niftiread([base filesep row.name{1} filesep niifile(1).name]);
    
    if strcmp(row.name{1}, '92.Baglioni')
        maschera(:,:,126) = int16(zeros(size(maschera,1), size(maschera,2)));
        maschera(:,:,127) = int16(zeros(size(maschera,1), size(maschera,2)));
        maschera(:,:,128) = int16(zeros(size(maschera,1), size(maschera,2)));
    end
    
    %----------------------------------------> caricamento volumi
    dcemriPath = dir([base filesep row.name{1} filesep 'DCEMRI_*']);
    
    if numel(dcemriPath) > 1
        dcemriPath = dcemriPath(2).name;
    else
        dcemriPath = dcemriPath(1).name;
    end
    
    disp(dcemriPath)
    numeroAcquisizioni = numel(dir([base filesep row.name{1} filesep dcemriPath filesep 'subvt_*'])); %numero di elementi sottrattivi
    disp(['    Extracting from ' dcemriPath ' num ' num2str(numeroAcquisizioni) ]);
    
    %----------------------------------------> selezione delle acquisizioni
    acquisizioni = 1:numeroAcquisizioni;
    disp(['        Acqusition Indeces   ' num2str(acquisizioni) ]);
    
    selectedAcquisition = [acquisizioni(1:2) ceil(median(acquisizioni(3:(end -1)))) acquisizioni(end)];
    disp(['        Selected Indeces   ' num2str(selectedAcquisition) ]);
    
    
    %----------------------------------------> creazione volume
    MRIvolume = int16(zeros(size(maschera,1), size(maschera,2),size(maschera,3), numel(selectedAcquisition)));
    for qq=1:numel(selectedAcquisition)
        filename = ['subvt_' num2str(selectedAcquisition(qq)) '.nii'];
        disp(['             Reading file   ' filename]);
        MRIvolume(:,:,:,qq) = niftiread([base filesep row.name{1} filesep dcemriPath filesep filename]);
    end
    
    %-------------------------------------> divisione dei volumi in base alla posizione della lesione
    [mask_basso, mask_alto] = SplitMaschera(maschera);
    
    totale = mask_basso | mask_alto;
    
    if ~isequal(totale, logical(maschera))
        error('DIVISIONE MASCHERE ERRATA');
    end
    
    %-------------------------------------> Salvataggio dei volumi estratti
    if sum(mask_basso(:)) > 30
        imdb = extractEachLesion(MRIvolume, mask_basso, specificPath, row.name{1}, 'b', dcemriPath, selectedAcquisition);
    end
    
    if sum(mask_alto(:)) > 30
        imdb = extractEachLesion(MRIvolume, mask_alto, specificPath, row.name{1}, 'a', dcemriPath, selectedAcquisition);
    end
    
    
end

function idz= DrawImages(subMri, subMask, subMask2, savepath)
if ~exist(savepath, 'dir')
    mkdir(savepath)
end

subMri = single(subMri);
subMri = (subMri - min(subMri(:)))/(max(subMri(:)) - min(subMri(:)));

[~,~,zz] = ind2sub(size(subMask), find(subMask == 1));
idz = unique(zz);

%dwi e mri comparata prima e dopo
for i=1:numel(idz)
    FigH = figure('Position', get(0, 'Screensize'));
    
    imshow(permute(subMri(:,:,idz(i), [1 2 3]), [1 2 4 3])); title('MRI');
    
    m = subMask(:,:,idz(i));
    [B,L] = bwboundaries(m,'noholes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'y', 'LineWidth', 0.5)
    end
    hold off
    
    m = subMask2(:,:,idz(i));
    [B,L] = bwboundaries(m,'noholes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 0.5)
    end
    hold off
    saveas(FigH, [savepath filesep 'img_' num2str(idz(i)) '.png']);
    
    close all
end
end

function [mask_basso, mask_alto ] = SplitMaschera(subMask)

CC = bwconncomp(subMask);
mask_basso = false(size(subMask));
mask_alto = false(size(subMask));

limite = size(subMask,1)/2;

for j=1:CC.NumObjects
    [xx,yy,zz] = ind2sub(size(subMask), CC.PixelIdxList{j});
    
    if all(xx > limite)
        mask_basso(CC.PixelIdxList{j}) = true;
        
    elseif all(xx <= limite)
        mask_alto(CC.PixelIdxList{j}) = true;
    else
        warning('   ----->  Lesione al centro!!!!');
        %quindi la mettiamo ad entrambi
        mask_basso(CC.PixelIdxList{j}) = true;
        mask_alto(CC.PixelIdxList{j}) = true;
    end
end

mask_alto = int16(mask_alto);
mask_basso = int16(mask_basso);
end

function [x1,x2] = controllo(x1,x2,shapex)
range = x2-x1+1;
if x1<=0
    x2 = x2 + abs(x1) +1;
    x1 = 1;
end


if x2 >shapex
    x1 = x1 - (x2 - shapex);
    x2 = shapex;
end

x1 = round(x1);
x2 = round(x2);

new_range = x2-x1+1;
if range~=new_range
    disp('errore nei controlli');
end

if range > shapex
    x1 = 1;
    x2 = shapex;
end
end

function imdb = extractEachLesion(volumeMRI, mask, specificPath, patientName, sub, dcemriPath, selectedAcquisition)
s= split(patientName,'.');
[shapex, shapey, shapez] = size(mask); %dimensioni della maschera

[xx,yy,zz] = ind2sub(size(mask), find(mask == 1));
xx_unique = unique(xx); % mi serve per il bounding box
yy_unique = unique(yy); % mi serve per il bounding box
zz_unique = unique(zz); % sono solo le slice che contengono lesione

x_min = min(xx_unique);
x_max = max(xx_unique);
ctrx = (x_max + x_min)/2;
rx = (x_max - x_min)/2;

y_min = min(yy_unique);
y_max = max(yy_unique);
ctry = (y_max + y_min)/2;
ry = (y_max - y_min)/2;

rxy = max([rx ry]);

%se è giusto non abbiamo bisogno del round
x_1 = ctrx - rxy;
x_2 = ctrx + rxy;

y_1 = ctry - rxy;
y_2 = ctry + rxy;

disp(['Azze x: ' num2str(x_1) ' ' num2str(x_2) ' Ctrx ' num2str(ctrx)]);
disp(['Azze y: ' num2str(y_1) ' ' num2str(y_2) ' Ctry ' num2str(ctry)]);

%controllo sugli assi
[x_1,x_2] = controllo(x_1,x_2,shapex);
[y_1,y_2] = controllo(y_1,y_2,shapey);

disp(['New Azze x: ' num2str(x_1) ' ' num2str(x_2)]);
disp(['New Azze y: ' num2str(y_1) ' ' num2str(y_2)]);

if (x_2 - x_1) ~= (y_2 - y_1)
    error('Non sono quadrate');
end

for y=1:numel(zz_unique)
    z_to_consider = zz_unique(y);
    
    %box del paziente
    submask = mask(x_1:x_2, y_1:y_2,z_to_consider);
    subVolume = permute(volumeMRI(x_1:x_2, y_1:y_2,z_to_consider,:),[1 2 4 3]);
    pix_only = sum(submask(:));
    
    
    imdb.volume = subVolume;
    imdb.mask = submask;
    imdb.pix_only = pix_only;
    imdb.dimLesione =sum(submask,[1 2]);
    imdb.lesionedSlice = z_to_consider;
    imdb.ZVolOriginal = zz_unique';
    imdb.dcemriPath = dcemriPath;
    imdb.selectedAcquisition = selectedAcquisition;
    save([specificPath filesep s{1} '_' s{2} '_' sub '_' num2str(imdb.dimLesione) '_' num2str(y) '.mat'], '-struct', 'imdb');
end
end