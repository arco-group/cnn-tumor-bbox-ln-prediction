t = readtable('D:\Michela\DOTTORATO_ICTH\Mammella_SODA\Deep_Approach_New\sintesiPazientiFinale.xlsx');
savePath = 'VolumiDCE_Isotropic';
name = t.name;
base = 'D:\Michela\DOTTORATO_ICTH\Mammella_SODA\Deep_Approach_New\EstrazioneDCEMRI\Estrazione_NEW\';

paziente ={};
pixel = [];
raggio = [];

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
    info_mri = niftiinfo([base filesep row.name{1} filesep dcemriPath filesep 'subvt_' num2str(selectedAcquisition(1)) '.nii']);
    new_size = ceil(info_mri.ImageSize .* info_mri.PixelDimensions);
    new_mask = imresize3(maschera, new_size, 'nearest');
    
    MRIvolume = int16(zeros(size(new_mask,1), size(new_mask,2), size(new_mask,3), numel(selectedAcquisition)));
    
    for qq=1:numel(selectedAcquisition)
        filename = ['subvt_' num2str(selectedAcquisition(qq)) '.nii'];
        disp(['             Reading file   ' filename]);
        MRIvolume(:,:,:,qq) = imresize3(niftiread([base filesep row.name{1} filesep dcemriPath filesep filename]), new_size);
    end
    
    if size(MRIvolume(:,:,:,1)) ~= size(new_mask)
        error('Attenzione ai resize!');
    end
    %-------------------------------------> divisione dei volumi in base alla posizione della lesione
    [mask_basso, mask_alto] = SplitMaschera(new_mask);
    
    totale = mask_basso | mask_alto;
    
    if ~isequal(totale, logical(new_mask))
        error('DIVISIONE MASCHERE ERRATA');
    end
    
    %-------------------------------------> Salvataggio dei volumi estratti
    if sum(mask_basso(:)) > 0
        [paziente, pixel, raggio] = extractEachLesion(MRIvolume, mask_basso, specificPath, row.name{1}, 'b', dcemriPath, selectedAcquisition, paziente, pixel, raggio);
    end
    
    if sum(mask_alto(:)) > 0
        [paziente, pixel, raggio] = extractEachLesion(MRIvolume, mask_alto, specificPath, row.name{1}, 'a', dcemriPath, selectedAcquisition, paziente, pixel, raggio);
    end
    
    
end

t = table(paziente', pixel', raggio', 'VariableNames', {'paziente', 'pixel', 'raggio'});
writetable(t, 'InfoAllLesionIsotropic.xlsx');

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

function [paziente, pixel, raggio] = extractEachLesion(volumeMRI, mask, specificPath, patientName, sub, dcemriPath, selectedAcquisition, paziente, pixel, raggio)
s= split(patientName,'.');
pathBase = ['Immagini_Isotropic_' filesep s{1} '_' s{2} '_' sub];

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

z_min = min(zz_unique);
z_max = max(zz_unique); %devo prendere tutto ciò che va da z_min a z_max
ctrz = (z_max + z_min)/2;
rz = (z_max - z_min)/2;

rxyz = max([rx ry rz]);
%controllo sugli assi
x_1 = ctrx - rxyz;
x_2 = ctrx + rxyz;
y_1 = ctry - rxyz;
y_2 = ctry + rxyz;
z_1 = ctrz - rxyz;
z_2 = ctrz + rxyz;

rangeX = x_2-x_1+1;
rangeY = y_2-y_1+1;
rangeZ = z_2-z_1+1;

if (rangeZ > shapez || rangeX > shapex || rangeY> shapey)
    diamToConsider = max([shapex shapey shapez]);
    warning(['Non possiamo usare il massimo raggio ' num2str(rxyz) ' , limitato al massimo range ' num2str(diamToConsider)]);
    rxyz = (diamToConsider - 1)/2;
    x_1 = ctrx - rxyz;
    x_2 = ctrx + rxyz;
    y_1 = ctry - rxyz;
    y_2 = ctry + rxyz;
    z_1 = ctrz - rxyz;
    z_2 = ctrz + rxyz;
end

disp(['Azze x: ' num2str(x_1) ' ' num2str(x_2) ' Ctrx ' num2str(ctrx)]);
disp(['Azze y: ' num2str(y_1) ' ' num2str(y_2) ' Ctry ' num2str(ctry)]);
disp(['Azze z: ' num2str(z_1) ' ' num2str(z_2) ' Ctrz ' num2str(ctrz)]);

[x_1,x_2] = controllo(x_1,x_2,shapex);
[y_1,y_2] = controllo(y_1,y_2,shapey);
[z_1,z_2] = controllo(z_1,z_2,shapez);

disp(['New Azze x: ' num2str(x_1) ' ' num2str(x_2)]);
disp(['New Azze y: ' num2str(y_1) ' ' num2str(y_2)]);
disp(['New Azze z: ' num2str(z_1) ' ' num2str(z_2)]);

%box del paziente
submask = mask(x_1:x_2, y_1:y_2,z_1:z_2);
subVolume = volumeMRI(x_1:x_2, y_1:y_2,z_1:z_2,:);
pix_only = sum(submask,[1,2]);
pix_only = pix_only(:)';

[~,~,newz] = ind2sub(size(submask), find(submask==1));
newz = unique(newz);


if find(pix_only>0) ~= newz'
    error('Problemi con asse z');
end
disp(size(subVolume));
imdb.volume = subVolume;
imdb.mask = submask;
imdb.pix_only = pix_only;
imdb.dimLesione =sum(submask,[1 2 3]);
imdb.lesionedSlice = newz'; %slice con lesione
imdb.ZVolOriginal = zz_unique';
imdb.z_range = z_1:z_2; %estrazione dal volume iniziale
imdb.dcemriPath = dcemriPath;
imdb.selectedAcquisition = selectedAcquisition;
save([specificPath filesep s{1} '_' s{2} '_' sub '_' num2str(imdb.dimLesione) '.mat'], '-struct', 'imdb');

paziente{end+1} = [s{1} '_' s{2} '_' sub];
pixel(end+1) = imdb.dimLesione;
raggio(end+1) = rxyz; 
end