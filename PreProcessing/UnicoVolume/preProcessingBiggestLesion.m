t = readtable('E:\CampusBiomedico\sintesiPazientiFinale.xlsx');

% dimx_vet =[];
% dimy_vet=[];
% dimz_vet=[];
% name_lesion = {};

savePath = 'VolumiDCE_Sub_BiggestLesion';
name = t.name;
base = 'E:\CampusBiomedico\Estrazione_NEW';
%numel(name)
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
    if sum(mask_basso(:)) > 0
        imdb = extractLesion(MRIvolume, mask_basso, specificPath, row.name{1}, 'b', dcemriPath, selectedAcquisition);
    end
    
    if sum(mask_alto(:)) > 0
        imdb = extractLesion(MRIvolume, mask_alto, specificPath, row.name{1}, 'a', dcemriPath, selectedAcquisition);
    end
    
    %-------------------------------------> Salvataggio dei volumi estratti
    %     if sum(mask_basso(:)) > 0
    %         [dimx, dimy, dimz] = dimensioneBiggestLesion(mask_basso);
    %         dimx_vet(end+1) =dimx;
    %         dimy_vet(end+1) =dimy;
    %         dimz_vet(end+1) =dimz;
    %         name_lesion{end+1} = [row.name{1} '_b'];
    %     end
    %
    %     if sum(mask_alto(:)) > 0
    %         [dimx, dimy, dimz] = dimensioneBiggestLesion(mask_alto);
    %         dimx_vet(end+1) =dimx;
    %         dimy_vet(end+1) =dimy;
    %         dimz_vet(end+1) =dimz;
    %         name_lesion{end+1} = [row.name{1} '_a'];
    %     end
    
end

%t = table( name_lesion', dimx_vet', dimy_vet', dimz_vet', 'VariableNames', { 'name_lesion', 'dimx_vet', 'dimy_vet', 'dimz_vet'});
%writetable(t, 'InfoBiggestLesion.xlsx');

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
end

function imdb = extractEachLesion(volumeMRI, mask, specificPath, patientName, sub,  x_size, z_size, dcemriPath, selectedAcquisition)
s= split(patientName,'.');
pathBase = ['Immagini_AllBox_BiggestLesion' filesep s{1} '_' s{2} '_' sub];


[shapex, shapey, shapez] = size(mask); %dimensioni della maschera

[xx,yy,zz] = ind2sub(size(mask), find(mask == 1));

xx_unique = unique(xx); % mi serve per il bounding box
yy_unique = unique(yy); % mi serve per il bounding box
zz_unique = unique(zz); % sono solo le slice che contengono lesione

x_min = min(xx_unique);
x_max = max(xx_unique);
ctrx = (x_max + x_min)/2;

y_min = min(yy_unique);
y_max = max(yy_unique);
ctry = (y_max + y_min)/2;

z_min = min(zz_unique);
z_max = max(zz_unique); %devo prendere tutto ciò che va da z_min a z_max
ctrz = (z_max + z_min)/2;

rxy = x_size/2;
rz = z_size/2;

%se è giusto non abbiamo bisogno del round
x_1 = ctrx - rxy;
x_2 = ctrx + rxy -1;

y_1 = ctry - rxy;
y_2 = ctry + rxy-1;

z_1 = ctrz - rz;
z_2 = ctrz + rz -1 ;

disp(['Azze x: ' num2str(x_1) ' ' num2str(x_2) ' Ctrx ' num2str(ctrx)]);
disp(['Azze y: ' num2str(y_1) ' ' num2str(y_2) ' Ctry ' num2str(ctry)]);
disp(['Azze z: ' num2str(z_1) ' ' num2str(z_2) ' Ctrz ' num2str(ctrz)]);

%controllo sugli assi
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

if size(subVolume,1) ~= x_size
    disp(size(subVolume));
    error('Attenzioneeee! Dimensioni Errate');
    
end

if size(subVolume,2) ~= x_size
    disp(size(subVolume));
    error('Attenzioneeee! Dimensioni Errate');
    
end

if size(subVolume,3) ~= z_size
    disp(size(subVolume));
    error('Attenzioneeee! Dimensioni Errate');
    
end

if find(pix_only>0) ~= newz'
    error('Problemi con asse z');
end

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

quadrato = false(size(mask));
quadrato(x_1:x_2, y_1:y_2,z_1:z_2) = 1;
DrawImages(volumeMRI, mask, quadrato, [pathBase filesep 'Standard']);

quadrato = ones(size(submask));
DrawImages(subVolume, submask, quadrato, [pathBase filesep 'Cropped']);

end

function imdb = extractLesion(volumeMRI, mask, specificPath, patientName, sub,dcemriPath, selectedAcquisition)
s= split(patientName,'.');
pathBase = ['Immagini_BiggestLesion' filesep s{1} '_' s{2} '_' sub];
mask = estraiBiggestLesion(mask);

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
z_max = max(zz_unique); 
ctrz = (z_max + z_min)/2;
rz = (z_max - z_min)/2;

rxy = max([rx ry]);
%se è giusto non abbiamo bisogno del round

x_1 = ctrx - rxy;
x_2 = ctrx + rxy;

y_1 = ctry - rxy;
y_2 = ctry + rxy;

z_1 = ctrz - rz;
z_2 = ctrz + rz;

disp(['Azze x: ' num2str(x_1) ' ' num2str(x_2) ' -> Ctrx ' num2str(ctrx)]);
disp(['Azze y: ' num2str(y_1) ' ' num2str(y_2) ' -> Ctry ' num2str(ctry)]);
disp(['Azze z: ' num2str(z_1) ' ' num2str(z_2) ' -> Ctrz ' num2str(ctrz)]);

%controllo sugli assi
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

if size(subVolume,1) ~= size(subVolume,2)
    error('Dimensione x e y non coerenti');
end

if find(pix_only>0) ~= newz'
    error('Problemi con asse z');
end

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

quadrato = false(size(mask));
quadrato(x_1:x_2, y_1:y_2,z_1:z_2) = 1;
DrawImages(volumeMRI, mask, quadrato, [pathBase filesep 'Standard']);

quadrato = ones(size(submask));
DrawImages(subVolume, submask, quadrato, [pathBase filesep 'Cropped']);

end

function [dimx, dimy, dimz] = dimensioneBiggestLesion(mask)
CC = bwconncomp(mask);
vett = [];

for i=1:CC.NumObjects
    vett(end+1) = numel(CC.PixelIdxList{i});
end

[~,idx] = max(vett);

sel_lesion = false(size(mask));
sel_lesion(CC.PixelIdxList{idx}) = true;

[xx,yy,zz] = ind2sub(size(sel_lesion), CC.PixelIdxList{idx});
zz_unique = unique(zz);
xx_unique = unique(xx); % mi serve per il bounding box
yy_unique = unique(yy); % mi serve per il bounding box

%bounding box in x,y,z
x_min = min(xx_unique);
x_max = max(xx_unique);

y_min = min(yy_unique);
y_max = max(yy_unique);

z_min = min(zz_unique);
z_max = max(zz_unique);

if numel([z_min:z_max]) > numel(zz_unique)
    disp('Alcune slice sono vuote');
end
submask = sel_lesion(x_min:x_max, y_min:y_max, z_min:z_max);
dimx = size(submask,1);
dimy = size(submask,2);
dimz = size(submask,3);
end

function sel_lesion = estraiBiggestLesion(mask)
CC = bwconncomp(mask);
vett = [];

for i=1:CC.NumObjects
    vett(end+1) = numel(CC.PixelIdxList{i});
end

[~,idx] = max(vett);

sel_lesion = false(size(mask));
sel_lesion(CC.PixelIdxList{idx}) = true;
end