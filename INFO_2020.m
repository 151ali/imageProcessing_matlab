% Test Traitement d'images

IM = imread('IM_03.jpg');
IM = rgb2gray(IM);
%% Egalisation 
Im1 = double(IM);
Im1 = 255*(Im1-min(Im1(:)))/(max(Im1(:))-min(Im1(:)));
Im1 = uint8(Im1);
%% Modification Histogramme
Imax = 128;
Imin = 64;
Im2 = double(IM);
Im2 = (Im2-min(Im2(:)))/(max(Im2(:))-min(Im2(:)));
Im2 = (Imax-Imin)*Im2 + Imin;
Im2 = uint8(Im2);

subplot(3,2,1);
imshow(IM);
subplot(3,2,2);
imhist(IM,255);

subplot(3,2,3);
imshow(Im1);
subplot(3,2,4);
imhist(Im1,255);

subplot(3,2,5);
imshow(Im2);
subplot(3,2,6);
imhist(Im2,255);

%% Négatif d'une image
figure;
imshow(255-IM)

%% Filtrage des images
IM = imread('IM_04.jpg');
IM = rgb2gray(IM);
IM = double(IM);
Hv = (1/6)*[-1 -1 -1;0 0 0; 1 1 1];
Hh = Hv';
Hd1 = (1/(4*sqrt(2)))*[0 0 0 0 -1;0 0 0 -1 0;0 0 0 0 0;0 1 0 0 0;1 0 0 0 0];
Hd2 = -1*flipud(Hd1);
Im1 = abs(imfilter(IM,Hv));
Im2 = abs(imfilter(IM,Hh));
Im3 = abs(imfilter(IM,Hd1));
Im4 = abs(imfilter(IM,Hd2));
Im5 = sqrt(0.25*(Im1.^2 + Im2.^2 + Im1.^3 + Im2.^4));
Im6 = sqrt(0.5*(Im1.^2 + Im2.^2));

figure;
subplot(1,3,1);
imagesc(Im5); colormap gray
subplot(1,3,2);
imagesc(Im6); colormap gray
subplot(1,3,3);
imagesc(max(Im3,2*mean(Im5(:)))); colormap gray


%% Filtre passe bas idéal
IM = imread('IM_04.jpg');
IM = rgb2gray(IM);
IM = double(IM);
N = (-63:63)/3;
M = (-63:63)/3;
[N,M] = meshgrid(N,M);
H = sinc(N).*sinc(M);
figure;
mesh(H);
figure;
imagesc(abs(fftshift(fft2(H))))
IM = imread('IM_04.jpg');
IM = rgb2gray(IM);
IM = double(IM);
Im1 = (imfilter(IM,H));
 imagesc(Im1); colormap gray; axis equal

 %% Filtre passe bas non idéal
 IM = imread('IM_04.jpg');
IM = rgb2gray(IM);
IM = double(IM);
N = (-63:63);
M = (-63:63);
[N,M] = meshgrid(N,M);
H = exp(-0.5*(N.^2 + M.^2)/(1^2));
HF = fftshift(fft2(H));
HF = HF/max(abs(HF(:)));
HF = (1-abs(HF));
HH = ifft2(ifftshift(HF)); 

figure;
subplot(2,1,1);
mesh(HH);
subplot(2,1,2);
mesh(abs(HF))
IM = imread('IM_04.jpg');
IM = rgb2gray(IM);
IM = double(IM);
Im1 = (imfilter(IM,HH));
figure;
 imagesc(IM); colormap gray; axis equal
figure;
 imagesc(Im1); colormap gray; axis equal




