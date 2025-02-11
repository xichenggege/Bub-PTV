function output = denoiser(image)

% Divide photo into two halves, top and bottom
image = imcrop(image, [200 1 800 1024]);
topHalf = image(3:250,:);
botHalf = image(251:end, :); 

%% Top half calculation
r = 50;
[U, S, V] = svd(cast(topHalf, 'double'));
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
toprecon =Ur*Sr*Vr';

top= medfilt2(toprecon, [3 3]);

% Can adjust binarize threshold 
top = imbinarize(top,30);


%% Bottom half calculation
r = 50;
[U, S, V] = svd(cast(botHalf, 'double'));
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);
botrecon =Ur*Sr*Vr';

bot = medfilt2(botrecon, [3 3]);
% Can adjust binarize threshold 
bot = imbinarize(bot,70);

% Stitch together top and bottom half
imageRecon = [top;bot];

output = imageRecon;

end


% Example of usage:
% image =  imread('temp701.png');
% image2 = denoiser(image);
% 
% figure()
% subplot(2,1,1)
% title('original')
% imshow(image)
% 
% subplot(2,1,2)
% title('Denoised')
% imshow(image2)





