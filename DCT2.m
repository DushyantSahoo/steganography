% ideal DCT with no change
% embed in lower frequency 
% less huffman length 
% ideal quantisation table

clear all; 
close all;

% p and q are 2 prime numbers 
% phi is totient of the product of p and q
% e is any number coprime to phi
% d, the modular multiplicative inverse of e (mod φ(n))

p = input('\nEnter value of p ( It should be prime number eg 31): ');
q = input('\nEnter value of q ( It should be prime number eg 19): ');
[Pk,Phi,d,e] = init(p,q);
M = input('\nEnter message: ','s');
x=length(M);
c=0;

for j= 1:x
    for i=0:122
        if strcmp(M(j),char(i))
            c(j)=i;
        end
    end
end

disp('ASCII Code of the entered Message:');
disp(c); 

% For Encryption

for j= 1:x
   cipher(j)= crypto(c(j),Pk,e); 
end

disp('Cipher Text of the entered Message:');
disp(cipher);

tempcipher=dec2bin(cipher);

% determine the size of message object to embed
Mm=size(tempcipher,1);	        %Height
Nm=size(tempcipher,2);	        %Width

% Compression

jk=8;
block_size=jk;
mask_size=jk;
    
I1 = imread('1.gif');
imwrite(I1,'32.jpg');
I1=im2double(I1); 
I1=I1*255;    
no1=(floor(size(I1,1)/(block_size)))*block_size;
no2=(floor(size(I1,2)/block_size))*block_size;
I1=imresize(I1,[no1,no2]);

    
I=I1;
Red= I(:,:,1);
T = dctmtx(block_size);
dct = @(block_struct) T* block_struct.data*T';
  
B1 = blockproc(I,[8,8],dct);
   
transformed=cat(1,B1);
    
%% table generation
quantization_table = [
 16  11  10  16  24   40   51   61;
 12  12  14  19  26   58   60   55;
 14  13  16  24  40   57   69   56;
 14  17  22  29  51   87   80   62;
 18  22  37  56  68   109  103  77;
 24  35  55  64  81   104  113  92;
 49  64  78  87  103  121  120  101;
 72  92  95  98  112  100  103  99
    ];

%quantization_table = ones(8,8);
%quantization_table(1,1) = 16;quantization_table(1,2) = 11;quantization_table(1,3) = 10;quantization_table(1,4) = 16;
%quantization_table(2,1) = 12;quantization_table(2,2) = 12;quantization_table(2,3) = 14;
%quantization_table(3,1) = 14;quantization_table(3,2) = 13;
%quantization_table(4,1) = 14;
    
mask=zeros(block_size,block_size);
for i=1:mask_size
    mask(i,1:(mask_size-i+1))=1;
end
mask;
    
%% quantization
quant=int16(zeros(size(I1,1),size(I1,1),1));
recon=double(zeros(size(I1,1),size(I1,1),1));
for k=1:1
    for i=1:block_size:size(I1,1)
        for j=1:block_size:size(I1,2)
            for ii=1:block_size
                for jj=1:block_size
                    aa=transformed(i+ii-1,j+jj-1,k);
                    quant(i+ii-1,j+jj-1,k)=(aa);
                    quant(i+ii-1,j+jj-1,k)=(aa/quantization_table(ii,jj));
                    
                end
            end
        end
    end
end

keymat=find(quant>10,Mm*Nm);

for ii = 1:Mm
    for jj = 1:Nm
       quant(keymat(Nm*(ii-1)+jj))=bitset(quant(keymat(Nm*(ii-1)+jj)),1,str2num(tempcipher(ii,jj)));
    end
end

% zigzag coding

outputas=[];
for i=1:block_size:size(I1,1)
    for j=1:block_size:size(I1,2)
        outputas=cat(2,outputas,zigzag(quant(i:i+7,j:j+7)));
    end
end

% hufman coding 
tb = tabulate(outputas);
pval = tb(:,3);
symbols=tb(:,1);
[dict,avglen] = huffmandict(symbols,(pval/100));
comp = huffmanenco(outputas,dict);

% hufman decoding
outputas = huffmandeco(comp,dict);


% inverse zigzag
tempout=[size(I1,1),size(I1,2)];
tempi=1;
tempj=1;
for i=1:64:size(outputas,2)
    tempout(tempi:tempi+7,tempj:tempj+7)=izigzag(outputas(i:i+63),8,8);
    tempj=tempj+8;
    if tempj > size(I1,2)
        tempi =tempi+8;
        tempj =1;
    end
    
end

quant=tempout;

%% dequantisation
for k=1:1
    for i=1:block_size:size(I1,1)
        for j=1:block_size:size(I1,2)
            for ii=1:block_size
                for jj=1:block_size
                    
                    recon(i+ii-1,j+jj-1,k)=(quant(i+ii-1,j+jj-1,k)*quantization_table(ii,jj));
                end
            end
        end
    end
end
    

%% reconstruction    
B1=(recon(:,:,1));
invdct = @(block_struct) T' * block_struct.data * T;
RE1 = blockproc(B1,[block_size block_size],invdct);

imwrite(RE1,'341.jpg');

%% decipher
re1 = blockproc(RE1,[8,8],dct);
quant1=int16(zeros(size(I1,1),size(I1,1),1));
for k=1:1
    for i=1:block_size:size(I1,1)
        for j=1:block_size:size(I1,2)
            for ii=1:block_size
                for jj=1:block_size
                    aa=re1(i+ii-1,j+jj-1,k);
                    quant1(i+ii-1,j+jj-1,k)=(aa);
                    quant1(i+ii-1,j+jj-1,k)=(aa)/quantization_table(ii,jj);
                end
            end
        end
    end
end

for ii = 1:Mm
    for jj = 1:Nm
        temp1cipher(ii,jj)=bitget(quant1(keymat(Nm*(ii-1)+jj)),1);
    end
end

I2=cat(1,RE1);
    
%% error
I1 = I1/255;
I2 = I2/255;
result1=0;
for i=1:size(I1,1)
    for j=1:size(I1,2)
        for k=size(I1,3)
            diff1=(I1(i,j,k)-I2(i,j,k));
            result1=result1+diff1*diff1;                
        end
     end
end
  
sizee=size(I2,1);    
MSE1 = result1/sizee;
PSNR1 = 10*log10(255*255/MSE1);    
PSNR1

temp1cipher=num2str(temp1cipher);
cipher = bin2dec(temp1cipher);
for j= 1:x
   message(j)= crypto(cipher(j),Pk,d); 
end

disp('Decrypted ASCII of Message:');
disp(message);
disp(['Decrypted Message is: ' message]);
    
    
