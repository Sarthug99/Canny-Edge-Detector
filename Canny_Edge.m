clc;
clear all;
close all;
I = imread('coins.png');
[m,n] = size(I);
%% GAUSSIAN FILTERING
sd = input("Enter standard deviation for kernel");
dim = input("Enter size of kernel");
a = (dim+1)/2;
Gaussian = zeros(dim,dim);
for x = 1:dim
    for y= 1:dim
        Gaussian(x,y) = (exp(((x-a)^2 + (y-a)^2)/(2*sd*sd)))^-1;    
    end
end
r = Gaussian(1,1);
Gaussian = 1/r.*Gaussian; 
G = round(Gaussian);
t =  sum(sum(G)); 
zero_padded = zeros(m+a,n+a);
for i= 1:m+a
    for j= 1:n+a
        if(i<m+1 && j<n+1)
            zero_padded(i+1,j+1) = I(i,j);
        end
    end
end
I_lpf = zeros(m,n);
for k = a:m+1
    for l= a:n+1
        for r = -(a-1):a-1
            for s = -(a-1):a-1 
                I_lpf(k-1,l-1) = I_lpf(k-1,l-1) + (zero_padded(k+r, l+s)*G(r+a,s+a))/t ;
            end
        end
    end
end
subplot 221
imshow(uint8(I_lpf));

for i= 1:m+a
    for j= 1:n+a
        if(i<m+1 && j<n+1)
            zero_padded(i+1,j+1) = I_lpf(i,j);
        end
    end
end

%% HGH PASS FILTERING
Gy = [-1 -2 -1; 0 0 0 ; 1 2 1];
Gx = [-1 0 1; -2 0 2; -1 0 1];
Mag = zeros(m,n);
Phase = zeros(m,n);
for k = 2:m+1
    for l= 2:n+1
        gradx = 0;
        grady = 0;
        for r = -1:1
            for s = -1:1 
                  gradx = gradx + zero_padded(k+r, l+s)*Gx(r+2,s+2) ;
                  grady = grady + zero_padded(k+r, l+s)*Gy(r+2,s+2) ;
            end
        end
        Mag(k-1,l-1) = ((gradx)^2 + (grady)^2)^(0.5);
        Phase(k-1,l-1) = atan(grady/gradx);
    end
end
subplot 222
imshow(uint8(Mag));
Phase = Phase + (pi/2);
Phase = Phase*4/pi;
Phase = floor(Phase);


%% QUANTIZATION AND SUPPRESSION OF NON MAXIMA
 for i = 2 : m-1
     for j = 2 : n-1
         if(Phase == 0)
             a = Mag(i-1,j);
             b = Mag(i,j);
             c = Mag(i+1,j);
         elseif(Phase == 1)
             a = Mag(i+1,j-1);
             b = Mag(i,j);
             c = Mag(i-1,j+1);
         elseif(Phase == 2)
             a = Mag(i,j-1);
             b = Mag(i,j);
             c = Mag(i,j+1);
         else
             a = Mag(i-1,j-1);
             b = Mag(i,j);
             c = Mag(i+1,j+1);
         end
          if( b < a | b < c )
              Mag(i,j) = 0;
          end
     end
 end
 Mag(:,1)=0;
 Mag(1,:)= 0 ;
 Mag(m,:)=0;
 Mag(:,n)=0;
 subplot 223
imshow(uint8(Mag)); 
 

 %% Thresholding
Thigh = 0.7*max(max(Mag));
Tlow = 0.105*max(max(Mag));
Normalized = Mag;
for i = 1:m
    for j = 1:n
        if(Mag(i,j) > Thigh)
            Normalized(i,j)=1;
        elseif(Normalized(i,j) < Tlow)
            Normalized(i,j)=0;
        end
    end
end
subplot 224
imshow(Normalized);
