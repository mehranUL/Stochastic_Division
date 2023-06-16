clear
tic
bp = 8;
N = 2^(bp);
% X1 = randi((N+1),1)-1; %pseudorandom integer from a uniform discrete distribution
%     %range -> [0,N]
% 
% X2 = randi((N+1),1)-1;
% X1 = 5;
% X2 = 105;

X1 = 56; %56
X2 = 198; %198

sobol_seq = net(sobolset(100), N);
sobol_seq_new = sobol_seq(:,1);

vd(:,1) = vdcorput(N-1,8);    %164,163,162,161,119,118,109,108,107,106
%vd(:,2) = vdcorput(N-1,3);   % 191,181,199,203

R2_ws = load("R2_1024.mat");
R2 = R2_ws.z;

Y = lhsdesign(N,20);
Y = Y(:,13:14);

alpha = (sqrt(2) - 1); % Silver ratio
weyl(:,1) = mod((1:N)*alpha, 1);
beta = pi;
weyl(:,2) = mod((1:N)*beta, 1);

ff = faure(N,2,7);
ff = ff';
ff = ff(1:N,:);

HT = net(haltonset(6),N);
HT = HT(:,5:6); %Bases 11,13

HH = Hammersley(N,3);
HH = HH';
HH = HH(:,2:3);

[ nd1, ~ ] = niederreiter2_generate ( 20, N, 2, 31 );
nd1 = nd1';
%nd = nd1(:,1:4);
nd(:,1) = nd1(:,1);
nd(:,2) = nd1(:,2);

pt = poissonDisc([100,100,100],7,N); %6 for 2k-- 5 for 4k --- 4 for 8k
pt(:,1) = pt(:,1)/100;
pt(:,2) = pt(:,2)/100;

% LF = LFSR__2([true false true true false false false],N);
% LF(2,:) = LFSR__2([true false true true false true false],N);

%LF = LFSR__2([true false true true true false true false],N);
%LF(2,:) = LFSR__2([true false true true false true true false],N);
%LF(1,:) = LF(1,:)/N;
%LF(2,:) = LF(2,:)/N;
%LF = LF';
%LF = rand(N,1);

X1_stream_sobol = zeros(2^bp, N);
X2_stream_sobol = zeros(2^bp, N);
X1_stream_vandercorput = zeros(2^bp, N);
X2_stream_vandercorput = zeros(2^bp, N);
X1_stream_R2 = zeros(2^bp, N);
X2_stream_R2 = zeros(2^bp, N);
% X1_stream_weyl = zeros(1, N);
% X2_stream_weyl = zeros(1, N);
% X1_stream_latin = zeros(1, N);
% X2_stream_latin = zeros(1, N);
% X1_stream_faure = zeros(1, N);
% X2_stream_faure = zeros(1, N);
% X1_stream_halton = zeros(1, N);
% X2_stream_halton = zeros(1, N);
% X1_stream_hammersley = zeros(1, N);
% X2_stream_hammersley = zeros(1, N);
% X1_stream_nd = zeros(1, N);
% X2_stream_nd = zeros(1, N);
% X1_stream_ps = zeros(1, N);
% X2_stream_ps = zeros(1, N);
X1_stream_LF = zeros(2^bp, N);
X2_stream_LF = zeros(2^bp, N);

J_stream_sobol = zeros(2^bp,2^bp,N);
K_stream_sobol = zeros(2^bp,2^bp,N);
J_stream_vd = zeros(2^bp,2^bp,N);
K_stream_vd = zeros(2^bp,2^bp,N);
J1_sobol = zeros(size(J_stream_sobol));

J_stream_LF = zeros(2^bp,2^bp,N);
K_stream_LF = zeros(2^bp,2^bp,N);

J_stream_R2 = zeros(2^bp,2^bp,N);
K_stream_R2 = zeros(2^bp,2^bp,N);

J_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
K_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
J_SS_vd = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
K_SS_vd = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
J_SS_LF = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
K_SS_LF = zeros(2^bp,2^bp,N);

Z_stream_sobol = ones(2^bp,2^bp,N);
Z_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
Z_CORDIV_sobol = zeros(2^bp,2^bp,N);

Z_stream_LF = zeros(2^bp,2^bp,N);
Z_CORDIV_LF = zeros(2^bp,2^bp,N);
Z_SS_LF = zeros(2^bp,2^bp,N);

Z_stream_R2 = ones(2^bp,2^bp,N);

Z_stream_vd = zeros(2^bp,2^bp,N);
Z_SS_vd = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
Z_CORDIV_vd = zeros(2^bp,2^bp,N);
% Z_stream_R2 = zeros(1,N);
% Z_stream_weyl = zeros(1,N);
% Z_stream_latin = zeros(1,N);
% Z_stream_faure = zeros(1,N);
% Z_stream_halton = zeros(1,N);
% Z_stream_hammersley = zeros(1,N);
% Z_stream_nd = zeros(1,N);
% Z_stream_ps = zeros(1,N);
% Z_stream_LF = zeros(1,N);

X_CORLD = zeros(2^bp,N);
X_CORLD_vd = zeros(2^bp,N);


Expected_div = zeros(2^bp,2^bp);

Err_sobol = zeros(2^bp,2^bp);
Err_sobol_sq = zeros(2^bp,2^bp);
Err_sobol_SS = zeros(2^bp,2^bp);
Err_sobol_CORDIV = zeros(2^bp,2^bp);

Err_LF = zeros(2^bp,2^bp);
Err_LF_SS = zeros(2^bp,2^bp);
Err_LF_CORDIV = zeros(2^bp,2^bp);

Err_R2 = zeros(2^bp,2^bp);

Err_vd = zeros(2^bp,2^bp);
Err_vd_sq = zeros(2^bp,2^bp);
Err_vd_SS = zeros(2^bp,2^bp);
Err_vd_CORDIV = zeros(2^bp,2^bp);


for i = 1:2^bp
    for j = 1:2^bp
        if i < j
            m = i;
            s = i;
            for z=1:N
%                 if (((i-1)/2^bp) > sobol_seq_new(z,1))   %(i-1)/N
%                     X1_stream_sobol(i,z) = 1;
%                     
%                 end
                

                if (((j-1)/2^bp) > sobol_seq_new(z,1))
                    X2_stream_sobol(j,z) = 1;                    
                end
%                 if (((i-1)/2^bp) > vd(z,1))
%                     X1_stream_vandercorput(i,z) = 1;
%                 end
                if (((j-1)/2^bp) > vd(z,1))
                    X2_stream_vandercorput(j,z) = 1;
                end
                if (((i-1)/2^bp) > R2(z,2))
                    X1_stream_R2(i,z) = 1;
                end
                if (((j-1)/2^bp) > R2(z,2))
                    X2_stream_R2(j,z) = 1;
                end
%         if ((X1/N) > weyl(z,1))
%             X1_stream_weyl(1,z) = 1;
%         end
%         if ((X2/N) > weyl(z,1))
%             X2_stream_weyl(1,z) = 1;
%         end
%         if ((X1/N) > Y(z,1))
%             X1_stream_latin(1,z) = 1;
%         end
%         if ((X2/N) > Y(z,1))
%             X2_stream_latin(1,z) = 1;
%         end
%         if ((X1/N) > ff(z,1))
%             X1_stream_faure(1,z) = 1;
%         end
%         if ((X2/N) > ff(z,1))
%             X2_stream_faure(1,z) = 1;
%         end
%         if ((X1/N) > HT(z,1))
%             X1_stream_halton(1,z) = 1;
%         end
%         if ((X2/N) > HT(z,1))
%             X2_stream_halton(1,z) = 1;
%         end
%         if ((X1/N) > HH(z,1))
%             X1_stream_hammersley(1,z) = 1;
%         end
%         if ((X2/N) > HH(z,1))
%             X2_stream_hammersley(1,z) = 1;
%         end
%         if ((X1/N) > nd(z,1))
%             X1_stream_nd(1,z) = 1;
%         end
%         if ((X2/N) > nd(z,1))
%             X2_stream_nd(1,z) = 1;
%         end
%         if ((X1/N) > pt(z,1))
%             X1_stream_ps(1,z) = 1;
%         end
%         if ((X2/N) > pt(z,1))
%             X2_stream_ps(1,z) = 1;
%         end
%                 if (((i-1)/2^bp) > LF(z,1))
%                     X1_stream_LF(i,z) = 1;
%                 end
%                 if (((j-1)/2^bp) > LF(z,1))
%                     X2_stream_LF(j,z) = 1;
%                 end

            end
            for d = 1:N
              if m>0
                X_CORLD(i,d) = X2_stream_sobol(j,d);
                if X2_stream_sobol(j,d) == 1
                    m = m-1;
                end
              else
                  X_CORLD(i,d) = 0;
              end
         
            end
            for d = 1:N
              if s>0
                X_CORLD_vd(i,d) = X2_stream_vandercorput(j,d);
                if X2_stream_vandercorput(j,d) == 1
                    s = s-1;
                end
              else
                  X_CORLD_vd(i,d) = 0;
              end
         
            end



%           Saturated Subtractor Divider Design
%             J_SS_sobol(i,j,:) = X1_stream_sobol(i,:);
%             K_SS_sobol(i,j,:) = and(not(X1_stream_sobol(i,:)),X2_stream_sobol(j,:));
            J_SS_sobol(i,j,:) = X_CORLD(i,:);
            K_SS_sobol(i,j,:) = and(not(X_CORLD(i,:)),X2_stream_sobol(j,:));

% 
            %J_SS_vd(i,j,:) = X1_stream_vandercorput(i,:);
            %K_SS_vd(i,j,:) = and(not(X1_stream_vandercorput(i,:)),X2_stream_vandercorput(j,:));
            %K_SS_vd(i,j,:) = xor(X1_stream_vandercorput(i,:),X2_stream_vandercorput(j,:));
            J_SS_vd(i,j,:) = X_CORLD_vd(i,:);
            K_SS_vd(i,j,:) = and(not(X_CORLD_vd(i,:)),X2_stream_vandercorput(j,:));

% 
% 
%             J_stream_sobol(i,j,:) = and(X1_stream_sobol(i,:),X2_stream_sobol(j,:));
%             K_stream_sobol(i,j,:) = xor(X1_stream_sobol(i,:),X2_stream_sobol(j,:));
            J_stream_sobol(i,j,:) = and(X_CORLD(i,:),X2_stream_sobol(j,:));
            K_stream_sobol(i,j,:) = xor(X_CORLD(i,:),X2_stream_sobol(j,:));

% 
%             J_stream_vd(i,j,:) = and(X1_stream_vandercorput(i,:),X2_stream_vandercorput(j,:));
%             K_stream_vd(i,j,:) = xor(X1_stream_vandercorput(i,:),X2_stream_vandercorput(j,:));
            J_stream_vd(i,j,:) = and(X_CORLD_vd(i,:),X2_stream_vandercorput(j,:));
            K_stream_vd(i,j,:) = xor(X_CORLD_vd(i,:),X2_stream_vandercorput(j,:));

            J_stream_R2(i,j,:) = and(X1_stream_R2(i,:),X2_stream_R2(j,:));
            K_stream_R2(i,j,:) = xor(X1_stream_R2(i,:),X2_stream_R2(j,:));

%             J_stream_LF(i,j,:) = and(X1_stream_LF(i,:),X2_stream_LF(j,:));
%             K_stream_LF(i,j,:) = xor(X1_stream_LF(i,:),X2_stream_LF(j,:));
% 
%             J_SS_LF(i,j,:) = X1_stream_LF(i,:);
%             K_SS_LF(i,j,:) = and(not(X1_stream_LF(i,:)),X2_stream_LF(j,:));

% J_stream_vd = X1_stream_vandercorput;
% K_stream_vd = and(not(X1_stream_vandercorput),X2_stream_vandercorput);
% 
% J_stream_R2 = X1_stream_R2;
% K_stream_R2 = and(not(X1_stream_R2),X2_stream_R2);
% 
% J_stream_weyl = X1_stream_weyl;
% K_stream_weyl = and(not(X1_stream_weyl),X2_stream_weyl);
% 
% J_stream_faure = X1_stream_faure;
% K_stream_faure = and(not(X1_stream_faure),X2_stream_faure);
% 
% J_stream_latin = X1_stream_latin;
% K_stream_latin = and(not(X1_stream_latin),X2_stream_latin);
% 
% J_stream_halton = X1_stream_halton;
% K_stream_halton = and(not(X1_stream_halton),X2_stream_halton);
% 
% J_stream_hammersley = X1_stream_hammersley;
% K_stream_hammersley = and(not(X1_stream_hammersley),X2_stream_hammersley);
% 
% J_stream_nd = X1_stream_nd;
% K_stream_nd = and(not(X1_stream_nd),X2_stream_nd);
% 
% J_stream_ps = X1_stream_ps;
% K_stream_ps = and(not(X1_stream_ps),X2_stream_ps);
% 
% J_stream_LF = X1_stream_LF;
% K_stream_LF = and(not(X1_stream_LF),X2_stream_LF);

%temp = reshape(J_stream_sobol(i,j,:),1,256);
%J1_sobol(i,j,3:end) = J_stream_sobol(i,j,1:end-2);
%J1_sobol(i,j,:) = circshift(J_stream_sobol(i,j,:),7);
            
            
            %J_stream_sobol(i,j,:) = circshift(J_stream_sobol(i,j,:),1);
%             Z_stream_sobol(i,j,:) = JK_behaviour(J_stream_sobol(i,j,:),K_stream_sobol(i,j,:),N);
%             Z_SS_sobol(i,j,:) = JK_behaviour(J_SS_sobol(i,j,:),K_SS_sobol(i,j,:),N);
            %Z_CORDIV_sobol(i,j,:) = or(and(not(X2_stream_sobol(j,:)),reshape(circshift(Z_CORDIV_sobol(i,j,:),1),1,N)),and(X2_stream_sobol(j,:),X1_stream_sobol(i,:)));

            %J_stream_LF(i,j,:) = circshift(J_stream_LF(i,j,:),1);
%             Z_stream_LF(i,j,:) = JK_behaviour(J_stream_LF(i,j,:),K_stream_LF(i,j,:),N);
%             Z_SS_LF(i,j,:) = JK_behaviour(J_SS_LF(i,j,:),K_SS_LF(i,j,:),N);

            %J_stream_vd(i,j,:) = circshift(J_stream_vd(i,j,:),1);
            Z_stream_vd(i,j,:) = JK_behaviour(J_stream_vd(i,j,:),K_stream_vd(i,j,:),N);
            Z_SS_vd(i,j,:) = JK_behaviour(J_SS_vd(i,j,:),K_SS_vd(i,j,:),N);

            Z_stream_R2(i,j,:) = JK_behaviour(J_stream_R2(i,j,:),K_stream_R2(i,j,:),N);

            for k = 1:N
%                 J_temp(i,j,1:end-2) = J_stream_sobol(i,j,3:end);
%                 Z_stream_sobol(i,j,k) = JK_FF(J_temp(i,j,k),K_stream_sobol(i,j,k));
                 %Z_SS_sobol(i,j,k) = JK_FF(J_SS_sobol(i,j,k),K_SS_sobol(i,j,k));
               % Z_CORDIV_sobol(i,j,k+1) = bitor(bitand(not(X2_stream_sobol(j,k)),Z_CORDIV_sobol(i,j,k)),and(X2_stream_sobol(j,k),X1_stream_sobol(i,k)));
                Z_CORDIV_sobol(i,j,k+1) = bitor(bitand(not(X2_stream_sobol(j,k)),Z_CORDIV_sobol(i,j,k)),and(X2_stream_sobol(j,k),X_CORLD(i,k)));
               

                %Z_CORDIV_LF(i,j,k+1) = bitor(bitand(not(X2_stream_LF(j,k)),Z_CORDIV_LF(i,j,k)),and(X2_stream_LF(j,k),X1_stream_LF(i,k)));

                %Z_stream_sobol(i,j,k) = JK_FF(J1_sobol(i,j,k),K_stream_sobol(i,j,k));
%                 Z_stream_vd(i,j,k) = JK_FF(J_stream_vd(i,j,k),K_stream_vd(i,j,k));
%                 Z_SS_vd(i,j,k) = JK_FF(J_SS_vd(i,j,k),K_SS_vd(i,j,k));
                %Z_CORDIV_vd(i,j,k+1) = bitor(bitand(not(X2_stream_vandercorput(j,k)),Z_CORDIV_vd(i,j,k)),and(X2_stream_vandercorput(j,k),X1_stream_vandercorput(i,k)));
                Z_CORDIV_vd(i,j,k+1) = bitor(bitand(not(X2_stream_vandercorput(j,k)),Z_CORDIV_vd(i,j,k)),and(X2_stream_vandercorput(j,k),X_CORLD_vd(i,k)));
%     Z_stream_R2(k) = JK_FF(J_stream_R2(k),K_stream_R2(k));
%     Z_stream_weyl(k) = JK_FF(J_stream_weyl(k),K_stream_weyl(k));
%     Z_stream_faure(k) = JK_FF(K_stream_faure(k),K_stream_faure(k));
%     Z_stream_latin(k) = JK_FF(J_stream_latin(k),K_stream_latin(k));
%     Z_stream_halton(k) = JK_FF(J_stream_halton(k),K_stream_halton(k));
%     Z_stream_hammersley(k) = JK_FF(J_stream_hammersley(k),K_stream_hammersley(k));
%     Z_stream_nd(k) = JK_FF(J_stream_nd(k),K_stream_nd(k));
%     Z_stream_ps(k) = JK_FF(J_stream_ps(k),K_stream_ps(k));
%     Z_stream_LF(k) = JK_FF(J_stream_LF(k),K_stream_LF(k));
            end
            

% if X1 > X2
    %if i > j
%     temp = X2;
%     X2 = X1;
%     X1 = temp;
%     Expexted_div = X1/X2
        %Expexted_div(i,j) = j/i;
    %else
%     Expexted_div = X1/X2
        %if j == 1
         %   Expexted_div(i,j) = 0;
        %else
            Expected_div(i,j) = (i-1)/(j-1);    %(i-1)/(j-1)
        %end
    %end
            Err_sobol(i,j) = abs((sum(Z_stream_sobol(i,j,:))/N) - Expected_div(i,j));
%             Err_sobol_sq(i,j) = power((sum(Z_stream_sobol(i,j,:))/N) - Expected_div(i,j),2);
            Err_sobol_SS(i,j) = abs((sum(Z_SS_sobol(i,j,:))/N) - Expected_div(i,j));
            Err_sobol_CORDIV(i,j) = abs((sum(Z_CORDIV_sobol(i,j,1:N))/N) - Expected_div(i,j));

%             Err_LF(i,j) = abs((sum(Z_stream_LF(i,j,:))/N) - Expected_div(i,j));
%             Err_LF_SS(i,j) = abs((sum(Z_SS_LF(i,j,:))/N) - Expected_div(i,j));
%             Err_LF_CORDIV(i,j) = abs((sum(Z_CORDIV_LF(i,j,1:N))/N) - Expected_div(i,j));
            Err_R2(i,j) = abs((sum(Z_stream_R2(i,j,:))/N) - Expected_div(i,j));

            Err_vd(i,j) = abs((sum(Z_stream_vd(i,j,:))/N) - Expected_div(i,j));
%             Err_vd_sq(i,j) = power((sum(Z_stream_vd(i,j,:))/N) - Expected_div(i,j),2);
            Err_vd_SS(i,j) = abs((sum(Z_SS_vd(i,j,:))/N) - Expected_div(i,j));
            Err_vd_CORDIV(i,j) = abs((sum(Z_CORDIV_vd(i,j,1:N))/N) - Expected_div(i,j));
% Err_R2 = abs((sum(Z_stream_R2)/N) - Expexted_div)
% Err_weyl = abs((sum(Z_stream_weyl)/N) - Expexted_div)
% Err_faure = abs((sum(Z_stream_faure)/N) - Expexted_div)
% Err_latin = abs((sum(Z_stream_latin)/N) - Expexted_div)
% Err_halton = abs((sum(Z_stream_halton)/N) - Expexted_div)
% Err_hammersley = abs((sum(Z_stream_hammersley)/N) - Expexted_div)
% Err_nd = abs((sum(Z_stream_nd)/N) - Expexted_div)
% Err_ps = abs((sum(Z_stream_ps)/N) - Expexted_div)
% Err_LF = abs((sum(Z_stream_LF)/N) - Expexted_div)
        end
    end
end

% MAE_S = sum(Err_sobol,"all")/((256)*(256))
% MAE_S_SS = sum(Err_sobol_SS,"all")/((256)*(256))

%MAE_Sobol_min_max = sum(Err_sobol,"all")/numel(find(Err_sobol ~= 0))
MAE_Sobol_min_max = sum(Err_sobol,"all")/(128*255)      % (2^(bp-1)*((2^bp)-1))
% MSE_Sobol = sum(Err_sobol_sq,"all")/numel(find(Err_sobol_sq ~= 0))
MAE_Sobol_Saturated_Sub = sum(Err_sobol_SS,"all")/(128*255)
MAE_Sobol_CORDIV = sum(Err_sobol_CORDIV,"all")/(128*255)

%MAE_VD = sum(Err_vd,"all")/((128)*(255))

MAE_VD_min_max = sum(Err_vd,"all")/(128*255)
% MSE_VD = sum(Err_vd_sq,"all")/numel(find(Err_vd_sq ~= 0))

%MAE_VD_SS = sum(Err_vd_SS,"all")/((128)*(255))

MAE_VD_Saturated_Sub = sum(Err_vd_SS,"all")/(128*255)
MAE_VD_CORDIV = sum(Err_vd_CORDIV,"all")/(128*255)

% MAE_LFSR_min_max = sum(Err_LF,"all")/(128*255)
% MAE_LFSR_Saturated_Sub = sum(Err_LF,"all")/(128*255)
% MAE_LFSR_CORDIV = sum(Err_LF,"all")/(128*255)

MAE_R2_min_max = sum(Err_R2,"all")/(128*255)
toc
