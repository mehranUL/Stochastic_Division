% Shared RNG Design Space Exploration (MAE)
clear
tic
bp = 8;
N = 2^(bp);



sobol_seq = net(sobolset(100), N);
sobol_seq_new = sobol_seq(:,5);


vd(:,1) = vdcorput(N-1,16);    %164,163,162,161,119,118,109,108,107,106
%vd(:,2) = vdcorput(N-1,3);   % 191,181,199,203



% LF = LFSR__2([true false true true false false false],N);
% LF(2,:) = LFSR__2([true false true true false true false],N);

% LF = LFSR__2([true false true true true false false false],N);
% % LF(2,:) = LFSR__2([true false true true false true true false],N);
% LF(1,:) = LF(1,:)/N;
% % LF(2,:) = LF(2,:)/N;
% LF = LF';
load('LF_256.mat');

X1_stream_sobol = zeros(2^bp, N);
X2_stream_sobol = zeros(2^bp, N);
X1_stream_vandercorput = zeros(256, N);
X2_stream_vandercorput = zeros(256, N);
% 
X1_stream_LF = zeros(256, N);
X2_stream_LF = zeros(256, N);


J_stream_sobol = zeros(2^bp,2^bp,N);
%J_temp = zeros(256,256,N);
K_stream_sobol = zeros(2^bp,2^bp,N);
J_stream_vd = zeros(256,256,N);
K_stream_vd = zeros(256,256,N);
%J1_sobol = zeros(size(J_stream_sobol));

J_stream_LF = zeros(2^bp,2^bp,N);
K_stream_LF = zeros(2^bp,2^bp,N);

J_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
K_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
J_SS_vd = zeros(256,256,N);  %Saturated Subtractor design
K_SS_vd = zeros(256,256,N);  %Saturated Subtractor design

J_SS_LF = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
K_SS_LF = zeros(2^bp,2^bp,N);

Z_stream_sobol = zeros(2^bp,2^bp,N);
Z_SS_sobol = zeros(2^bp,2^bp,N);  %Saturated Subtractor design
Z_CORDIV_sobol = zeros(2^bp,2^bp,N);

Z_stream_vd = zeros(256,256,N);
Z_SS_vd = zeros(256,256,N);  %Saturated Subtractor design
Z_CORDIV_vd = zeros(256,256,N);

Z_stream_LF = zeros(2^bp,2^bp,N);
Z_CORDIV_LF = zeros(2^bp,2^bp,N);
Z_SS_LF = zeros(2^bp,2^bp,N);

Expected_div = zeros(2^bp,2^bp);

Err_sobol = zeros(2^bp,2^bp);
%Err_sobol_sq = zeros(2^bp,2^bp);
Err_sobol_SS = zeros(2^bp,2^bp);
Err_sobol_CORDIV = zeros(2^bp,2^bp);

Err_vd = zeros(256,256);
% Err_vd_sq = zeros(256,256);
Err_vd_SS = zeros(256,256);
Err_vd_CORDIV = zeros(256,256);

Err_LF = zeros(2^bp,2^bp);
Err_LF_SS = zeros(2^bp,2^bp);
Err_LF_CORDIV = zeros(2^bp,2^bp);


for i = 1:2^bp
    for j = 1:2^bp
        if i < j
            for z=1:N
                if (((i-1)/N) > sobol_seq_new(z,1))
                    X1_stream_sobol(i,z) = 1;
                end
                if (((j-1)/N) > sobol_seq_new(z,1))
                    X2_stream_sobol(j,z) = 1;
                end
                if (((i-1)/N) > vd(z,1))
                    X1_stream_vandercorput(i,z) = 1;
                end
                if (((j-1)/N) > vd(z,1))
                    X2_stream_vandercorput(j,z) = 1;
                end
                if (((i-1)/2^bp) > LF(z,1))
                    X1_stream_LF(i,z) = 1;
                end
                if (((j-1)/2^bp) > LF(z,1))
                    X2_stream_LF(j,z) = 1;
                end


            end

%             Saturated Subtractor Divider Design
            J_SS_sobol(i,j,:) = X1_stream_sobol(i,:);
            K_SS_sobol(i,j,:) = and(not(X1_stream_sobol(i,:)),X2_stream_sobol(j,:));

            
            J_stream_sobol(i,j,:) = and(X1_stream_sobol(i,:),X2_stream_sobol(j,:));
            K_stream_sobol(i,j,:) = xor(X1_stream_sobol(i,:),X2_stream_sobol(j,:));

            J_SS_vd(i,j,:) = X1_stream_vandercorput(i,:);
            K_SS_vd(i,j,:) = and(not(X1_stream_vandercorput(i,:)),X2_stream_vandercorput(j,:));

            J_stream_vd(i,j,:) = and(X1_stream_vandercorput(i,:),X2_stream_vandercorput(j,:));
            K_stream_vd(i,j,:) = xor(X1_stream_vandercorput(i,:),X2_stream_vandercorput(j,:));

            %J_stream_LF(i,j,:) = circshift(J_stream_LF(i,j,:),2);
            J_stream_LF(i,j,:) = and(X1_stream_LF(i,:),X2_stream_LF(j,:));
            K_stream_LF(i,j,:) = xor(X1_stream_LF(i,:),X2_stream_LF(j,:));

            J_SS_LF(i,j,:) = X1_stream_LF(i,:);
            K_SS_LF(i,j,:) = and(not(X1_stream_LF(i,:)),X2_stream_LF(j,:));

            %J_stream_sobol(i,j,:) = circshift(J_stream_sobol(i,j,:),1);
            Z_stream_sobol(i,j,:) = JK_behaviour(J_stream_sobol(i,j,:),K_stream_sobol(i,j,:),N);
            Z_SS_sobol(i,j,:) = JK_behaviour(J_SS_sobol(i,j,:),K_SS_sobol(i,j,:),N);

            %J_stream_vd(i,j,:) = circshift(J_stream_vd(i,j,:),1);
            Z_stream_vd(i,j,:) = JK_behaviour(J_stream_vd(i,j,:),K_stream_vd(i,j,:),N);
            Z_SS_vd(i,j,:) = JK_behaviour(J_SS_vd(i,j,:),K_SS_vd(i,j,:),N);

            %J_stream_LF(i,j,:) = circshift(J_stream_LF(i,j,:),1);
            Z_stream_LF(i,j,:) = JK_behaviour(J_stream_LF(i,j,:),K_stream_LF(i,j,:),N);
            Z_SS_LF(i,j,:) = JK_behaviour(J_SS_LF(i,j,:),K_SS_LF(i,j,:),N);


            for k = 1:N
%                 J_temp(i,j,1:end-2) = J_stream_sobol(i,j,3:end);
%                 Z_stream_sobol(i,j,k) = JK_FF(J_temp(i,j,k),K_stream_sobol(i,j,k));
%                 Z_SS_sobol(i,j,k) = JK_FF(J_SS_sobol(i,j,k),K_SS_sobol(i,j,k));
                Z_CORDIV_sobol(i,j,k+1) = bitor(bitand(not(X2_stream_sobol(j,k)),Z_CORDIV_sobol(i,j,k)),and(X2_stream_sobol(j,k),X1_stream_sobol(i,k)));

                %Z_stream_sobol(i,j,k) = JK_FF(J1_sobol(i,j,k),K_stream_sobol(i,j,k));
%                 Z_stream_vd(i,j,k) = JK_FF(J_stream_vd(i,j,k),K_stream_vd(i,j,k));
%                 Z_SS_vd(i,j,k) = JK_FF(J_SS_vd(i,j,k),K_SS_vd(i,j,k));
                Z_CORDIV_vd(i,j,k+1) = bitor(bitand(not(X2_stream_vandercorput(j,k)),Z_CORDIV_vd(i,j,k)),and(X2_stream_vandercorput(j,k),X1_stream_vandercorput(i,k)));

                Z_CORDIV_LF(i,j,k+1) = bitor(bitand(not(X2_stream_LF(j,k)),Z_CORDIV_LF(i,j,k)),and(X2_stream_LF(j,k),X1_stream_LF(i,k)));

            end


            Expected_div(i,j) = (i-1)/(j-1);
        %end
    %end
%             
            Err_sobol(i,j) = abs((sum(Z_stream_sobol(i,j,:))/N) - Expected_div(i,j));
%             Err_sobol_sq(i,j) = power((sum(Z_stream_sobol(i,j,:))/N) - Expected_div(i,j),2);
            Err_sobol_SS(i,j) = abs((sum(Z_SS_sobol(i,j,:))/N) - Expected_div(i,j));
            Err_sobol_CORDIV(i,j) = abs((sum(Z_CORDIV_sobol(i,j,1:N))/N) - Expected_div(i,j));

            Err_vd(i,j) = abs((sum(Z_stream_vd(i,j,:))/N) - Expected_div(i,j));
%             Err_vd_sq(i,j) = power((sum(Z_stream_vd(i,j,:))/N) - Expected_div(i,j),2);
            Err_vd_SS(i,j) = abs((sum(Z_SS_vd(i,j,:))/N) - Expected_div(i,j));
            Err_vd_CORDIV(i,j) = abs((sum(Z_CORDIV_vd(i,j,1:N))/N) - Expected_div(i,j));

            Err_LF(i,j) = abs((sum(Z_stream_LF(i,j,:))/N) - Expected_div(i,j));
            Err_LF_SS(i,j) = abs((sum(Z_SS_LF(i,j,:))/N) - Expected_div(i,j));
            Err_LF_CORDIV(i,j) = abs((sum(Z_CORDIV_LF(i,j,1:N))/N) - Expected_div(i,j));


        end
    end
end

MAE_Sobol_min_max = sum(Err_sobol,"all")/(128*255)      % (2^(bp-1)*((2^bp)-1))
% MSE_Sobol = sum(Err_sobol_sq,"all")/numel(find(Err_sobol_sq ~= 0))
MAE_Sobol_Saturated_Sub = sum(Err_sobol_SS,"all")/(128*255)
MAE_Sobol_CORDIV = sum(Err_sobol_CORDIV,"all")/(128*255)

MAE_VD_min_max = sum(Err_vd,"all")/(128*255)
MAE_VD_Saturated_Sub = sum(Err_vd_SS,"all")/(128*255)
MAE_VD_CORDIV = sum(Err_vd_CORDIV,"all")/(128*255)

MAE_LFSR_min_max = sum(Err_LF,"all")/(128*255)
MAE_LFSR_Saturated_Sub = sum(Err_LF_SS,"all")/(128*255)
MAE_LFSR_CORDIV = sum(Err_LF_CORDIV,"all")/(128*255)

toc
