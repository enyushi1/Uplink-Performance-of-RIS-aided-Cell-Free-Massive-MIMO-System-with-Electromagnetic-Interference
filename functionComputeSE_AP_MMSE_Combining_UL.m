function [SE_MMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat,H,C,tau_c,tau_p,nbrOfRealizations,N,K,M,p)
%Compute uplink SE for Cell-free mMIMO for the four different receiver
%cooperation levels, using either MR or MMSE/L-MMSE combining
%
%This function was developed as a part of the paper:
%
%Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%Competitive With MMSE Processing and Centralized Implementation,"
%IEEE Transactions on Wireless Communications, To appear.
%
%Download article: https://arxiv.org/abs/1903.10611
%
%This is version 1.0 (Last edited: 2019-03-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
%INPUT:
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%H                 = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the true collective channel from all
%                    BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences and number of UEs per cell
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%p1                = (Optional) Same as p but only for Level 1
%
%OUTPUT:
%SE_MR     = K x 4 matrix where the (k,n):th element is the uplink SE of 
%            UE k achieved with MR combining at cooperation level n
%SE_MMMSE  = Same as SE_MR but with MMSE or L-MMSE combining
%sumSE_SIC = Scalar with the sum SE achieved using MMSE-SIC combining



%If only one transmit power is provided, use the same for all the UEs
if length(p) == 1
   p = p*ones(K,1);
end

%If no specific Level 1 transmit powers are provided, use the same as for
%the other levels
if nargin<12
    p1 = p;
end


%Store identity matrices of different sizes
eyeN = eye(N);

%Compute the prelog factor
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MMSE = zeros(K,1);

%Compute sum of all estimation error correlation matrices at every BS
% C_tot = zeros(N,N,M);
C_tot = sum(C,4);
% for k = 1:K
%     
%     C_tot = C_tot + p(k)*(R(:,:,:,k)-B(:,:,:,k));
% 
% end


%Diagonal matrix with transmit powers and its square root
Dp = diag(p);
Dp12 = diag(sqrt(p));


%Prepare to save simulation results

signal_MMSE = zeros(M,K);
scaling_MMSE = zeros(M,K);
G_MMSE = zeros(M,M,K);
A = zeros(M,M,K);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    
    %-----------------Levels 1-3
    gp_MMSE = zeros(M,K,K);
    
    
    %Go through all APs
    for m = 1:M
        
        %Extract channel realizations from all UEs to AP l
        Hallj = reshape(H(1+(m-1)*N:m*N,n,:),[N K]);
        
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(1+(m-1)*N:m*N,n,:),[N K]);
        
        
        %Compute MR combining
        V_MMSE = ((Hhatallj*Dp*Hhatallj')+C_tot(:,:,m)+eyeN)\(Hhatallj*Dp);
        

        %Go through all UEs
        for k = 1:K
            
           
            
            
            %%MMSE combining
            v = V_MMSE(:,k); %Extract combining vector
            
            
            %Level 2 and Level 3
            signal_MMSE(m,k) = signal_MMSE(m,k) + (v'*Hallj(:,k))/nbrOfRealizations;
            gp_MMSE(m,:,k) = gp_MMSE(m,:,k) + (v'*Hallj)*Dp12;
            scaling_MMSE(m,k) = scaling_MMSE(m,k) + norm(v).^2/nbrOfRealizations;
            
            
   
           
        end
        
    end
    
    %Compute averaging of terms for Level 2 and Level 3
    for k = 1:K
        
        G_MMSE(:,:,k) = G_MMSE(:,:,k) + gp_MMSE(:,:,k)*gp_MMSE(:,:,k)'/nbrOfRealizations;
        
    end
    
end



%Compute SE for Level 2 and Level 3
for k = 1:K
    

    %With L-MMSE combining
    b = signal_MMSE(:,k);
    A(:,:,k) = G_MMSE(:,:,k) + diag(scaling_MMSE(:,k)) - p(k)*(b*b');
    SE_MMSE(k) = prelogFactor*real(log2(1+p(k)*b'*(A(:,:,k)\b)));  
    


end
