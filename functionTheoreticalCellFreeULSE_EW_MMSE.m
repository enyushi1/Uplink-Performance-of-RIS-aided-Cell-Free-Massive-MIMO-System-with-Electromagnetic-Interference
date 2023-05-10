function [SE_CC] = functionTheoreticalCellFreeULSE_EW_MMSE(HMean_Withoutphase,R_AP,Lk,Dk,Sigma,Phi_EW,M,K,N,tau_p,tau_c,p,Pset)       

%---This function is used to computes the theoretical uplink SE for
%Element-wise MMSE estimator.
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-05-12)

% M = CellFreeParameter.M;
% K = CellFreeParameter.K;
% N = CellFreeParameter.N;
% tau_p = CellFreeParameter.tau_p;
% tau_c = CellFreeParameter.tau_c;

%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K  where(:,:,m,k) is
%                       the spatial correlation matrix between AP m and UE k 
%                       in setup n, normalized by the noise power 
%                       normalized by the noise power
%Lk                   = Matrix with dimension MN x MN x K 
%  
%Phi                  = Matrix with dimension MN x MN x K
%                    
%Omega                = Matrix with dimension MN x MN x K
%
%A                    = Diagonal matrix with dimension M x M x K where (:,:,k)
%                       is the LSFD coefficients of UE k (when MMSE
%                       estimator is used.)
%M                   = Number of APs
%K                   = Number of UEs 
%p                   = 1xK vector, uplink power at each UE
%tau_p               = Pilot length
%tau_c               = Length of the coherence block
%Pset                = Pilot allocation set
%
%
%OUTPUT:
%
%SE_CC              = Vector with dimension K x 1 where (k) is the SE of UE k




%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store the results
Zk = zeros(M,M,K); %Zk
Ksi = zeros(M,M,K,K); %Ksi
X_p1 = zeros(M,M,K,K); %X(1)
Sigmak = zeros(M,M,K);
A_EW_MMSE_part1 = zeros(M,M,K,K);
CCterm1 = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm3 = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm2_P1 = zeros(K,K);
LK = zeros(M,M,K);
SE_CC = zeros(K,1);%Store the result
A = zeros(M,M,K);

%Go through all UEs
for k = 1:K
    
    for m = 1:M
        
        Zk(m,m,k) = trace(p(k)*tau_p*Dk(:,:,m,k)/Phi_EW(:,:,m,k)*Dk(:,:,m,k)+Lk(:,:,m,k));

        for l = 1:K  %Non-coherent interference (i=k')
            
            Ksi(m,m,k,l) = trace(R_AP(:,:,m,l)*Sigma(:,:,m,k))...
                + (HMean_Withoutphase((m-1)*N+1:m*N,k))'*R_AP(:,:,m,l)*HMean_Withoutphase((m-1)*N+1:m*N,k)...
                + (HMean_Withoutphase((m-1)*N+1:m*N,l))'*Sigma(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l)...
                + abs((HMean_Withoutphase((m-1)*N+1:m*N,k))'*HMean_Withoutphase((m-1)*N+1:m*N,l))^2;
                
           
           if any(l==Pset(:,k))&& l~=k %Coherent interference (If there is pilot contamination) 

               X_p1(m,m,k,l) = sqrt(p(k)*p(l))*tau_p*trace(Dk(:,:,m,k)/Phi_EW(:,:,m,k)*Dk(:,:,m,l));

           
           end
          
           
        end
        
        LK(m,m,k) = trace(Lk(:,:,m,k));
        Sigmak(m,m,k) = trace(Sigma(:,:,m,k)) + trace(Lk(:,:,m,k));
        
    end
    
end


%------------Calculation of LSFD coefficients
for k = 1:K 
    
    for l = 1:K  %Non-coherent interference (i=k')
        
        A_EW_MMSE_part1(:,:,k,l) = p(l)*Ksi(:,:,k,l);
        
        if any(l==Pset(:,k)) && l~=k
            
            Xkl = diag(X_p1(:,:,k,l));
            
            A_EW_MMSE_part1(:,:,k,l) = A_EW_MMSE_part1(:,:,k,l)...
                +p(l)*Xkl*(Xkl)';
            
        end
            %Check the case when l == k
        if l == k
            
            A_EW_MMSE_part1(:,:,k,k)=  A_EW_MMSE_part1(:,:,k,k) - p(k)* LK(:,:,k)*LK(:,:,k);
            
        end

            
    end

        
end
    


for k = 1:K
    
    A_EW_MMSE = sum(A_EW_MMSE_part1(:,:,k,:),4) + Sigmak(:,:,k);
    A(:,:,k) = diag(A_EW_MMSE\diag(Zk(:,:,k)));
    
end
%---------------------------------------------------------------------------------


for k = 1:K
    
    CCterm1(k) = trace(A(:,:,k)*Zk(:,:,k));
    CCterm3(k) = trace(A(:,:,k)'*Sigmak(:,:,k)*A(:,:,k));
    
    for l = 1:K  %Non-coherent interference (i=k')
    
        CCterm2_P1(k,l) = p(l)*trace(A(:,:,k)'*Ksi(:,:,k,l)*A(:,:,k));
        
        if any(l==Pset(:,k)) && l~=k
            
            CCterm2_P1(k,l)=  CCterm2_P1(k,l) + p(l)*abs(trace(A(:,:,k)*X_p1(:,:,k,l)))^2; 

            
        end
        
        if l == k
            
            CCterm2_P1(k,l)= CCterm2_P1(k,l) -  p(k)*trace(A(:,:,k)'*LK(:,:,k)*LK(:,:,k)*A(:,:,k));
        
        end
        
    end
    
end
    
CCterm2 = sum(CCterm2_P1,2);

for k=1:K
    
    SE_CC(k)= prelogFactor*real(log2(1+(p(k)*abs(CCterm1(k))^2)/((CCterm2(k))+CCterm3(k))));
    
end





