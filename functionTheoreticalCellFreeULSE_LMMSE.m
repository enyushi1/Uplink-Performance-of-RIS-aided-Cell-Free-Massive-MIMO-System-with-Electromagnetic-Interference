function [SE_CC_LMMSE] = functionTheoreticalCellFreeULSE_LMMSE( R_AP,Rp,Phip,Omegap,HMean_Withoutphase,M,K,N,p,tau_p,tau_c,Pset)
%---This function is used to computes the theoretical uplink SE for
%LMMSE estimator.
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-20)


%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K  where(:,:,m,k) is
%                       the spatial correlation matrix between AP m and UE k 
%                       in setup n, normalized by the noise power 
%                       normalized by the noise power
%Lk                   = Matrix with dimension MN x MN x K 
%  
%Rp                   = Matrix with dimension MN x MN x K
%
%Phip                 = Matrix with dimension MN x MN x K
%   
%Omegap               = Matrix with dimension MN x MN x K
%M                   = Number of APs
%K                   = Number of UEs 
%N                   = Number of antennas per AP
%p                   = 1xK vector, uplink power at each UE
%tau_p               = Pilot length
%tau_c               = Length of the coherence block
%Pset                = Pilot allocation set
%
%
%OUTPUT:
%
%SE_CC_LMMSE         = Vector with dimension K x 1 where (k) is the SE of UE k




%Compute the pre-log factor
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c);

%Prepare to store the results
Zk_LMMSE = zeros(M,M,K);     %Zk_LMMSE
Ksi_LMMSE = zeros(M,M,K,K);  %Ksi
X_p1_LMMSE = zeros(M,M,K,K); %X(1)
X_p2_LMMSE = zeros(M,M,K,K); %X(2)
X_p3_LMMSE = zeros(M,M,K,K); %X(3)
A_LMMSE_part1 = zeros(M,M,K,K);
A_LMMSE = zeros(M,M,K);
CCterm1_LMMSE = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm3_LMMSE = zeros(K,1); %Store E{v^H_k h_k}= E{|v_k|^2}
CCterm2_P1_LMMSE = zeros(K,K);
SE_CC_LMMSE = zeros(K,1);%Store the result


for k = 1:K
    
    for m = 1:M
        
        Zk_LMMSE(m,m,k) = p(k)*tau_p*trace(Omegap(:,:,m,k));
        
        
        for l = 1:K  %Non-coherent interference (i=k')
            

            Ksi_LMMSE(m,m,k,l) = p(k)*tau_p*(trace(R_AP(:,:,m,l)*Omegap(:,:,m,k)) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Omegap(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l));

            
            if any(l==Pset(:,k))
                
                Smk = Rp(:,:,m,k)/Phip(:,:,m,k);
                Gama1 = tau_p*Smk*(Phip(:,:,m,k) - p(l)*tau_p*Rp(:,:,m,l))*Smk';
                Gama2 = Smk*R_AP(:,:,m,l)*Smk';

                 X_p1_LMMSE(m,m,k,l) = p(k)*(trace(R_AP(:,:,m,l)*Gama1) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Gama1*HMean_Withoutphase((m-1)*N+1:m*N,l))...
                    + p(k)*p(l)*tau_p^2*(abs(trace((Gama2')^(1/2)*(R_AP(:,:,m,l))^(1/2)))^2 + trace(R_AP(:,:,m,l)*Gama2) + abs(HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk'*HMean_Withoutphase((m-1)*N+1:m*N,l))^2 ...
                    + 2*real(trace((Gama2')^(1/2)*(R_AP(:,:,m,l))^(1/2))*HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))...
                    + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk'*R_AP(:,:,m,l)*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Gama2*HMean_Withoutphase((m-1)*N+1:m*N,l))...
                    - p(k)*tau_p*(trace(R_AP(:,:,m,l)*Omegap(:,:,m,k)) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Omegap(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l));
                
                    
                

%                  X_p1_LMMSE(m,m,k,l) = p(k)*(trace(R_AP(:,:,m,l)*Gama) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Gama*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     + p(k)*p(l)*tau_p^2*(abs(trace(Smk*R_AP(:,:,m,l)))^2 + (trace(R_AP(:,:,m,l)*Smk*R_AP(:,:,m,l)*Smk))+ abs(HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))^2 ...                   
%                     + 2*real(trace(Smk*R_AP(:,:,m,l))*HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     + 2*HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*R_AP(:,:,m,l)*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     - p(k)*tau_p*(trace(R_AP(:,:,m,l)*Omegap(:,:,m,k)) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Omegap(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l));
%                 
                
                
%                 X_p1_LMMSE(m,m,k,l) = p(k)*(trace(R_AP(:,:,m,l)*Gama) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Gama*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     + p(k)*p(l)*tau_p^2*((trace(Smk*R_AP(:,:,m,l)))^2 + trace(R_AP(:,:,m,l)*Smk*R_AP(:,:,m,l)*Smk)+ (HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))^2....
%                     + 2*(trace(Smk*R_AP(:,:,m,l))*HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     + 2*HMean_Withoutphase((m-1)*N+1:m*N,l)'*Smk*R_AP(:,:,m,l)*Smk*HMean_Withoutphase((m-1)*N+1:m*N,l))...
%                     - p(k)*tau_p*(trace(R_AP(:,:,m,l)*Omegap(:,:,m,k)) + HMean_Withoutphase((m-1)*N+1:m*N,l)'*Omegap(:,:,m,k)*HMean_Withoutphase((m-1)*N+1:m*N,l));
    
                
                
                
%                 X_p1_LMMSE(m,m,k,l) = p(k)*(trace(R_AP(:,:,m,l)*Gama) + trace(Lk(:,:,m,l)*Gama))...
%                     + p(k)*p(l)*tau_p^2*(abs(trace(Smk*R_AP(:,:,m,l)))^2 + trace(R_AP(:,:,m,l)*Smk*R_AP(:,:,m,l)*Smk)+abs(trace(Lk(:,:,m,l)*Smk))^2....
%                     + 2*real(trace(Smk*R_AP(:,:,m,l))*trace(Lk(:,:,m,l)*Smk)) + trace(Lk(:,:,m,l)*Smk*Smk*R_AP(:,:,m,l))...
%                     + trace(Lk(:,:,m,l)*Smk*R_AP(:,:,m,l)*Smk))- p(k)*tau_p*trace(Rp(:,:,m,l)*Omegap(:,:,m,k)) ;

                
                
                X_p2_LMMSE(m,m,k,l) = sqrt(p(l)*p(k))*tau_p*trace(Rp(:,:,m,l)*Rp(:,:,m,k)/Phip(:,:,m,k));
                
                X_p3_LMMSE(m,m,k,l) = p(l)*p(k)*tau_p^2*trace(Rp(:,:,m,l)*Rp(:,:,m,k)/Phip(:,:,m,k))*trace(Rp(:,:,m,l)*Rp(:,:,m,k)/Phip(:,:,m,k));
                
                clear Smk Gama
            
            end

        end
           
    end

end

%------------Calculation of LSFD coefficients
for k = 1:K 
    
    for l = 1:K  %Non-coherent interference (i=k')
        
        A_LMMSE_part1(:,:,k,l) = p(l)*Ksi_LMMSE(:,:,k,l);
        
        if any(l==Pset(:,k))
            
            dkl_LMMSE = diag(X_p2_LMMSE(:,:,k,l));
            A_LMMSE_part1(:,:,k,l) = A_LMMSE_part1(:,:,k,l)...
                + p(l)*X_p1_LMMSE(:,:,k,l) + p(l)*dkl_LMMSE*(dkl_LMMSE)' - p(l)*X_p3_LMMSE(:,:,k,l);
            
            
        end
        
        
    end
    
end

for k = 1:K
    
    bkl_LMMSE = diag(Zk_LMMSE(:,:,k));
    A_LMMSEp1 = sum(A_LMMSE_part1(:,:,k,:),4) - p(k)*bkl_LMMSE*(bkl_LMMSE)' + Zk_LMMSE(:,:,k);
    A_LMMSE(:,:,k) = diag(A_LMMSEp1\bkl_LMMSE);
    
end
%---------------------------------------------------------------------------------

for k = 1:K
    
    CCterm1_LMMSE(k) = trace(A_LMMSE(:,:,k)*Zk_LMMSE(:,:,k));
    CCterm3_LMMSE(k) = trace(A_LMMSE(:,:,k)'*Zk_LMMSE(:,:,k)*A_LMMSE(:,:,k));
    
    for l=1:K  %Non-coherent interference (i=k')
    
        CCterm2_P1_LMMSE(k,l) = p(l)*trace(A_LMMSE(:,:,k)'*Ksi_LMMSE(:,:,k,l)*A_LMMSE(:,:,k));
        
        if any(l==Pset(:,k))
            
            CCterm2_P1_LMMSE(k,l) = CCterm2_P1_LMMSE(k,l) + p(l)*trace(A_LMMSE(:,:,k)'*X_p1_LMMSE(:,:,k,l)*A_LMMSE(:,:,k))...
                + p(l)*abs(trace(A_LMMSE(:,:,k)*X_p2_LMMSE(:,:,k,l)))^2 - p(l)* trace(A_LMMSE(:,:,k)'*X_p3_LMMSE(:,:,k,l)*A_LMMSE(:,:,k));
            
        end
       
        
    end
    
end
    
CCterm2_LMMSE = sum(CCterm2_P1_LMMSE,2);

for k=1:K
    
    SE_CC_LMMSE(k)= prelogFactor*real(log2(1+(p(k)*abs(CCterm1_LMMSE(k))^2)/((CCterm2_LMMSE(k)) - p(k)*abs(CCterm1_LMMSE(k))^2 + CCterm3_LMMSE(k))));
    
end





