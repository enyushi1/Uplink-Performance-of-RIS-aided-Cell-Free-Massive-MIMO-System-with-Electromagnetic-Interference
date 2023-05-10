function [Hhat_LMMSE] = functionChannelEstimates_LMMSE(R_AP,HMean_Withoutphase,H,nbrOfRealizations,M,K,N,tau_p,p,Pset)

%LMMSE channel estimator for Cell-Free setup. The estimation is locally
%performed at the APs.Note that HMean should be without phase HMeanWithoutPhase
%And each AP is equipped with N antennas.
%This is version 1.0 (Last edited: 2020-04-17)


%INPUT:
%R_AP                 = Matrix with dimension N x N x M x K where (:,:,m,k) is
%                       the spatial correlation matrix between AP l and UE k 
%                       in setup n, normalized by the noise power
%HMeanWithoutPhase    = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,i,k) is the i^th realization of the channel mean
%                       between the n^th antenna of AP m and UE k(without random phase shifts)
%H                    = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,i,k) is the i^th channel realization
%                       between the n^th antenna of AP m and UE k (only used for generating y_p)
%nbrOfRealizations    = Number of realizations
%M                    = Number of APs
%K                    = Number of UEs 
%N                    = Number of antennas per AP
%p                    = 1xK vector, uplink power at each UE
%tau_p                = Pilot length
%Pset                 = Pilot allocation set
%
%OUTPUT:
%Hhat_LMMSE           = Matrix with dimension MN x nbrOfRealzations x K
%                       where (mn,i,k) is the i^th  realization of LMMSE channel estimate
%                       between the n^th antenna of AP m and UE k


%Prepare to store LMMSE channel estimates
Hhat_LMMSE = zeros(M*N,nbrOfRealizations,K);

%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(N,nbrOfRealizations,M,K) + 1i*randn(N,nbrOfRealizations,M,K));

%Store identity matrix of size N x N
eyeN = eye(N);

Lk = zeros(N,N,M,K);
for m = 1:M
    for k = 1:K
        
         Lk(:,:,m,k) = HMean_Withoutphase((m-1)*N+1:m*N,k)*(HMean_Withoutphase((m-1)*N+1:m*N,k))';
         
    end
end
 

Rp = R_AP + Lk;


for m = 1:M
    for k = 1:K
        
        yp = zeros(N,nbrOfRealizations);
        PsiInv_LMMSE = zeros(N,N);
        inds = Pset(:,k); 
        
        for z = 1:length(inds)
            
            yp = yp + sqrt(p(inds(z)))*tau_p*H((m-1)*N+1:m*N,:,inds(z));
            PsiInv_LMMSE = PsiInv_LMMSE + p(inds(z))*tau_p*Rp(:,:,m,inds(z));
            
        end
        
        yp = yp + sqrt(tau_p)*Np(:,:,m,k);
        PsiInv_LMMSE = PsiInv_LMMSE + eyeN;
        
        for z = 1:length(inds) 
            
            RPsi = Rp(:,:,m,inds(z))/PsiInv_LMMSE;
            Hhat_LMMSE((m-1)*N+1:m*N,:,inds(z)) = sqrt(p(inds(z)))*RPsi*yp;
                
        end
    end
end

       