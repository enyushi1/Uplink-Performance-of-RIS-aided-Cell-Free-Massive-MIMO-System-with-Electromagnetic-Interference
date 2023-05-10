

%Empty workspace and close figures
clear
close all
clc

%Define the range of number of Access Points (APs)
M = 40:10:100;


%Number of UEs
K = 40;

%Number of antennas per AP
N = 4;

%Pilot length
tau_p = 10;


%Select length of coherence block
tau_c = 200;


%Uplink transmit power per UE (W)
p = 0.2; %200 mW
%Create the power vector for all UEs (The uplink power is the same
%(p)at each UE)
pv = p*ones(1,K );


%Select the number of setups with random AP/UE locations
nbrOfSetups = 50;
%Select the number of channel realizations per setup
nbrOfRealizations = 1000;


%Prepare to save simulation results 
userSE_MMSE_LSFD_MR_Combining= zeros(length(M),nbrOfSetups);
userSE_EW_MMSE_LSFD_MR_Combining= zeros(length(M),nbrOfSetups);
userSE_LMMSE_LSFD_MR_Combining= zeros(length(M),nbrOfSetups);

userSE_MMSE_LSFD_MMSE_Combining= zeros(length(M),nbrOfSetups);
userSE_EW_MMSE_LSFD_MMSE_Combining= zeros(length(M),nbrOfSetups);
userSE_LMMSE_LSFD_MMSE_Combining= zeros(length(M),nbrOfSetups);

userSE_TheoreticalMMSE_LSFD_MR_Combining = zeros(length(M),nbrOfSetups);
userSE_TheoreticalEW_MMSE_LSFD_MR_Combining = zeros(length(M),nbrOfSetups);
userSE_TheoreticalLMMSE_LSFD_MR_Combining = zeros(length(M),nbrOfSetups);


%Prepare to save simulation results 

SE_MMSE_LSFD_MR_Combining= zeros(K,nbrOfSetups,length(M));
SE_EW_MMSE_LSFD_MR_Combining= zeros(K,nbrOfSetups,length(M));
SE_LMMSE_LSFD_MR_Combining= zeros(K,nbrOfSetups,length(M));

SE_MMSE_LSFD_MMSE_Combining= zeros(K,nbrOfSetups,length(M));
SE_EW_MMSE_LSFD_MMSE_Combining= zeros(K,nbrOfSetups,length(M));
SE_LMMSE_LSFD_MMSE_Combining= zeros(K,nbrOfSetups,length(M));

SE_MMSE_LSFD_Theoretical_MR_Combining = zeros(K,nbrOfSetups,length(M));
SE_EW_MMSE_LSFD_Theoretical_MR_Combining = zeros(K,nbrOfSetups,length(M));
SE_LMMSE_LSFD_Theoretical_MR_Combining = zeros(K,nbrOfSetups,length(M));
       
%Go through the number of APs
for m = 1:length(M)

    A_singleLayer = reshape(repmat(eye(M(m)),1,K),M(m),M(m),K);
    for n=1:nbrOfSetups   
        
       %Deploy UEs and generate the covariance and mean matrices
       [R_AP,HMean_Withoutphase,~] = RandomAP_generateSetup_Rician_Multi_Antenna(M(m),K,N,1,1);

       %Create channel generations for each UE-AP pair
        [H,HMean] = functionChannelGeneration( R_AP,HMean_Withoutphase,M(m),K,N,nbrOfRealizations);
        disp(['ChannelGeneration of setup ' num2str(n)]);

       [Pset] = functionPilotAllocation( R_AP,HMean_Withoutphase,A_singleLayer,M(m),K,N,tau_p,pv);
       disp(['PilotAllocation of setup ' num2str(n)]);

       [Rp,Dk,Lk,Phi,Phip,Phi_EW,Omega,Omegap,Sigma,C_MMSE,C_EW_MMSE,C_LMMSE,C_LS] = functionMatrixGeneration(R_AP,HMean_Withoutphase,M(m),K,N,tau_p,pv,Pset);
       disp(['MatrixGeneration of setup ' num2str(n)]);


       [Hhat_MMSE] = functionChannelEstimates_MMSE(R_AP,HMean,H,nbrOfRealizations,M(m),K,N,tau_p,pv,Pset);
       disp(['MMSE Estimate of setup ' num2str(n)]);
       
       [Hhat_EW_MMSE] = functionChannelEstimates_EW_MMSE(R_AP,HMean,H,nbrOfRealizations,M(m),K,N,tau_p,pv,Pset);
       disp(['EW_MMSE Estimate of setup ' num2str(n)]);

       [Hhat_LMMSE] = functionChannelEstimates_LMMSE(R_AP,HMean_Withoutphase,H,nbrOfRealizations,M(m),K,N,tau_p,pv,Pset);
       disp(['LMMSE Estimate of setup ' num2str(n)]);
       
      

      
       %Second layer decoding Spectral Efficiency (SE) computation
       %SE with MMSE estimator

       [SE_MR_MMSE] = functionComputeMonteCarloSE_UL(Hhat_MMSE,H,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MR Combining Monte Compution of MMSE Estimate of setup ' num2str(n)]);
       [SE_MMSE_Combining_MMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat_MMSE,H,C_MMSE,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MMSE Combining Monte Compution of MMSE Estimate of setup ' num2str(n)]);
       [SE_CC_MMSE] = functionTheoreticalCellFreeULSE_MMSE(HMean_Withoutphase,R_AP,Lk,Omega,Phi,M(m),K,N,tau_p,tau_c,pv,Pset); 
       disp(['MR Combining Theoretical SE of MMSE Estimate of setup ' num2str(n)]);

       [SE_MR_EW_MMSE] = functionComputeMonteCarloSE_UL(Hhat_EW_MMSE,H,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MR Combining Monte Compution of EW_MMSE Estimate of setup ' num2str(n)]);
       [SE_MMSE_Combining_EW_MMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat_EW_MMSE,H,C_EW_MMSE,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MMSE Combining Monte Compution of EW_MMSE Estimate of setup ' num2str(n)]);
       [SE_CC_EW_MMSE] = functionTheoreticalCellFreeULSE_EW_MMSE(HMean_Withoutphase,R_AP,Lk,Dk,Sigma,Phi_EW,M(m),K,N,tau_p,tau_c,pv,Pset);
       disp(['MR Combining Theoretical SE of EW_MMSE Estimate of setup ' num2str(n)]);

       [SE_MR_LMMSE] = functionComputeMonteCarloSE_UL(Hhat_LMMSE,H,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MR Combining Monte Compution of LMMSE Estimate of setup ' num2str(n)]);
       [SE_MMSE_Combining_LMMSE] = functionComputeSE_AP_MMSE_Combining_UL(Hhat_LMMSE,H,C_LMMSE,tau_c,tau_p,nbrOfRealizations,N,K,M(m),pv);
       disp(['MMSE Combining Monte Compution of LMMSE Estimate of setup ' num2str(n)]);
       [SE_CC_LMMSE] = functionTheoreticalCellFreeULSE_LMMSE(R_AP,Rp,Phip,Omegap,HMean_Withoutphase,M(m),K,N,pv,tau_p,tau_c,Pset);
       disp(['MR Combining Theoretical SE of LMMSE Estimate of setup ' num2str(n)]);

     

       userSE_MMSE_LSFD_MR_Combining(m,n) = mean(SE_MR_MMSE,1);
       userSE_EW_MMSE_LSFD_MR_Combining(m,n) = mean(SE_MR_EW_MMSE,1);
       userSE_LMMSE_LSFD_MR_Combining(m,n) = mean(SE_MR_LMMSE,1);
 
       
       userSE_MMSE_LSFD_MMSE_Combining(m,n) = mean(SE_MMSE_Combining_MMSE,1);
       userSE_EW_MMSE_LSFD_MMSE_Combining(m,n) = mean(SE_MMSE_Combining_EW_MMSE,1);
       userSE_LMMSE_LSFD_MMSE_Combining(m,n) = mean(SE_MMSE_Combining_LMMSE,1);
  
       
       userSE_TheoreticalMMSE_LSFD_MR_Combining(m,n) = mean(SE_CC_MMSE,1);
       userSE_TheoreticalEW_MMSE_LSFD_MR_Combining(m,n) = mean(SE_CC_EW_MMSE,1);
       userSE_TheoreticalLMMSE_LSFD_MR_Combining(m,n) = mean(SE_CC_LMMSE,1);
   
       
       SE_MMSE_LSFD_MR_Combining(:,n,m) = SE_MR_MMSE;
       SE_EW_MMSE_LSFD_MR_Combining(:,n,m) = SE_MR_EW_MMSE;
       SE_LMMSE_LSFD_MR_Combining(:,n,m) = SE_MR_LMMSE;

       
       SE_MMSE_LSFD_MMSE_Combining(:,n,m)= SE_MMSE_Combining_MMSE;
       SE_EW_MMSE_LSFD_MMSE_Combining(:,n,m) = SE_MMSE_Combining_EW_MMSE;
       SE_LMMSE_LSFD_MMSE_Combining(:,n,m)= SE_MMSE_Combining_LMMSE;



       SE_MMSE_LSFD_Theoretical_MR_Combining(:,n,m) = SE_CC_MMSE;
       SE_EW_MMSE_LSFD_Theoretical_MR_Combining(:,n,m) = SE_CC_EW_MMSE;
       SE_LMMSE_LSFD_Theoretical_MR_Combining(:,n,m) = SE_CC_LMMSE;
 
   
        disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    
    end 
    clear R_AP H HMean Hhat_MMSE Hhat_LMMSE Hhat_LS
%     clear R_AP H HMean Hhat_LMMSE
    disp([num2str(M(m)) ' APs out of ' num2str(M(end))]);
end

figure()
m1 = plot(M,mean(userSE_MMSE_LSFD_MMSE_Combining,2),'r-.');
hold on
m2 = plot(M,mean(userSE_EW_MMSE_LSFD_MMSE_Combining,2),'b--.');
hold on
m3 = plot(M,mean(userSE_LMMSE_LSFD_MMSE_Combining,2),'k:.');
hold on
plot(M,mean(userSE_MMSE_LSFD_MR_Combining,2),'r-.');
hold on
plot(M,mean(userSE_EW_MMSE_LSFD_MR_Combining,2),'b--.');
hold on
plot(M,mean(userSE_LMMSE_LSFD_MR_Combining,2),'k:.');
hold on
m4 = plot(M,mean(userSE_TheoreticalMMSE_LSFD_MR_Combining,2),'rx');
hold on
plot(M,mean(userSE_TheoreticalEW_MMSE_LSFD_MR_Combining,2),'bx');
hold on
plot(M,mean(userSE_TheoreticalLMMSE_LSFD_MR_Combining,2),'kx');
hold on



grid on
grid minor
xlabel('Number of APs (M)','Interpreter','Latex');
ylabel('Average SE [bit/s/Hz]','Interpreter','Latex');
set(gca,'FontSize',12);
legend([m1,m2,m3,m4],{'MMSE estimator','EW-MMSE estimator','LMMSE estimator','Analytical results (MR)'},'Interpreter','Latex','Location','Southeast')


