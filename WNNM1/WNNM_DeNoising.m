function [E_Img]   =  WNNM_DeNoising( N_Img, O_Img, Par )


E_Img           = N_Img;                                                        % Estimated Image
[Height Width]  = size(E_Img);   
TotalPatNum     = (Height-Par.patsize+1)*(Width-Par.patsize+1);                 %Total Patch Number in the image
Dim             = Par.patsize*Par.patsize;  


[Neighbor_arr Num_arr Self_arr] =	NeighborIndex(N_Img, Par);                  % PreCompute the all the patch index in the searching window 
            NL_mat              =   zeros(Par.patnum,length(Num_arr));          % NL Patch index matrix
            CurPat              =	zeros( Dim, TotalPatNum );
            Sigma_arr           =   zeros( 1, TotalPatNum);            
            EPat                =   zeros( size(CurPat) );     
            W                   =   zeros( size(CurPat) );          
            
for iter = 1 : Par.Iter        
    E_Img             	=	E_Img + Par.delta*(N_Img - E_Img);
    [CurPat Sigma_arr]	=	Im2Patch( E_Img, N_Img, Par );                      % image to patch and estimate local noise variance            
    
    if (mod(iter-1,Par.Innerloop)==0)
        Par.patnum = Par.patnum-10;                                             % Lower Noise level, less NL patches
        NL_mat  =  Block_matching(CurPat, Par, Neighbor_arr, Num_arr, Self_arr);% Caculate Non-local similar patches for each 
        if(iter==1)
            Sigma_arr = Par.nSig * ones(size(Sigma_arr));                       % First Iteration use the input noise parameter
        end
    end       

    [EPat, W]  =  PatEstimation( NL_mat, Self_arr, Sigma_arr, CurPat, Par );   % Estimate all the patches
    E_Img      =  Patch2Im( EPat, W, Par.patsize, Height, Width );             
    E_Img(E_Img > 255) = 255; 
    E_Img(E_Img < 0) = 0; 
    if Par.bEvaluate
        fprintf( 'Ch = %2.3f, Iter = %2.3f, PSNR = %2.4f, SSIM = %2.4f \n', Par.ch, iter, psnr_MultiCh( E_Img, O_Img ), ssim( E_Img, O_Img ) );
    end
end
return;


