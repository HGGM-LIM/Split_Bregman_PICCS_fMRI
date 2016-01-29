function RAll   =   genRAll(dims,spars,P,radius,dir_out,rname,h2)

        [route nam ext]     =   fileparts(dir_out);
        RAll                =   zeros(dims);    
        NAll                =   size(RAll);
        N                   =   NAll(1:2);
        DN                  =   [N(1),1];     % data Size
        for time=1:dims(3)        
                time,
                pdf             =   genPDF(DN,P,spars,2,radius,0);    % generates the sampling PDF
                temp            =   genSampling_LIM(pdf,20,1);            % generates each time-point sampling pattern
				
                RAll(:,:,time)  =   temp;
                
                temp            =   fftshift(temp);
                figure(h2); subplot(1,2,1); 
                spy(temp), 100*nnz(temp)/N(1)
                title(['Sampling for time point ' num2str(time)]);
        end

        
end