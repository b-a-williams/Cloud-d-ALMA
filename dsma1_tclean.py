### FOR EXECUTION IN CASA ###

targ='d_sma1'

do_ALL_CONT_IMG=True
#===== IMAGE    
if do_ALL_CONT_IMG:


        myVis = 'calibrated_final_all.ms'
        bigAggSPWstring='0:0~700;730~1390,1:0~220;300~479,2:0~188;255~479,3:0~220;250~295;327~479,4:0~210;250~479,5:0~958,6:0~958,7:0~700;730~1390,8:0~220;300~479,9:0~188;255~479,10:0~220;250~295;327~479,11:0~210;250~479,12:0~958,13:0~958,14:0~700;730~1390,15:0~220;300~479,16:0~188;255~479,17:0~220;250~295;327~479,18:0~210;250~479,19:0~958,20:0~958,21:0~700;730~1390,22:0~220;300~479,23:0~188;255~479,24:0~220;250~295;327~479,25:0~210;250~479,26:0~958,27:0~958'

        #--- DATA Notes
        """
        central freq 225.5GHz (ish)
        resolution ~0.09arcsec
        FoV ~ 28arcsec

        cell size = 0.025arcsec
        imsize = 1120.0

        """

        #--- imaging parameters
        useImsize = 1280
        useCell = 0.025
        useRestFreq = '225.5GHz'
        useThreshold = '10.0uJy'


        imNameC = targ+'_ALL_LINEFREE_CONT'                                          
        tclean(vis = myVis,                                      
             imagename = imNameC,    
             field = targ,                                                      
             stokes = 'I',                                                      
             spw = bigAggSPWstring,                                        
             outframe = 'LSRK',                                                 
             restfreq = useRestFreq,                                           
             specmode = 'mfs',                                                  
             imsize = [useImsize, useImsize],                                               
             cell = str(useCell)+'arcsec',                                                
             deconvolver = 'multiscale',                                        
             scales = [0, 4, 12], #in pixels the size of scales to probe (this = point source, syth beam size & 3 times synth beam)                    
             niter = 100000,                                                     
             weighting = 'briggs',                                              
             robust = 0.5,                                                      
             usemask = 'auto-multithresh',                                      
             sidelobethreshold = 3.0,                                           
             lownoisethreshold = 2.0,                                           
             minbeamfrac =  0.15,                                               
             noisethreshold = 4.5,                                              
             gridder = 'standard',                                              
             pbcor = True,     
             threshold = useThreshold,                                             
             interactive = True,                                                
             restoringbeam = 'common'                                           
             ) 
