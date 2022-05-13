### FOR EXECUTION IN CASA ###


#=== Step Switches ===#
do_step1 = False #--- split out wanted SPW
do_step2 = False #--- uvcontsub
do_step3 = True #--- Image cube
do_step4 = False #--- immoments
#=====================#


#=== USER INPUTS ===#
targ='d_sma1'

"""H2CO321220 VALUES
wantSPW = '4,11,18,25'
LineFreeChans = '0:0~210;250~479,1:0~210;250~479,2:0~210;250~479,3:0~210;250~479'
targetLine ='H2CO321220'

H2CO322221 VALUES"""
wantSPW = '3,10,17,24'
LineFreeChans = '0:0~220;250~295;327~479,1:0~220;250~295;327~479,2:0~220;250~295;327~479,3:0~220;250~295;327~479'
targetLine ='H2CO322221'

"""H2CO303202 VALUES

wantSPW = '2,9,16,23'
LineFreeChans = '0:0~188;255~479,1:0~188;255~479,2:0~188;255~479,3:0~188;255~479'
targetLine ='H2CO'

CN3CN VALUES
wantSPW = '0,7,14,21'
LineFreeChans = '0:0~700;730~1390,1:0~700;730~1390,2:0~700;730~1390,3:0~700;730~1390'
targetLine ='CH3CN'

SiO VALUES
wantSPW = '1,8,15,22'
LineFreeChans = '0:0~220;300~479,1:0~220;300~479,2:0~220;300~479,3:0~220;300~479'
targetLine ='SiO'
"""

#--- imaging parameters (step 3)
useImsize = 1280
useCell = 0.025
useRestFreq = '218.475632GHz' #H2CO321220 == '218.760066GHz' #H2CO322221 =='218.475632GHz' #H2CO303202 == '218.222192GHz' #CH3CN == '220.709GHz' #SiO =='217.10498GHz'
useThreshold = '0.06mJy'
useNchan = 45
useStart = '30km/s'
useWid = '-0.671km/s'

#--- immoments parameters (step 4)
useChans = '40~60'
#myRegion = 'chans=44'

#===================#


#=== Step 1 ====#    
if do_step1:
    """ Split out target SPWs """
    split(vis='calibrated_final_all.ms', 
        spw=wantSPW, 
        outputvis='calibrated_final_all_'+targetLine+'.ms', 
        datacolumn='data')


    listobs('calibrated_final_all_'+targetLine+'.ms', 
        listfile='calibrated_final_all_'+targetLine+'.ms.listobs')

#=== Step 2 ====#
if do_step2:
    """ uvcontsub """
    uvcontsub(vis='calibrated_final_all_'+targetLine+'.ms', fitspw=LineFreeChans)

#=== Step 3 ====#
if do_step3:
    """ tclean """
    myVis = 'calibrated_final_all_'+targetLine+'.ms.contsub'

    imNameC = targ+'_'+targetLine                                       
    tclean(vis = myVis,                                      
             imagename = imNameC,    
             field = targ,                                                      
             stokes = 'I',                                                      
             spw = '',                                        
             outframe = 'LSRK',                                                 
             restfreq = useRestFreq,                                           
             specmode = 'cube',
             nchan = useNchan,
             width = useWid,
             start = useStart,                                             
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


#=== Step 4 ====#
if do_step4:
    """ immoments """
    imStats = imstat(imagename=targ+'_'+targetLine+'.image.fits')

    immoments(imagename=targ+'_'+targetLine+'.image.fits',
            moments=[0,1,2],
            axis='spectral',
            chans=useChans,
            includepix=[0.001,100],#!!!!! CHANGE THIS AS REQUIRED!
            outfile=targ+'_'+targetLine+'.image.mom'
            )

