function [emptypar, PARvars,PARSliceVars] = defaultPARVars()
emptypar = struct('info',[],'scn',[],'max',[],'orient',[],'special',[],'cardiac',[],'diffusion',[],'img',[]);
PARvars = {...
    ... % Variables from Prerelease 1. Feb 1, 2005
    '.    Patient name                       ','info','patient_name',0;
    '.    Examination name                   ','scn','exam_name',0;
    '.    Protocol name                      ','scn','protocol_name',0;
    '.    Examination date/time              ','info','exam_datetime',0;
    '.    Acquisition nr                     ','scn', 'acquisitin_num',1;
    '.    Reconstruction nr                  ','scn','recon_num',1;
    '.    Scan Duration [sec]                ','scn','scan_dur',1;
    '.    Max. number of cardiac phases      ','max','card_phs',1;
    '.    Max. number of echoes              ','max','num_echo',1;
    '.    Max. number of slices/locations    ','max','num_slices',1;
    '.    Max. number of dynamics            ','max','num_dynamics',1;
    '.    Max. number of mixes               ','max','num_mixes',1;
    '.    Image pixel size [8 or 16 bits]    ','scn','pix_bits',1;
    '.    Technique                          ','scn','technique',0;
    '.    Scan mode                          ','scn','scan_mode',0;
    '.    Scan resolution  (x, y)            ','scn','scan_res',1;
    '.    Scan percentage                    ','scn','scan_pct',1;
    '.    Recon resolution (x, y)            ','scn','recon_res',1;
    '.    Number of averages                 ','scn','NEX',1;
    '.    Repetition time [msec]             ','scn','rep_time',1;
    '.    FOV (ap,fh,rl) [mm]                ','scn','fov',1;
    '.    Slice thickness [mm]               ','scn','slicethk',1;
    '.    Slice gap [mm]                     ','scn','slicegap',1;
    '.    Water Fat shift [pixels]           ','scn','water_fat_shift',1;
    '.    Angulation midslice(ap,fh,rl)[degr]','orient','ang_midslice',1;
    '.    Off Centre midslice(ap,fh,rl) [mm] ','orient','off_ctr_midslice',1;
    '.    Flow compensation <0=no 1=yes> ?   ','special','flow_comp',1;
    '.    Presaturation     <0=no 1=yes> ?   ','special','presatuaration',1;
    '.    Cardiac frequency                  ','cardiac','cardiac_freq',1;
    '.    Min. RR interval                   ','cardiac','min_rr_int',1;
    '.    Max. RR interval                   ','cardiac','max_rr_int',1;
    '.    Phase encoding velocity [cm/sec]   ','cardiac','phase_enc_vel',1;
    '.    MTC               <0=no 1=yes> ?   ','special','mtc',1;
    '.    SPIR              <0=no 1=yes> ?   ','special','spir',1;
    '.    EPI factor        <0,1=no EPI>     ','special','epi_factor',1;
    '.    TURBO factor      <0=no turbo>     ','special','turbo_factor',1;
    '.    Dynamic scan      <0=no 1=yes> ?   ','special','dynamic_scan',1;
    '.    Diffusion         <0=no 1=yes> ?   ','diffusion','diffusion',1;
    '.    Diffusion echo time [msec]         ','diffusion','diffusion_echo',1;
    '.    Inversion delay [msec]             ','special','inversion_delay',1;
    ... % Variables for May 31, 2005
    '.    Series Type                        ','scn','series_type',0;
    '.    Patient position                   ','orient','patient_pos',0;
    '.    Preparation direction              ','orient','prep_dir',0;
    '.    Repetition time [ms]               ','scn','rep_time',1;
    '.    Diffusion echo time [ms]           ','special','diffusion_echo_time',1;
    ... % Variables for December 29, 2006 (release 2.1)
    '.    Max. number of diffusion values    ','special','max_num_diffusion_values',1;
    '.    Max. number of gradient orients    ','special','max_num_gradient_orients',1;
    ... % Variables for August 16, 2007 (release 2.5)
    '.    Number of label types   <0=no ASL> ','special','number_of_label_types',1;
    };

PARSliceVars  = {...
    '#  slice number                             (integer)',1,'info','slice_num';
    '#  echo number                              (integer)',1,'info','echo_num';
    '#  dynamic scan number                      (integer)',1,'info','dynamic_scan_num';
    '#  cardiac phase number                     (integer)',1,'info','cardiac_phase_num';
    '#  image_type_mr                            (integer)',1,'info','image_type_mr';
    '#  scanning sequence                        (integer)',1,'info','scan_seq';
    '#  index in REC file (in images)            (integer)',1,'info','idx_rec';
    '#  image pixel size (in bits)               (integer)',1,'info','pix_bits';
    '#  scan percentage                          (integer)',1,'info','scanpct';
    '#  recon resolution (x y)                   (2*integer)',2,'info','recon_res';
    '#  rescale intercept                        (float)',1,'vis','rescale_intercept';
    '#  rescale slope                            (float)',1,'vis','rescale_slope';
    '#  scale slope                              (float)',1,'vis','scale_slope';
    '#  window center                            (integer)',1,'notused','window_center';
    '#  window width                             (integer)',1,'notused','window_width';
    '#  image angulation (ap,fh,rl in degrees )  (3*float)',3,'orient','img_angulation';
    '#  image offcentre (ap,fh,rl in mm )        (3*float)',3,'orient','img_offcentre';
    '#  slice thickness (in mm )                 (float)',1,'info','slicethk';
    '#  slice gap (in mm )                       (float)',1,'info','slicegap';
    '#  image_display_orientation                (integer)',1,'orient','display_orientation';
    '#  slice orientation ( TRA/SAG/COR )        (integer)',1,'orient','slice_orientation';
    '#  fmri_status_indication                   (integer)',1,'special','fmri_status_indication';
    '#  image_type_ed_es  (end diast/end syst)   (integer)',1,'special','image_type_ed_es';
    '#  pixel spacing (x,y) (in mm)              (2*float)',2,'special','pix_spacing';
    '#  echo_time                                (float)',1,'special','echo_time';
    '#  dyn_scan_begin_time                      (float)',1,'special','dyn_scan_begin_time';
    '#  trigger_time                             (float)',1,'special','trigger_time';
    '#  diffusion_b_factor                       (float)',1,'special','diffusion_b_factor';
    '#  number of averages                       (integer)',1,'special','NEX';
    '#  image_flip_angle (in degrees)            (float)',1,'special','image_flip_angle';
    '#  cardiac frequency   (bpm)                (integer)',1,'special','cardiac_freq'
    '#  minimum RR-interval (in ms)              (integer)',1,'special','minRR';
    '#  maximum RR-interval (in ms)              (integer)',1,'special','maxRR';
    '#  TURBO factor  <0=no turbo>               (integer)',1,'special','turbo';
    '#  Inversion delay (in ms)                  (float)',1,'special','inversion_delay';
    ... % new columns for December 29, 2006 (release 2.1)
    '#  diffusion b value number    (imagekey!)  (integer)',1,'special','number_of_diffusion_b_factor'; 
    '#  gradient orientation number (imagekey!)  (integer)',1,'special','gradient_orientation_number';
    '#  contrast type                            (string)',1,'special','contrast_type';
    '#  diffusion anisotropy type                (string)',1,'special','diffusion_aniostropy_type';
    '#  diffusion (ap, fh, rl)                   (3*float)',3,'special','diffusion_ap_fh_lr';
    ... % new columns for August 16, 2007 (release 2.5)
    '#  label type (ASL)            (imagekey!)  (integer)',1,'special','label_type';
    };
