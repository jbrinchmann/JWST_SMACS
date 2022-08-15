pro _time_c_vs_nonc, s, models = models, $
                      dim = dim, oh = oh, results = results, $
                      etaha=etaha, etao2=etao2, mugas = mugas, $
                      te=te, dgr=dgr, area=area_in, $
                      lines=lines, info=info, lgm=lgm


   inds = [8, 10, 11]

   print, 'Doing a run-in run'
   gas_fit_one_cl01, s, 8, models=models, dim=dim, oh=oh, $
                     res=res, etaha=etaha, etao2=etao2, mugas=mugas, $
                     te=te, dgr=dgr, area=area_in, lines=lines, $
                     info=info, prob=p, /steep, chisq=chisq, /ks_plane, $
                     lgm=lgm, /kss_plane, /ksoh_plane
   

   print, 'Doing the non-C version'
   starttime = systime(/sec)
   for i=0L, n_elements(inds)-1 do begin
      gas_fit_one_cl01, s, inds[i], models=models, dim=dim, oh=oh, $
                        res=res, etaha=etaha, etao2=etao2, mugas=mugas, $
                        te=te, dgr=dgr, area=area_in, lines=lines, $
                        info=info, prob=p, /steep, chisq=chisq, /ks_plane, $
                        lgm=lgm, /kss_plane, /ksoh_plane
   endfor
   endtime = systime(/sec)
   time_nonC = endtime-starttime
   print, format='("That took ",F0.2," seconds")', time_nonC

   print, 'Doing the C-version'
   starttime = systime(/sec)
   for i=0L, n_elements(inds)-1 do begin
      gas_fit_one_cl01, s, inds[i], models=models, dim=dim, oh=oh, $
                        res=res, etaha=etaha, etao2=etao2, mugas=mugas, $
                        te=te, dgr=dgr, area=area_in, lines=lines, $
                        info=info, prob=p, /steep, chisq=chisq, /ks_plane, $
                        lgm=lgm, /kss_plane, /ksoh_plane, /c_code
   endfor
   endtime = systime(/sec)
   time_C = endtime-starttime
   print, format='("That took ",F0.2," seconds")', time_C


end

pro jwst_fit_one_cl01, z, data, err, models = models, $
                      dim = dim, oh = oh, results = results, $
                      etaha=etaha, etao2=etao2, mugas = mugas, $
                      te=te, $
                      dgr=dgr, area=area_in, $
                      lines=lines, info=info, $
                      probability=p, apply_xsi_prior=apply_xsi_prior, $
                      mu_dependent_n=mu_dependent_n, n_slope=n_slope, $
                      steep_cf00=steep_cf00, non_sdss=non_sdss, chisq=chisq, $
                      apply_mu_prior=apply_mu_prior, ks_plane=ks_plane, $
                      lgm=lgm, kss_plane=kss_plane, ksoh_plane=ksoh_plane, $
                      c_code=c_code, fast_planes=fast_planes, $
                       no_theoretical_error=no_theoretical_error, $
                       apply_tauv_prior=apply_tauv_prior
;+
; NAME: GAS_FIT_ONE_CL01
;
; PURPOSE:
;    This is a routine - branched off GENERIC_CHANGER_OF_FITS for
;    running a CL01 fit. It keeps the possibility to create multi-D
;    PDFs if needed.
;
;
; CATEGORY: Model fitting routines
;
; CALLING SEQUENCE: 
;       GAS_FIT_ONE_CL01, S, IND [KEYWORDS=KEYWORDS]
;
; INPUTS: 
;
;   IND: The index of the object to fit in the SDSS structure. 
;     S: The SDSS structure as returned from READ_SDSS_INFO 
; 
;
; OPTIONAL INPUTS:
;     None
; KEYWORD PARAMETERS: 
;
;             DIM: The dimension structure belonging to the
;                  models. For all practial uses this is the
;                  interpolated model grid dimensions
;           ETAHA: Eta(Ha) for the models
;           ETAO2: Eta([O II]) for the models
;              OH: The 12 + Log O/H values for each model
;           MUGAS: The log mu_gas values for each model
;             DGR: The log dust-to-gas ratio for each model
;            INFO: Information structure relevant for the different
;                  methods. See the discussion of the individual
;                  methods under PROCEDURE below.
;           LINES: A string array with the names of the individual
;                  lines going into the fit.
;          MODELS: A model array from cl01_read_models
;         RESULTS: The result structure. Discussed under OUTPUTS below
;          METHOD: A string array which controls how the loops are
;                  constrained. The possible values are:
;                  'ChangingLines', 'AGN', 'ModelAdd'. The different
;                  methods are discussed in detail under PROCEDURE.
;               P: Set this keyword to a variable that will keep the
;                  probability distribution from the fits.
; APPLY_XSI_PRIOR: Set this keyword if you wish to apply a Gaussian
;                  prior on the CL01 xsi parameter. The default is to
;                  center the Gaussian on xsi=0.3 with a width of
;                  0.08. If APPLY_XSI_PRIOR is set to a structure with
;                  keys MEAN_XSI adn SIGMA_XSI, these keys can be used
;                  to set the prior to the specified values. These are
;                  not recorded in the result structure so keeping
;                  track of the values used rests on the user.
;  APPLY_MU_PRIOR: The same for MU. The default is a central value of
;                  0.65 with a sigma of 0.1
;  MU_DEPENDENT_N: Set this keyword to apply a fit where the slop
;                  depends on the value of mu. This is a lot slower -
;                  a factor of 8 slower in fact, because then all mu
;                  values must be included in the fit. The exact
;                  equation used is defined in
;                  CL01_MU_DEPENDENT_N. The default is to use a mu
;                  dependent n from Chevallard et al. However
;                  MU_DEPENDENT_N can also be set to a structure and
;                  if this has a key WILD07 set to 1, then the simpler
;                  equation from Wild et al (2007) is used where the
;                  _effective_ slope is mu dependent but not that of
;                  the ISM alone.
;     FAST_PLANES: In this case only the covariance matrix and center
;                  of the PDF is calculated, approximating it with a
;                  Gaussian. 
;
; OUTPUTS:
;   RESULTS: The result structure contains all the results of the
;            fits. It contains the following keys
;
;     BIN_Z    : The likelihoods for Log Z [N_REPETITIONS x N_L
;     BIN_U    : The likelihoods for Log U
;     BIN_XSI  : The likelihoods for Xsi
;     BIN_TAUV : The likelihoods for TauV
;     BIN_OH   : The likelihoods for 12 + Log O/H. These are
;                calculated in bins of width 0.1 dex from 7.5 to 9.5
;                whose centers are stroed in XOH
;     BIN_ETAHA: The likelihoods for Log Eta(Ha). Calculated in bins
;                of width 0.1 from 5.5 to 7.9
;     BIN_ETAO2: The same for [O II]
;     MIN_CHISQ: The minimum chi^2 found for the fits.
;
;           XOH: The x-position of the bins in 12 + Log O/H
;        XETAHA: The x-position of the bins in Log Eta(Ha)
;        XETAO2: The same for [O II]
;
; OPTIONAL OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:
;
;       Jul 31, 2022, jarle (jarle@strw.leidenuniv.nl)
;	  Created from gas_fit_one_cl01
;
;       Feb 14, 2013, jarle (jarle@astro.up.pt)
;	  Created from GENERIC_CHANGER_OF_FITS
;
;       Jan 31, 2006, Jarle Brinchmann (jarle@astro.up.pt)
;	  Added support for variable number of lines.
;
;       Nov 20, 2004, Jarle Brinchmann (jarle@astro.up.pt)
;	   Documented routine
;
;
;-

   if (n_elements(lines) eq 0) then begin
      lines = ['LOII_3727', 'LHb_4861', 'LHg_4340', 'LOIII_5007'] 
      
   endif

   ;;-----------------------
   ;; Read in the models
   ;;-----------------------

   if (n_elements(models) eq 0) Then begin
      if (keyword_set(mu_dependent_n)) then begin
         if (size(mu_dependent_n, /TNAME) eq 'STRUCT') then $
          wild07 = mu_dependent_n.wild07

         models=cl01_read_models(lines, dim=dim, mu_dependent_n=1, wild07=wild07)
      endif else begin
         models=cl01_read_models(lines, dim=dim, logmuind=1, steep=steep_cf00)
;         stop
      endelse
   endif
   if (size(models, /type) ne 5) then models = double(models)

   if (n_elements(oh) eq 0) then $
    oh=hdf_to_struct('/data2/jarle/HDFFiles/OH_20p0.hdf')

   if (keyword_set(mu_dependent_n)) then begin
      yoh =oh.time___1[*,*,*,*,*]
   endif else begin
      yoh=reform(oh.time___1[*,*,*,1,*])
   endelse


   ;;--------------------------------------------------------------------
   ;; In reading below we need to make sure that our data are the same
   ;; shape as our variables (mugas, etaha, etao2). This is basically
   ;; a check to see whether the mu dimension is the same or not. The
   ;; reading above will keep all mu values - but this is only
   ;; necessary when using a mu-depdendent n.
   ;;
   ;; Note also that we apply 'reform' to remove the single dimension
   ;; because that is done as part of the chisq calculation as well. 
   ;;
   ;;--------------------------------------------------------------------
   dims_m = size(models, /dimen)
   if (keyword_set(mu_dependent_n) and dims_m[3+1] eq 1) then begin
      print, 'Seemingly incorrect input models - check your input!'
      return
   endif

   if (n_elements(etaha) eq 0) then begin
      etaha=get_interpolated_parameter('EtaHa', timeind=1, sfhind=2, dim=dim, /new) 
      dims_v = size(etaha, /dimen)

      ;; Logmu is dimension #4, ie index 3 but recall that the lines
      ;; end up as the first dimension of models - hence the +1
      if (dims_v[3] ne dims_m[3+1]) then begin
         if (dims_m[3+1] eq 1) then begin
            ;; Easy peasy.
            etaha = reform(etaha[*, *, *, 1, *])
         endif else begin
            print, 'I cannot figure out the dimensions!'
            stop
         endelse
      endif 
   endif

   if (n_elements(etao2) eq 0) then begin
      etao2=get_interpolated_parameter('EtaOII', timeind=1, sfhind=2, dim=dim, /new) 

      dims_v = size(etao2, /dimen)
      ;; Logmu is dimension #4, ie index 3 but recall that the lines
      ;; end up as the first dimension of models - hence the +1
      if (dims_v[3] ne dims_m[3+1]) then begin
         if (dims_m[3+1] eq 1) then begin
            ;; Easy peasy.
            etao2 = reform(etao2[*, *, *, 1, *])
         endif else begin
            print, 'I cannot figure out the dimensions!'
            stop
         endelse
      endif 
   endif

   if (keyword_set(mu_dependent_n) and n_elements(n_slope) eq 0) then begin
      cl01_dim_to_array, dim, adim
      mu = 10.0^adim.logmu
      ;; Chevallard et al (in prep) - this suppresses the scatter. 
      n_slope = 2.0/(1 + mu*adim.tauv)
;      stop
   endif


   if (n_elements(mugas) eq 0) then begin
      mugas=hdf_to_struct('/data/jarle/HDFFiles/Mugas_20p0.hdf')
      mugas=alog10(mugas.time___1[*,*,*,*,*])
      bad = where(finite(mugas) eq 0, n_bad)
      if (n_bad gt 0) then mugas[bad] = -99.

      dims_v = size(mugas, /dimen)
      ;; Logmu is dimension #4, ie index 3 but recall that the lines
      ;; end up as the first dimension of models - hence the +1
      if (dims_v[3] ne dims_m[3+1]) then begin
         if (dims_m[3+1] eq 1) then begin
            ;; Easy peasy.
            mugas = reform(mugas[*, *, *, 1, *])
         endif else begin
            print, 'I cannot figure out the dimensions!'
            stop
         endelse
      endif 

   endif
   dims_v = size(mugas, /dimen)


   if (n_elements(te) eq 0) then begin
      te = hdf_to_struct('/data/jarle/HDFFiles/Te_20p0.hdf')
      te=te.time___1[*,*,*,*,*]
      bad = where(finite(te) eq 0, n_bad)
      if (n_bad gt 0) then te[bad] = -99.

      dims_v = size(te, /dimen)
      ;; Logmu is dimension #4, ie index 3 but recall that the lines
      ;; end up as the first dimension of models - hence the +1
      if (dims_v[3] ne dims_m[3+1]) then begin
         if (dims_m[3+1] eq 1) then begin
            ;; Easy peasy.
            te = reform(te[*, *, *, 1, *])
         endif else begin
            print, 'I cannot figure out the dimensions!'
            stop
         endelse
      endif 

   endif
   dims_v = size(te, /dimen)


   if (n_elements(dgr) eq 0) then begin
      dgr=hdf_to_struct('/data/jarle/HDFFiles/DustToGas_20p0.hdf')
      dgr=alog10(dgr.time___1[*,*,*,*,*])
      bad = where(finite(dgr) eq 0, n_bad)
      if (n_bad gt 0) then dgr[bad] = -99.

      dims_v = size(dgr, /dimen)
      ;; Logmu is dimension #4, ie index 3 but recall that the lines
      ;; end up as the first dimension of models - hence the +1
      if (dims_v[3] ne dims_m[3+1]) then begin
         if (dims_m[3+1] eq 1) then begin
            ;; Easy peasy.
            dgr = reform(dgr[*, *, *, 1, *])
         endif else begin
            print, 'I cannot figure out the dimensions!'
            stop
         endelse
      endif 

   endif
   dims_v = size(dgr, /dimen)


   

   
   ;;------------------------------------------
   ;; Create the data array for the lines
   ;;------------------------------------------

   n_lines = n_elements(lines)

   lumname = ['LOII_3727', 'LNeIII_3869', 'LHd_4101', 'LOIII_4363', 'LHg_4340', 'LHb_4861', 'LOIII_4959', $
              'LOIII_5007', 'LOI_6300', 'LHa_6563', 'LNII_6584', 'LSII_6716', 'LSII_6731', 'LSII_6724']
   tagname = ['oii3727', 'neiii3869', 'hd', 'oiii4363', 'hg', 'hb', 'oiii4959', 'oiii5007', $
              'oi6300', 'ha', 'nii6584', 'sii6716', 'sii6731', 'sii6724']
   wavelengths = [3727, 3869, 4101, 4363, 4340, 4861, 4959, 5007, 6300, 6563, 6584, 6716, 6731, 6724.]

   if (n_elements(data) eq 0) then begin
      data = dblarr(n_lines)
      err = dblarr(n_lines)
      fill_data = 1
   endif else fill_data = 0

;   tn = tag_names(s)
   lambdas = fltarr(n_lines)
   for i = 0, n_lines-1 do begin
      ii = where(lumname eq lines[i], num)

      if (num eq 0) then stop
      lambdas[i] = wavelengths[ii[0]]
      if (num eq 0) then begin
         print, 'I could not find a match for '+lines[i]
         return
      endif
 
   endfor         

   ;;---------------------------------------------------------------
   ;; Did we have 5007 included?
   ;; Here we have chance to check for bad 4959/5007 and replace as
   ;; needed. NOT CURRENTLY IMPLEMENTED
   ;;---------------------------------------------------------------
   is_5007 = where(lambdas eq 5007, n_o3)
   

   ;;---------------------------------------------
   ;; Convenient variables for [O III] fitting
   ;;---------------------------------------------
   o3flux = data[is_5007]
   o3ind = where(lines  eq 'LOIII_5007', n_o3)

   ;;---------------------------
   ;; Add a theoretical error
   ;;---------------------------
   if (keyword_set(no_theoretical_error)) then $
    theo_err = 0.0 $
   else $
    theo_err = 0.04
   err_bak = err
   err = sqrt(err^2 + (theo_err*data)^2)
      
   bad = where(data eq 0.0, n_bad)
   if (n_bad gt 0) then err[bad] = 1.0d30




   ;;-------------------------------------------------------------------
   ;; The data & the models have been read in. Now we can focus on
   ;; the actual fitting loop.
   ;;-------------------------------------------------------------------

   dims = size(models, /dimen)

   bin_Z = dblarr(n_elements(dim.logz))
   bin_U = dblarr(n_elements(dim.logU))
   bin_xsi = dblarr(n_elements(dim.xsi))
   bin_tauv = dblarr(n_elements(dim.tauv))
   if (keyword_set(mu_dependent_n)) then $
    bin_logmu = dblarr(n_elements(dim.logmu))


   ;; Calculate Chi^2!
   c = cl01_chisq(models, data, err)

;   stop

   tmp = !except                ; Turn off underflow warnings!
   !except = 0
   p = exp(-c.chi2/2.0)
   p = p/total(p)
   
   dum = 1.0+1.0
   !except = tmp

   if (keyword_set(apply_xsi_prior)) then begin
      ;; Create a prior in xsi but 1s everywhere else.
      prior = p*0.0+1
      if (size(apply_xsi_prior, /tname) eq 'STRUCT') then begin
         ;; The user has provided the settings for the xsi prior. 
         mean_xsi = struct_var(apply_xsi_prior, 'MEAN_XSI')
         sigma_xsi = struct_var(apply_xsi_prior, 'SIGMA_XSI')
      endif else begin
         mean_xsi = 0.3
         sigma_xsi = 0.08
      endelse
      
      ;; Define the prior as a Gaussian.
      prior_xsi = exp(-(dim.xsi-mean_xsi)^2/(2.0*sigma_xsi^2))/(2*sqrt(!pi*sigma_xsi))
      
      if (keyword_set(mu_dependent_n)) then begin
         prior_xsi = transpose(rebin(prior_xsi, 9, 33, 24, 8, 24), [1, 0, 2, 3, 4])
;            stop
      endif else begin
         prior_xsi = transpose(rebin(prior_xsi, 9, 33, 24, 24), [1, 0, 2, 3])
      endelse
      p = p*prior_xsi
      p = p/total(p)
   endif


   if (keyword_set(apply_tauv_prior)) then begin
      ;; Create a prior in xsi but 1s everywhere else.
      prior = p*0.0+1
      if (size(apply_tauv_prior, /tname) eq 'STRUCT') then begin
         ;; The user has provided the settings for the tau prior. 
         mean_tauv = struct_var(apply_tauv_prior, 'MEAN_TAUV')
         sigma_tauv = struct_var(apply_tauv_prior, 'SIGMA_TAUV')

         ;; Define the prior as a Gaussian.
         prior_tauv = exp(-(dim.tauv-mean_tauv)^2/(2.0*sigma_tauv^2))/(2*sqrt(!pi*sigma_tauv))

      endif else begin
         ;; Flat prior up to tauv=1
         prior_tauv = fltarr(n_elements(dim.tauv))
         ok = where(dim.tauv le 1.0)
         prior_tauv[ok] = 1.0
      endelse
      
;      stop
      
      if (keyword_set(mu_dependent_n)) then begin
         prior_tauv = transpose(rebin(prior_tauv, 24, 33, 9, 8, 24), [1, 2, 3, 0, 4])
;            stop
      endif else begin
         prior_tauv = transpose(rebin(prior_tauv, 24,33,9,24), [1, 2, 0,3])
      endelse
      p = p*prior_tauv
      p = p/total(p)
   endif

   

   ;;
   ;; We can also have a prior on log mu - useful for
   ;; variable N calculations where it might not be easy to
   ;; constrain. 
   ;;
   if (keyword_set(apply_mu_prior)) then begin
      if (not keyword_set(mu_dependent_n)) then begin
         print, 'GENERIC_CHANGER_OF_FITS: There is no point setting APPLY_MU_PRIOR'
         print, 'unless you also have VARIABLE_N set. Not proceeding!'
         return
      endif
      
      ;; Create a prior on mu but 1s everywhere else.
      prior = p*0.0+1
      if (size(apply_mu_prior, /tname) eq 'STRUCT') then begin
         ;; The user has provided the settings for the xsi prior. 
         mean_mu = struct_var(apply_mu_prior, 'MEAN_XSI')
         sigma_mu = struct_var(apply_mu_prior, 'SIGMA_XSI')
      endif else begin
         ;; This is the mean and 1 sigma scatter from the NG +
         ;; NGC628 sample.
         mean_mu = 0.66
         sigma_mu = 0.1 
      endelse
      
      ;; Define the prior as a Gaussian.
      prior_mu = exp(-(10.0^dim.logmu-mean_mu)^2/(2.0*sigma_mu^2))/(2*sqrt(!pi*sigma_mu))
;         stop
      prior_mu = transpose(rebin(prior_mu, 8, 33, 9, 24, 24), [1, 2, 3, 0, 4])
;         stop
      p = p*prior_mu
      p = p/total(p)
   endif
   
   min_chisq = min(c.chi2)


   ;;---------------------------------------------
   ;; Marginalise the probability distributions
   ;;---------------------------------------------

   if (keyword_set(mu_dependent_n)) then begin
      tmp_Z = total(total(total(total(p, 1), 1), 1), 1) 
      bin_z=tmp_Z/total(tmp_Z)
      tmp_U = total(total(total(total(p, 5), 4), 3), 2)
      bin_U=tmp_U/total(tmp_U)
      tmp_xsi = total(total(total(total(p, 5), 4), 3), 1) 
      bin_xsi=tmp_xsi/total(tmp_xsi)
      tmp_tauv = total(total(total(total(p, 5), 4), 2), 1) 
      bin_tauv=tmp_tauv/total(tmp_tauv)
      tmp_logmu = total(total(total(total(p, 5), 3), 2), 1) 
      bin_logmu=tmp_logmu/total(tmp_logmu)
   endif else begin
      tmp_Z = total(total(total(p, 1), 1), 1) & bin_Z=tmp_Z/total(tmp_Z)
      tmp_U = total(total(total(p, 4), 3), 2) & bin_U=tmp_U/total(tmp_U)
      tmp_xsi = total(total(total(p, 4), 3), 1) & bin_xsi=tmp_xsi/total(tmp_xsi)
      tmp_tauv = total(total(total(p, 4), 1), 1) & bin_tauv=tmp_tauv/total(tmp_tauv)
   endelse

   ;; 12 + Log O/H
   bin_oh = whist1d(yoh[*], p[*], binsize=0.1, min=7.5, max=9.5, $
                    obin=obin, omax=omax, omin=omin, density=density)
   xoh = obin

   ;; Log Eta(Ha)
   bin_etaha = whist1d(etaha[*], p[*], binsize=0.1, min=5.5, max=7.9, $
               obin=obin_eha, omax=omax, omin=omin, density=density)
   xetaha = obin

   ;; Log Eta([O II])
   bin_etao2 = whist1d(etao2[*], p[*], binsize=0.1, min=5.5, max=7.9, $
               obin=obin_eo2, omax=omax, omin=omin, density=density)
   xetao2 = obin_eo2

   ;; Log Mugas
   bin_mugas = whist1d(mugas[*], p[*], binsize=0.1, min=-1, max=3.8, $
                       obin=obin, omax=omax, omin=omin, density=density)
   xmugas = obin

;   stop
   
   ;; And T_e
   bin_te = whist1d(te[*], p[*], binsize=150, min=1750, max=18000, $
               obin=obin_te, omax=omax, omin=omin, density=density)
   xte = obin_te


   ;; The area input is assumed to be in kpc^2.
   if (n_elements(area_in) eq 0) then begin
      ;; If not provided, assume the SDSS fibre size.
      area = !dpi*arcsec_to_kpc(1.5, z, h=h, omega=omega0, $
                                lam=lambda0)^2
   endif else $
    area = area_in[ind]



   ;;---------------------------------------------------------
   ;; Was a stellar mass structure provided. If so, let us
   ;; calculate RGAS and FGAS
   ;;---------------------------------------------------------
   if (n_elements(lgm) gt 0) then begin
      if (tag_exist(lgm, 'LOGMU')) then $
       log_mustar = lgm.logmu[ind] $
      else begin
         log_mustar = lgm.lgm[ind] - alog10(area) ; Msun/kpc^2
      endelse
      log_mustar = log_mustar-6 ; pc^-2

      rgas = mugas[*]-log_mustar
      bin_rgas = whist1d(rgas, p[*], binsize=0.05, $
                         min=-3, max=2, $
                         obin=obin_rgas, omax=omax, omin=omin, density=density)
      xrgas = obin_rgas

      fgas = mugas[*]-alog10(10.0^mugas[*]+10.0^log_mustar)
      bin_fgas = whist1d(fgas, p[*], binsize=0.05, $
                         min=-3, max=0.1, $
                         obin=obin_fgas, omax=omax, omin=omin, density=density)
      xfgas = obin_fgas

;      stop
      
   endif
   

   ;;------------------------------------------------------------
   ;; SFR - this assumes that the input flux is in luminosities!
   ;;------------------------------------------------------------
   i_ha = where(lines eq 'LHa_6563', n_ha) 

   i_ha = i_ha[0]
   l_ha = data[i_ha]            ; I am assuming luminosity input

;   stop
   
   ;; SFR from combination with EtaHa.
;      log_SFR = alog10(l_ha) - etaha
   ;; And this is using the normalisation - this is what the
   ;; model fitting normally uses.
   log_SFR = alog10(c.a[*])

   bin_sfr = whist1d(log_SFR, p[*], binsize=0.05, min=-5, max=4., $
                     obin=obin_sfr, omax=omax, omin=omin, density=density)
   xsfr = obin_sfr
   
   ;; And mugas/SFR.
   
   mu_sfr = (log_sfr - alog10(area)) ; log Msun/yr/kpc^2
      

   log_tR = mugas - (mu_sfr - 6) ; Convert area to pc^2 
   
   bin_tr = whist1d(log_tR[*], p[*], binsize=0.075, min=6.5, max=11.5, $
                    obin=obin_tr, omax=omax, omin=omin, density=density)
   xtr = obin_tr




   ;;------------------------------
   ;; Multi-D PDFs. If requested.
   ;;------------------------------
   if (keyword_set(ks_plane)) then begin
      ;; 
      ;; The user wants MuGas vs MuSFR. - The KS law basically.
      ;;
      
      use = where(mugas gt -2 and mugas lt 3 and mu_sfr gt -3 and mu_sfr lt 2)
      
      if (keyword_set(fast_planes)) then begin
         xx = [[mugas[use]], [mu_sfr[use]]]
         Cw = weighted_covariance(xx, p[use])            
         xpos = total(mugas[use]*p[use])/total(p[use])
         ypos = total(mu_sfr[use]*p[use])/total(p[use])
         
         ;; Now create a structure like that returned by
         ;; JB_FIT_2DGAUSS:
         invC = invert(Cw)

         ;; Arbitrary peak value because we do not need that. 
         res = [1.0, xpos, ypos, invC[0,0], invC[1,0], invC[0, 1], invC[1,1]]
         ks_res = {parameters: res, $
                   x_mean: res[1], $
                   y_mean: res[2], $
                   precision: invC, $
                   covariance: Cw, $
                   amplitude: res[0], $
                   niter: 0}
         
;            stop
         
      endif else begin
         ;; Calculate 2D fits to the data.
         whist2d, mugas, mu_sfr, h_ks, [-1, 3.0], [-3, 1], 100, 100, weight=p[*], /linbin
         
         ;; To avoid rounding errors I scale p here to a max of 1d9
         jb_fit_2dgauss, mugas[use], mu_sfr[use], p[use], weight=p[use], $
                         res=ks_res, yfit=yfit, c_code=c_code
         
         print, 'NITER=', ks_res.niter
      endelse
   endif
   
   if (keyword_set(ksoh_plane)) then begin
      ;;
      ;; For the fit I enforce a min(chi2)=1. Otherwise the
      ;; fit doesn't generally converge.
      ;;
      use = where(mugas gt -2 and mugas lt 3 and yoh gt 7)
      
      if (keyword_set(fast_planes)) then begin
         xx = [[mugas[use]], [yoh[use]]]
         Cw = weighted_covariance(xx, p[use])            
         xpos = total(mugas[use]*p[use])/total(p[use])
         ypos = total(yoh[use]*p[use])/total(p[use])
         ;; Now create a structure like that returned by
         ;; JB_FIT_2DGAUSS:
         invC = invert(Cw)
         
         ;; Arbitrary peak value because we do not need that. 
         res = [1.0, xpos, ypos, invC[0,0], invC[1,0], invC[0, 1], invC[1,1]]
         ohp_res = {parameters: res, $
                    x_mean: res[1], $
                    y_mean: res[2], $
                    precision: invC, $
                    covariance: Cw, $
                    amplitude: res[0], $
                    niter: 0}
      endif else begin
         jb_fit_2dgauss, mugas[use], yoh[use], p[use], weight=p[use]^(1.0/min_chisq), $
                         res=ohp_res, yfit=yfit, c_code=c_code
;         h = jb_2dgauss_image(ohp_res.parameters, {nx: 512, ny: 512, xrange: [-1, 3], yrange:[8.5, 9.5]})
;         stop
      endelse
   endif

   if (keyword_set(kss_plane)) then begin
      ;; Specific SFR within the fibre
      ssfr = log_sfr - lgm.lgm[ind]
      
      use = where(mugas gt -2 and mugas lt 3 and ssfr gt -13 and ssfr lt -7, n_use)
      if (n_use gt 0) then begin
         if (keyword_set(fast_planes)) then begin
            xx = [[mugas[use]], [ssfr[use]]]
            Cw = weighted_covariance(xx, p[use])            
            xpos = total(mugas[use]*p[use])/total(p[use])
            ypos = total(ssfr[use]*p[use])/total(p[use])
            
            ;; Now create a structure like that returned by
            ;; JB_FIT_2DGAUSS:
            invC = invert(Cw)
            
            ;; Arbitrary peak value because we do not need that. 
            res = [1.0, xpos, ypos, invC[0,0], invC[1,0], invC[0, 1], invC[1,1]]
            kss_res = {parameters: res, $
                       x_mean: res[1], $
                       y_mean: res[2], $
                       precision: invC, $
                       covariance: Cw, $
                       amplitude: res[0], $
                       niter: 0}
               
         endif else begin
            jb_fit_2dgauss, mugas[use], ssfr[use], p[use], weight=p[use], $
                            res=kss_res, yfit=yfit, c_code=c_code
         endelse
      endif else kss_res = {n_use: 0}
   endif
   


   ;; And Dust-to-gas ratio
   bin_dgr = whist1d(dgr[*], p[*], binsize=0.05, min=-3.7, max=-1.375, $
                     obin=obin_dgr, omax=omax, omin=omin, density=density)
   xdgr = obin_dgr

   ;;
   ;; In the case of a mu-dependent n we also record PDFs for n. 
   ;;
   if (keyword_set(mu_dependent_n)) then begin
      bin_n = whist1d(n_slope, p, binsize=0.025, min=0.4, max=1.4, $
                      obin=obin_n, omax=omax, omin=omin, density=density)
         
      xn = obin_n
   endif

   results = {bin_z: bin_z, bin_u: bin_u, bin_xsi: bin_xsi, bin_tauv: bin_tauv, $
              bin_oh: bin_oh, xoh: xoh, min_chisq: min_chisq, $
              bin_etaha: bin_etaha, xetaha: xetaha, bin_etao2: bin_etao2, $
              xetao2: xetao2, xmugas: xmugas, bin_mugas: bin_mugas, $
              xdusttogas: xdgr, bin_dusttogas: bin_dgr, bin_te: bin_te, xte: xte}

   if (n_elements(bin_tr) gt 0) then begin
      results = create_Struct(results, 'bin_logSFR', bin_SFR, 'xlogSFR', xsfr, $
                              'bin_log_tR', bin_tr, 'xlog_tr', xtr)
   endif

   if (n_elements(lgm) gt 0) then begin
      results = create_struct(results, 'bin_rgas', bin_rgas, 'xrgas', xrgas, $
                              'bin_fgas', bin_fgas, 'xfgas', xfgas)
   endif


   if (keyword_set(ks_plane)) then $
    results = create_struct(results, 'ksp', ks_res)
   if (keyword_set(ksoh_plane)) then $
    results = create_struct(results, 'ohp_minchi1', ohp_res)
   if (keyword_set(kss_plane)) then $
    results = create_struct(results, 'kssp', kss_res)


   if (keyword_set(mu_dependent_n)) then begin
      ;; And add some more output if a mu-dependent n is assumed. 
      results = create_struct(results, 'bin_logmu', bin_logmu, 'bin_n', bin_n, $
                              'xn', xn)
   endif

   ;; Finally add the entropy and the Kullback-Leibner divergence.
   todo = ['Z', 'U', 'XSI', 'TAUV', 'OH', 'ETAHA', 'ETAO2', 'MUGAS', 'DUSTTOGAS']
   if (keyword_set(mu_dependent_n)) then $
    todo = [todo, 'LOGMU', 'N']
   
   for i=0L, n_elements(todo)-1 do begin
      y = struct_var(results, 'BIN_'+todo[i])
      entropy = entropy_of_pdf(y)
      results = create_struct(results, 'ENTROPY_'+todo[i], entropy)
   endfor
      
;   stop
end
