;+
; NAME: 
;      build_ring
;
; PURPOSE: 
;      
;     Reconstructs images of galaxies from a PCA basis and a set of components
;     (These sets and basis have to be found as fits file respectively
;     called 'Eigen_vec_rmring.fits' and 'Eigen_val_rmring.fits') with
;     a limited number of coefficients. The reconstructed images are
;     then subtracted to original images. The resulting images are
;     thus free from their central galaxy. These residuals are written
;     as a cube in a fitsfile
;
; CALLING:
;
;     build_ring, path, fitsfilein, fitsfileout, Training_set =
;     training_set, T_set = T_set, nvecini, base_cleaning =
;     base_cleaning, qfactor 
;
; INPUTS: 
;     
;     PATH: path where to look for PCA basis and coeefficients as well
;     as original images for comparison. Residuals will also be
;     written at that location.
;
;     FITSFILEIN: name of the fitsfile containing a cube of original
;     images for which the PCA basis has been built. They are also the
;     original images from which central galaxies will be removed.
;          
;     FITSFILEOUT: name of the fitsfile where to write the cube of
;     images resulting from subtraction.
;
;     NVECINI: number of PCA coefficients to be used when
;     reconstructed the images. This number will highly impact the
;     quality of the galaxy subtraction. A higher number will result
;     in "overreconstructing" the galaxies by reconstructing noise
;     and/or galaxies' neighbous that we want to keep unchanged.
;
; KEYWORDS:
;
;     \TRAINING_SET: If set, FITSFILEIN will be ignored and T_SET is
;     expected to be provided as a cube of images that will be
;     decomposed on the PCA basis and used for subtraction.
;
;     T_SET: if /TRAINING_SET is set T_SET is expected to provide a
;     set of images with central galaxies to be subtracted on an
;     already existing PCA basis.
; 
;     /BASE_CLEANING: if set, a second iteration of PCA will be
;     performed on cleaner images of galaxies. Indeed, once central
;     galaxies have been reconstructed and subtracted, the resulting
;     residuals are smoothed (by convolving the images with a gaussian
;     profile) and subtracted from the original images living only
;     the central galaxies intact. These newly formed images are used
;     to form a new basis for PCA analysis. Then the
;     reconstruction/subtraction process is repeated producing images
;     free from their central galaxy with a limited number of coefficients.     
;
; OUTPUTS: 
;     
;      QFACTOR: quality factor for reconstruction. This factor is the
;      median value of chi squarre between original and reconstructed
;      images computed in the central region of images. Values between
;      1 and 1.2 occure to provide good reconstructions.
;      
; EXAMPLE: 
;    
;      build_ring, './', 'Galaxies_3_4_20.fits', 'Residuals_3_4_20_q=1.03.fits', 50, 
;      /base_cleaning, qfactor
;
; HISTORY:
;     
;-
;--------------------------------------------------------------------------
pro build_ring, path, fitsfilein, fitsfileout, Training_set = training_set, T_set = T_set, nvecini, base_cleaning = base_cleaning, qfactor

E_N2 = readfits(path+'Eigen_vec_rmring.fits')
alpha = readfits(path+'Eigen_val_rmring.fits')
print, 'eigen vectors and values loaded'
;;;;;# The original images are cut in patches of the same size as the
;;;;;reconstructed images
n = n_elements(E_N2(*,0))
print, n
i=0

if keyword_set(Training_set) then begin 
   S_set = dblarr(sqrt(n)*sqrt(n), n_elements(T_set(0,0,*)))
   S_set(*,*) = reform(T_set)
   alpha = matrix_multiply(E_N2,S_set, /ATRANSPOSE)
   galtest = T_set
endif else begin
   galtest = readfits(path+fitsfilein)
   print, 'original images for subtraction loaded'
endelse

;# rebuilds images with a limited number of coefficients
imvec = matrix_multiply(E_N2(*, 0:nvecini), alpha(0:nvecini, *))

npsf = n_elements(alpha(0,*))
imarray = dblarr(sqrt(n),sqrt(n), npsf)
;# If asked, the PCA is performed a second time on a set in which
;poluting companions have been removed
lisse = psf_gaussian(npixel = 10, fwhm = 0.5)
lisse = lisse/total(lisse)
nvecfinal = 1000000.
newset = imarray
if keyword_set(base_cleaning) then begin
  ; for i = 0, niter do begin; 
      nvecini = min([nvecini, nvecfinal])

      newset(*,*,*) = reform(imvec)
      imarray(*,*,*) = galtest-newset
      quality_rec, imarray, Q
      locup = where(Q lt 1.1)
      Q = Q(locup)
      locdown = where(Q gt 0.9)
      print, 'Size of the reduced basis :', n_elements(locdown)
      ;stop
      ;# creates a clean version of the set of images (no poluting companion)
      ;by subtracting the residuals once smoothed to the original image
      ;lisse = transpose(psf_gaussian(ndimen = 2, npixel = 10, fwhm = 0.5))
      ;
      ;# smoothing the residuals
      for k = 0, npsf-1 do imarray(*,*,k) = convol(imarray(*,*,k),lisse, /edge_truncate )
      ;# builds images free from poluting companions
      imres = galtest-imarray
      imvec(*,*)=reform(imres)
     
      ;# performs PCA on the set of galaxies without companions
      mk_PCA_ring, '', '', /Training_set, T_set = imvec(*,locdown), input_set = imvec, prefix = 'New_'+strmid(strcompress(i),1)+'_', /iset, /select

      if keyword_set(Training_set) then alpha = matrix_multiply(E_N2,S_set, /ATRANSPOSE)
      nvecfinal = min( n_elements(locdown)-1, nvecini)
      if nvecfinal ne nvecini then print, 'Number of coefficients in second iteration can not be as high as in the first iteration'
      print, 'number of PCA coefficients in the last iteration: ', nvecfinal
      imvec = matrix_multiply(E_N2(*, 0:nvecfinal), alpha(0:nvecfinal, *))
  ;    ;stop
  ; endfor
endif

;# refoms images from the set of reconstructed vectors
imarray(*,*,*) = reform(imvec)

;# Final subtraction. Produces clean images without central galaxies
residuals = galtest-imarray

quality_rec, residuals, q
print, 'quality factor =', median(q)

qfactor = median(q)

writefits, fitsfileout, residuals;, header

end
