;+
; NAME: 
;      mk_PCA_ring
;
; PURPOSE: 
;      Performs PCA on a set of galaxies. Galaxies with fewer
;      companions or parasites in their patch are selected to form the
;      basis. Images are reshaped into 1D vectors. The matrix formed
;      by the concatenation of these
;      vectors is used to formed a PCA Basis. Other galaxies are
;      decomposed on this basis. Both PCA basisand coefficients from
;      decomposition are written in .fits
;      files. Respectively:'Eigen_vec_rmring.fits' and 'Eigen_val_rmring.fits'
;
; CALLING:
;
;      mk_PCA_ring, path, fitsfilein,Training_set= training_set, T_set
;      =  T_set,E_N2, alpha, gauss = gauss, fwhm = fwhm, select =
;      select, iset = iset, input_set = S_base, prefix = prefix 
;
; INPUTS: 
;     
;     PATH: path where to look for images to build the PCA basis and
;     to write the output files
;
;     FITSFILEIN: name of the file containing a cube of images to
;     treat. 
;    
; KEYWORDS:
;     
;     /TRAINING_SET: if set, the program expects T_SET to be a set
;     images shapes as an array of vector (each vector being a
;     reshaped galaxy). This set will be used to form the basis. If
;     not set, mk_pca_ring will look for suitable candidates to
;     form the basis in the input file. 
;
;     T_SET: set to be used to build the PCA basis
;
;     /GAUSS: set this keyword to multiply omages used for the basis
;     by a gaussian profile (allows to clean the basis).
;
;     FWHM: fwhm of the gaussian to be used if /GAUSS is set.
;
;     /SELECT : if set, only galaxies with fewer
;     companions or parasites in their patch are selected to form the
;     basis.
;
;     /ISET: similar as /TRAINING_SET. In this case S_BASE is used
;     both as a basis for the PCA and a set to decompose.
;
;     INPUT_SET=S_BASE: if /ISET is set, S_base is the set to use to
;     perform PCA analysis.
;
;     PREFIX: adds a prefix to Eigen_val and Eigen_vecs files when writing
;
; OUTPUTS: 
;     E_N2 : vectorrs of the generated PCA basis.
;     ALPHA : coefficients of the decomposition of the vectors in
;      
; EXAMPLE: 
;     # For a classical use, only this is necessary: 
;     mk_PCA_ring, path'./', 'Galaxies_3_4_20.fits',E_N2, alpha, select =
;     select
;
; HISTORY:
;     
;-
;-------------------------------------------------------------------------------

pro mk_PCA_ring, path, fitsfilein,Training_set= training_set, T_set =  T_set,E_N2, alpha, gauss = gauss, fwhm = fwhm, select = select, iset = iset, input_set = S_base, prefix = prefix

init = SYSTIME( 1, /SECONDS )

;# Builds the images used to form the basis
if not keyword_set(iset) then begin
   print, 'making the set'
   mk_set_ring, path, fitsfilein, S_base, S_set, gauss = gauss, fwhm = fwhm, select = select
print, 'set built'
endif else begin
   S_set = S_base
   print, 'The set already exists'
endelse

if keyword_set(Training_set) then S_base = T_set


;# covaraince matrix
;cov = transpose(S_base)#(S_base)
cov = matrix_multiply(S_base, S_base, /ATRANSPOSE)
;# decomposes the covariance matrix to build the basis
la_svd, cov, W, U, V

print, 'svd done'

sz = size(W)

;# W are the eigenvalues of the covariance matrix
;#  W : eigenvalues 
;#  U : eigenvector

;if m gt n then begin
;# eigenvectors
  ;E_N2 = diag_matrix(double(1./sqrt(W)))##transpose(U)##(S_base)
E_N2 = matrix_multiply(matrix_multiply(S_base,U, /BTRANSPOSE), diag_matrix(double(1./sqrt(W))))
print, 'eigenvectors built'

;# eigenvalues
;alpha = (S_set)##transpose(E_N2)
alpha = matrix_multiply(E_N2,S_set, /ATRANSPOSE)

n1 = sqrt(n_elements(E_N2(*,0)))
n = n_elements(E_N2(0,*))
base = dblarr(n1, n1, n)
base(*,*,*) = reform(E_N2)

print, 'eigen values calculated'
if keyword_set(prefix) then begin
   writefits, path+prefix+'Base.fits', base
   writefits, path+prefix+'Eigen_vec_rmring.fits', E_N2
   writefits, path+prefix+'Eigen_val_rmring.fits', alpha 
endif else begin
   writefits, path+'New_Base.fits', base
   writefits, path+'Eigen_vec_rmring.fits', E_N2
   writefits, path+'Eigen_val_rmring.fits', alpha 
endelse



fin = SYSTIME( 1, /SECONDS )
print, 'total execution time: ',fin-init, ' seconds' 
end
