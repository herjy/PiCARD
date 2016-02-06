;+
; NAME: 
;       select_gal
;
; PURPOSE: 
;       Selects galaxies from Field of Views or patches in order to
;       build a set of galaxies suitable for ring extraction by PCA
;       method.
;       It will rotate the galaxies so that each of them has the same
;       alignement and it writes the result in a .fits file
;
; CALLING:
;
;       select_gal, path, img=img, n1=n1, n2=n2, patchs = patchs
;
; INPUTS:
;      PATH : the path where to look for the FoV or the patches to be processed.
;
; KEYWORDS:
;    
;      IMG : Name of the field of view as it is in path in case of the
;      data are as a form of a field of view
;
;      N1, N2 : Sizes of the stamps to extract from the field of
;      view. These stamps must be at least sqrt(2)*NPCASTAMP,
;      NPCASTAMP being the size of the stanps used in the final
;      PCA. Indeed, rotation unallows us to take the full image we
;      extract from the FoV
;
;      PATCHS : If set, the program will look for a listfits.txt file
;      with the names of every stamps to use in the PCA.
;
; OUTPUTS:
;          S_BASE : a set of images in rows designed to be used as a
;          basis for the PCA
;
;          S_SET : Set of image to analyze. If SELECT is set, then
;          S_SET is different from S_BASE
; EXTERNAL CALLS:
;           None
;
; EXAMPLE:
;           
;
;
; HISTORY:
;	
;-


pro mk_set_ring, path, fitsfilein, S_base, S_set, gaussian = gaussian, fwhm = fwhm, select = select

;filename = 'listfits.txt'

;readcol, path+filename, patchname, format = 'A'
;#cube containing the images to use for the basis
galaxies = readfits(path+fitsfilein)
sz = n_elements(galaxies(0,0,*))
n1 = n_elements(galaxies(*,0,0))
n2 = n_elements(galaxies(0,*,0))

print, 'ok''performing PCA'
readcol, path+ 'select_gals.txt', num, Theta, x, y, a, b, format = 'A'

mult = dblarr(n1, n2, sz)
S_base = dblarr(n1*n2)
S_set = dblarr(n1*n2, sz)

;# multiplication by a gaussian profile if required
if keyword_set(gaussian) then begin
   gauss = psf_Gaussian( NPIXEL=[n1,n2], FWHM=fwhm )
   for p = 0, sz-1 do galaxies(*,*,p) = galaxies(*,*,p)*gauss
endif else fwhm = 20

S_set(*,*) = reform(galaxies) 

if keyword_set(select) then begin
   for k = 0, sz-1 do begin
      gal = galaxies(*,*,k)
      ;# Writes/overwrite a temporary file to be treated by sextractor
      writefits, path+'psf_out.fits', gal
      spawn, 'sex '+path+'psf_out.fits -c correct_temp.sex' ; psfout.fits remplace patchnames[k] pour essai
      readcol,'temp_object.txt', num, x, y, format = 'A'

      ;# Selects or rejects candidates for the basis
      is_suitable = 1
      if n_elements(num)-1 gt 0 then begin
         for j = 0, n_elements(num)-1 do begin
            r = ((n1/2-x[j])^2+(n2/2-y[j])^2)
            if r lt fwhm && gal(x[j], y[j]) gt gal[n1/2,n2/2]/2  then is_suitable = 0 
         endfor       
      endif
      if is_suitable eq 1 then S_base = [[S_base], [S_set(*,K)]]
   endfor
      S_base = S_base(*, 1:*)
      file_delete, path+'psf_out.fits'
endif else begin
   S_base = S_set
endelse

;for k = 0, sz-1 do begin
;   print, k
;   ;psf = readfits(path+patchname[k])
;   
;   szpatch = size(gal)
;   nn1 = (szpatch(1))
;   nn2 = (szpatch(2))
;   if not keyword_set(gauss) then begin
;      psfout = psf((nn1-n1)/2.:(nn1+n1)/2.-1, (nn2-n2)/2.:(nn2+n2)/2.-1)
;   endif else begin 
;      psfout = psf((nn1-n1)/2.:(nn1+n1)/2.-1, (nn2-n2)/2.:(nn2+n2)/2.-1)*gauss
      ;star = star2d(psf((nn1-n1)/2.:(nn1+n1)/2.-1, (nn2-n2)/2.:(nn2+n2)/2.-1))
      ;psfout = psfout+star(*,*,0)*(1-gauss)
;   endelse
   ;Set that is to be decomposed:
;   S_set(*,k) = reform(psfout(*,*))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;         Selection of apropriate galaxies for building the basis
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;writefits, path+'psf_out.fits', psfout

;;;;;;;;comment to try
;   if keyword_set(select) then begin
;      spawn, 'sex '+path+'psf_out.fits' ; psfout.fits remplace patchnames[k] pour essai
;      readcol,'test.txt', num, x, y, format = 'A'
;      file_delete, path+'psf_out.fits'
      ;r0 = ((x-n1/2)^2+(y-n2/2)^2); nn1 et nn2 remplace par n1 et n2 pour essai
     ;center = where(r0 eq min(r0))
;      is_suitable = 1
;      if n_elements(num)-1 gt 0 then begin
;         for j = 0, n_elements(num)-1 do begin
;            r = ((n1/2-x[j])^2+(n2/2-y[j])^2)
;            if r lt fwhm && psf(x[j], y[j]) gt psf[n1/2,n2/2]/2  then is_suitable = 0 
;         endfor       
;      endif
;      if is_suitable eq 1 then S_base = [[S_base], [S_set(*,K)]]
;   endif;

;   i+=1
;endfor

;S_set = S_set(*, 0:i-1)

;# Builds the final basis
if keyword_set(select) then begin
   S_base = S_base(*, 1:*)
endif else begin
   S_base = S_set
endelse

;print, size(S_base)
end
