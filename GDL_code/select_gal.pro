;+
; NAME: 
;       select_gal
;
; PURPOSE: 
;       Extracts galaxies that are in a
;       specified range of sizes from fields of view or patches in order to
;       build a set of galaxies suitable for ring extraction by PCA
;       method.
;       Selected galaxies will be rotated so that each of them has the same
;       alignement. The resulting images are writen in a .fits file.
;       When runing this process on a Field for the first time, the
;       /havesex keyword is necessary for thus it will build the
;       catalog (field_objects.txt) of every object detected as a
;       Galaxy. Once this catalog is built, it is not necessary to run
;       it again even for other sizes of galaxies. 
;
; CALLING:
;
;       select_gal, path, fitsfileout, low, high, magmax, img=img, n1, n2, patches = patches, prefix = prefix, havesex = havesex, output_dir = output_dir, ngal.
;
; INPUTS:
;      PATH: the path where to look for the FoV or the patches to be
;      processed as well as the sextractor configuration files.
;      
;      FITSFILEOUT: string that will be used to name the .fits file
;      where the images of selected galaxies will be written.
;
;      LOW, HIGH: range of size in which ti select galaxies
;
;      MAGMAX: Maximal magnitude for a galaxy to be selected
;
;      N1, N2 = size of the windows to extract galaxies. The final
;      images will be 2^n*2^n images with (N1,N2) > sqrt(2)*(2^n,2^n)
;
;      
;
; KEYWORDS:
;    
;      IMG: If the galaxies to be treated are in a field of view, IMG
;      is a string that stands for the name of the field of view.
;
;      /PATCHES: If set, the program expects a set of patch images.  A
;      'listfits.txt' file is expected in PATH, which contains the
;      names of every stamps of galaxies to be processed.
;
;      /HAVESEX: if set (only for a field of view), have sextractor
;      running on IMG. It has to be ran only once per field. Even when
;      using select_gal to extract several different sets of galaxies
;      with different sizes, sextractor is not required as long as the
;      'field_objects.txt' file exists in PATH.
;
;      /HISTOPLOT: plots histogram of sizes and magnitudes as
;      calculated in sextractor.
;      
; OUTPUTS:
;
;          NGAL: number of galaxies found in the required range of sizes
;
;          Several files are written at the end of the execution: 
;          'field_objects.txt' contains the characteristics of each sextracted
;          objects.
;          'select_gals.txt' contains the characteristics of each
;          selected galaxies for the considered range of sizes.
;          FITSFILEOUT, a cube of images with the selected galaxies.
;
; EXTERNAL CALLS:
;           calls sextractor
;
; EXAMPLE:
;   # To extract two sets of galaxies with different sizes
;        
;   select_gal, './', 'Galaxies_3_4_20.fits', 3, 4, 20,
;   img='field.fits', 100,100, /havesex, ngal 
;
;   select_gal, './', 'Galaxies_4_6_20.fits', 4, 6, 20,
;   img='field.fits', 100,100, ngal 
;-


pro select_gal, path, fitsfileout, low, high, magmax, img=img, n1, n2, patches = patches, havesex = havesex, histoplot = histoplot, ngal

init = SYSTIME( 1, /SECONDS )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                               In case each galaxy is in a separated patch
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
n1temp = n1/sqrt(2)
n2temp = n2/sqrt(2)
n1t = 2^floor((alog(n1temp)/alog(2)))
n2t = 2^floor((alog(n2temp)/alog(2)))
galaxy = dblarr(n1t, n2t)

header = [' SIMPLE  =                    T / Written by IDL:  Thu Jul 11 17:19:55 2013',$
'BITPIX  =                  -64 /',$
'NAXIS   =                    2 /',$
'NAXIS1  =                   '+strmid(strcompress(n1t),1)+' /',$
'NAXIS2  =                   '+strmid(strcompress(n1t),1)+' /', $ 
'EXTEND  =                    T / File May Contain Extensions',    $
'CRVAL1  =          209.8892083 /',  $                                                
'CRVAL2  =          53.57527778 /',  $                                                
'CRPIX1  =          9677.500000 /',  $                                                
'CRPIX2  =          9677.500000 /',  $                                                
'CD1_1   =     -0.0000516666679 /',  $                                                
'CD1_2   =          0.000000000 /',  $                                                
'CD2_1   =          0.000000000 /',  $                                                
'CD2_2   =      0.0000516666679 /',  $
'END']  

if keyword_set(patches) then begin
if not keyword_set(prefix) then prefix = 'GAL'

readcol, path+'listfits.txt', names, format = 'A'
n = long(n_elements(names))

if keyword_set(histoplot)then,
 plothist, A, gaussfita, /noplot
 plothist, B, gaussfitb, /noplot
;stop
;# distribution of sizes and magnitudes
 window, 0
 plothist, A,  bin = 0.2, xrange = [0,25]
 plothist, B,  bin = 0.2, color = 120, /overplot
 plothist, mag, bin = 0.1, color = 60, /overplot
endif

;recnames = ['']
angles = [0]
select = [0]
xs = [0]
ys = [0]
As = [0]
Bs = [0]
For k = 0L, n-1 do begin
   print, k
   ;# execute sextractor on each patch
   spawn, 'sex '+path+names[k]+' -c correct_temp.sex'
   readcol, './temp_object.txt', num, x, y, A, B, Theta, format = 'A'
   gal = readfits(path+names[k])
   n0 = n_elements(gal(0,*))
   ;# select the wrigth object in the patch
   if n_elements(num)gt 1 then begin
      r = ((x-n0/2)^2+(y-n0/2)^2)
      minr = where(r eq min(r))
      theta0 = Theta[minr]
      B = B[minr]
      A = A[minr]
      endif else begin
      theta0 = Theta[0]
      num0 = num[0]
      B = B[0]
      A = A[0]
   endelse
   print, A, B
   ;if A lt high && A gt low && B lt high && B gt low then begin
      if B gt A then begin 
         galrot = rot(gal, theta0[0]+90)
      endif else begin
         galrot = rot(gal, theta0[0])
      endelse
      galaxy = [[[galaxy]],[[galrot[(n1-n1t)/2:(n1+n1t)/2-1, (n1-n1t)/2:(n1+n1t)/2-1]]]]
      ;writefits, path+'Galaxy_set/'+names[k], galrot, header
      select = append(select,k)
      angles = append(angles, -theta0)
      xs = append(xs, x[minr])
      ys = append(ys, y[minr])
      As = append(As, A)
      Bs = append(Bs, B)
      ;recnames = append(recnames,[names[k]])
   ;endif
endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;                                       in case the data are a FoV
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif else begin
;# executes sextractor over the FoV

if keyword_set(havesex) then spawn, 'sex '+path+img+' -c field.sex'

inter = SYSTIME( 1, /SECONDS )
print, 'execution time after sextractor (or not): ',inter-init, ' seconds' 

if not keyword_set(patches) then $
readcol, path+'field_objects.txt', num, mag, x, y, alpha,delta, A, B, Theta, format = 'A'
if keyword_set(patches) then readcol, path+'field_objects.txt', num, mag, x, y, A, B, Theta, format = 'A'

n = long(n_elements(num))
field = readfits(img)
sz = size(field)

As = [0]
Bs = [0]
mags = [0]
angles = [0]
select = [0]
xs = [0]
ys = [0]
alphas = [0]
deltas = [0]

;#chooses the galaxies that are far enough from edges in order to form
;an homogeneous set of galaxies in terms of size for the pca 
for i = 0L, n-1 do begin
   borders = [n2/2, sz[2]-y[i],y[i], x[i], sz[1]-x[i], n1/2]
   ;bordermin = n1/2.; 
   bordermin = min(abs(borders))
   galtemp = dblarr(n1,n2)
   ;# Cuts patches around the found galaxies that matches the size criterion
   galtemp = field(x[i]-bordermin:x[i]+bordermin-1,y[i]-bordermin:y[i]+bordermin-1)
   spot = where(galtemp ge 5)
   nspots = n_elements(spot)
   if A(i) lt high && A(i) gt low && B[i] gt A[i]/2 && mag[i] lt magmax && nspots lt 2 then begin
                                ;# Runs sextractor on small patches to
                                ;find the exact location of the
                                ;centroid
      writefits, 'Galtemp.fits', galtemp
      spawn, 'sex Galtemp.fits -c correct_temp.sex'
      readcol, 'temp_object.txt', numtemp, xtemp, ytemp, format = 'A'
      numtemp = uint(numtemp)
      ;if n_elements(numtemp gt 0) then begin
         xdiff_tab = n1/2.-xtemp
         ydiff_tab = n2/2.-ytemp
         xdiff = min([abs(xdiff_tab)])
         ydiff = min([abs(ydiff_tab)])
         x[i] = x[i]-xdiff*mk_sign(xdiff_tab(where(abs(xdiff_tab) eq xdiff)))
         y[i] = y[i]-ydiff*mk_sign(ydiff_tab(where(abs(ydiff_tab) eq ydiff)))
         gal = field(x[i]-bordermin:x[i]+bordermin-1,y[i]-bordermin:y[i]+bordermin-1)

;#rotates the galaxies so that every patch will have galaxies with 
;the same inclination
         if B[i] gt A[i] then begin 
            galrot = rot(gal, Theta[i]+90.,/INTERP)
         endif else begin
            galrot = rot(gal, Theta[i],/INTERP)
         endelse
         galaxy = [[[galaxy]],[[galrot[(n1-n1t)/2:(n1+n1t)/2-1,(n1-n1t)/2:(n1+n1t)/2-1]]]] 
         select = append(select, num[i])
         angles = append(angles, theta[i])
         xs = append(xs, x[i])
         ys = append(ys, y[i])
         As = append(As, A[i])
         Bs = append(Bs, B[i])
         alphas = append(alphas, alpha[i])
         deltas = append(deltas, delta[i])
         mags = append(mags, mag[i])
         ;plt_image, galrot
         print, i
         ;endif
      ;endif

   endif
endfor

endelse

if keyword_set(patches) then begin
   infos = [[select[1:*]], [angles[1:*]],[xs[1:*]],[ys[1:*]], [As[1:*]],[Bs[1:*]]]
endif else begin
   infos = [[select[1:*]], [angles[1:*]],[xs[1:*]],[ys[1:*]], [As[1:*]],[Bs[1:*]], [mags[1:*]], [alphas[1:*]],[deltas[1:*]]]
endelse
print, 'encore plus gros zizi'

writefits, PATH+fitsfileout, galaxy(*,*,1:*)

;#writes the properties of the selected galaxies in a file : position,
;number, original alignement....

OPENW,1,PATH +'select_gals.txt' 
PRINTF,1,transpose(infos)
CLOSE,1

;file_delete, './listfits.txt'
;OPENW,2,path+'./listfits.txt' 
;PRINTF,2,transpose(recnames)
;CLOSE,2

print, 'The actual size of the patches is now ', n1t,'x',n2t
fin = SYSTIME( 1, /SECONDS )
print, 'execution time after sextractor (or not): ',inter-init, ' seconds' 
print, 'total execution time: ',fin-init, ' seconds' 

print, 'Selectiom done'
print, n_elements(select)-1, ' galaxies have been extracted.'

ngal =  n_elements(select)-1
end
