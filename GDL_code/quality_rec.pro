pro quality_rec, cuberec, chi2	


;readcol, pathrec+listrec, recnames, format = 'A'
n = n_elements(cuberec(0,0,*))

rad = 10
chi2 = dblarr(n)
sig = dblarr(n)

;rec = readfits(pathrec+recnames[0])
n1rec = n_elements(cuberec[*,0,0])
coord, n1rec, xrec, yrec
rrec = sqrt((xrec-n1rec/2)^2+(yrec-n1rec/2)^2)

xselrec = xrec(where(rrec lt rad))
yselrec = yrec(where(rrec lt rad))

maskrec = dblarr(n1rec, n1rec)
maskrec(xselrec, yselrec) = 1

for k = 0, n-1 do begin
    ;name = recnames[k]
    ;rec = readfits(pathrec+name)
   rec = cuberec(*,*,k)

   noise = [rec(0,*), rec(n1rec-1,*), transpose(rec(*,0)), transpose(rec(*,n1rec-1)), rec(1,*), rec(n1rec-2,*), transpose(rec(*,1)), transpose(rec(*,n1rec-2))]

   sig = stddev(noise)

   chi2[k] = total(((rec*maskrec)/sig)^2)/n_elements(xselrec)
    
 endfor


end
