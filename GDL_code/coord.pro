pro coord, size, xin, yin

  xin = dblarr(size*size)
  yin = dblarr(size*size)
  for inc = 0,size-1 do begin
    xin(size*inc:(inc+1)*size-1) = inc
    yin(inc*size:(inc+1)*size-1) = indgen(size)
  endfor
end