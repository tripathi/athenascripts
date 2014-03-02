FUNCTION alphab, T
a = 2.59e-13*(T/1.0d4)^(-0.7)
RETURN,a
END

PRO PLOTCOMP, ps=ps, plotvel=plotvel, plottemp=plottemp, plotpres=plotpres, plotdens = plotdens, plotni = plotni, plotifrac = plotifrac, plotni2=plotni2, plotmass=plotmass
;Plotting related to ionization fronts for ionized HII regions
;Updated 2/27/14 - for mass plotting


;Plotting setup
plotsym, 0, .75, /fill
set_plot, 'x'
if (keyword_set(ps)) then begin
    !p.font=0
    !p.charsize=1.15
    !p.charthick=4
    !p.thick=5
    !x.thick=5
    !y.thick=5
;    !x.margin = [2,2]
;    !y.margin = [2,2]
endif

;Postscript setup
;Gas velocity plot
if (keyword_set(plotvel)) then begin
   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='velplot.ps', xsize=8., ysize=10.5, xoffset =.13, yoffset = .13, /inches, /portrait
      !p.multi = [0,5,5]
   endif
endif
;Temperature plot
if (keyword_set(plottemp)) then begin
   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='tempplot.ps', xsize=8., ysize=10.5, xoffset =.13, yoffset = .13, /inches, /portrait
      !p.multi = [0,5,5]
   endif
endif
;Pressure plot
if (keyword_set(plotpres)) then begin
   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='pressureplot.ps', xsize=8., ysize=10.5, xoffset =.2, yoffset = .13, /inches, /portrait
      !x.charsize = 1
      !p.multi = [0,3,5]
   endif
endif
;Density plot
if (keyword_set(plotdens)) then begin
   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='dens.ps', xsize=8., ysize=10.5, xoffset =0.13, yoffset = 0.13, /inches, /portrait
      !p.multi = [0,2,3]
   endif
endif
;Ionization fraction plot
if (keyword_set(plotifrac)) then begin
   if (keyword_set(ps)) then begin
      set_plot, 'ps'
      device, filen='ifrac.ps', xsize=8., ysize=10.5, xoffset =.13, yoffset = .13, /inches, /portrait
      !p.multi = [0,5,5]
   endif
endif

if (keyword_set(plotni2)) then !p.multi=[0,3,4]


;Physical values
kb = 1.38d-16 ;Boltzmann's constant in cgs
n0 = 63.       ;Neutral number density [cm^-3]
flux = 3.67d9 ;Photoionizing number flux
alphac = 3.0d-3                 ;Carbon
m_H = 2.3d-24                   ;Mean mass per hydrogen nucleus [g]
mu_n = 2.1d-24                  ;Mean mass per particle of neutral atomic gas [g]

mu_i = m_H/2.                   ;Mean mass per particle of ionized atomic gas [g] - Should REALLY be mu_n/2.
;; alpha = 2.59e-13*(T/1.0d4)^(-0.7) ;Case-B recombination coefficient



;Analytic calculations
lstrom = flux/alphab(1.d4)/n0/n0 ;Stromgren length for 10^4K gas


;Read in all binary files in directory
list=findfile('*.bin')
nfiles = n_elements(list)

;Initialize variables
times = dblarr(nfiles)
dens = dblarr(nfiles)
pos = dblarr(nfiles)
xend = dblarr(nfiles)
cs2 = dblarr(nfiles)
testing = dblarr(nfiles)
temperature = dblarr(nfiles)
niprod = dblarr(nfiles)
maxpy = dblarr(nfiles)
maxpz = dblarr(nfiles)
ni_sim = dblarr(nfiles)
massini = dblarr(nfiles)
massnow = dblarr(nfiles)
massnow2 = dblarr(nfiles)

checkedsize = 0

for i = 0, nfiles-1 do begin
    ;Clear variables
    x = 0 & y = 0 & z = 0
    px = 0 & py = 0 & pz = 0
    d = 0 & energy = 0 & esp = 0
    
    ;Use Mark's script to read in files
    read_bin_ath, list[i], d, px, py, pz, ion=1, e=energy, time = t, xvec = x, yvec=y, zvec=z, d_n = dneutral


    if (checkedsize lt 1) then begin
       griddim = size(d)
       ncells = griddim[2]            ;Number of cells in box
       checkedsize = 999
    endif
 
dnan = finite(d, /nan)
print, where(dnan)

    vel = px/d                  ;Gas x-velocity
    xpc = x/(3.086d18)          ;x-distance in pc

    ;Checking that other directions have 0 momentum
    maxpy[i] = max(py)
    maxpz[i] = max(pz)

    ;Calculate ionization fraction
    nh = dneutral/m_h           ;Number density of neutral H (atomic)
    nhp = (d - dneutral)/m_h    ;Number density of ionized H
    nel = nhp + d * alphac / (14. * m_h) ;Electron number density
    ifrac = nel/(nh + nhp)               ;Ionization fraction

   ;Calculate T the way Mark does in ionrad_3d.c
    eperm = (energy - .5 * (px*px+py*py+pz*pz)/d)/d ;Thermal energy per unit mass
    mu_temp = (ifrac * mu_i) + (1.-ifrac)*mu_n
;print, mean(mu_temp), mu_i, mu_n
    temp = eperm * mu_temp / kb ;Temperature of gas, based on ionization fraction [K]

    ;Sound speed, Pressure, and Number Density
    cs2all = eperm; kb *temp / mu
    pressure = d * cs2all
    nall = d/mu_temp

    ;Set overall timestep values
    times[i] =t
    dens[i] = max(d[*,0,0], maxdloc)
    pos[i] = x[maxdloc] - x[0]
    massini[i] = m_H * n0 * ((x[maxdloc]-x[0])-lstrom) ;Not sure why lstrom was included here

    ;Find average sound speed and number density in ionized region
    if (maxdloc le 5) then begin
        cs2[i] = cs2all[floor(maxdloc/2),0,0]
        ni_sim[i] = nall[floor(maxdloc/2),0,0]
        temperature[i] = temp[floor(maxdloc/2),0,0]
    endif else begin
        csioniz = cs2all[0:maxdloc-5,0,0]
        nioniz = nall[0:maxdloc-5,0,0]
        tempioniz = temp[0:maxdloc-5,0,0]
        cs2[i] = mean(csioniz[sort(csioniz)])
        ni_sim[i] = mean(nioniz)
        temperature[i]=mean(tempioniz)
    endelse

;if (maxdloc eq 0) then begin
;    stop
;    cs2[i]=cs2all[0]
;endif else begin
;    csioniz = cs2all[0:maxdloc]
;    niprodorig=d[0:maxdloc]/mu *sqrt(x[0:maxdloc])
;    cs2[i] = median(csioniz[sort(csioniz)], /even)
;    niprod[i] = median(niprodorig(sort(niprodorig)), /even)
;endelse

    ;If the front is still in the box and the box has 64 cells:
if (keyword_set(plotmass)) then begin

;Find edge of ionization front
subx = d[maxdloc+1:*,0,0]
minsubx = where (d[*,0,0] eq min(subx))
xend[i] = min(minsubx)



    if (maxdloc lt ncells-1) then begin
       if (maxdloc eq 0) then begin
          massnow[i] = 0
       endif else begin
          massnow[i] = int_tabulated(x[maxdloc:xend[i]-1]-x[0], d[maxdloc:xend[i]-1,0,0])
          massnow2[i] = mean(d[maxdloc:xend[i]-1,0,0]) * abs(x[xend[i]-1]-x[maxdloc])
       endelse
    endif
 endif

    ;Make time step plots

    ;Ionization fraction plot
    if (keyword_set(plotifrac)) then begin
        if (t lt 1.35d15) then begin
            plot, xpc, ifrac, xsty = 1, ysty = 2, title = t, /nodata, xtitle = 'X [pc]', ytitle= 'Ionization fraction'
            oplot, [-50,50], [0,0], linestyle =2, color=fsc_color("gray")
            oplot, [-50,50], [1,1], linestyle =2, color=fsc_color("gray")
            oplot, xpc, ifrac, color=fsc_color("navy")
        endif
     endif

    ;Gas velocity plot
    if (keyword_set(plotvel)) then begin
        !x.omargin=[0,0]
        !y.omargin=[0,0]
        if (i mod 32 eq 0) then begin
            plot, x, vel, title = t, xsty=1, ysty=1, /nodata
            oplot, x, vel, color=fsc_color("blue")
        endif
    endif

    ;Pressure plot
    if (keyword_set(plotpres)) then begin
        !x.omargin=[0,0]
        !y.omargin=[0,0]
        if (i mod 15 eq 0 AND t lt 1.3d15) then begin
           if (i eq 0) then begin
              plot, xpc, pressure, title = t,  /nodata, xsty=1, yra=[5e-13,5e-10], /ylog, xtitle = 'X [pc]', ytitle ='Pressure [microbar]'
           endif else begin
              plot, xpc, pressure, title = t,  /nodata, xsty=1, yra=[5e-13,5e-10], /ylog, xtitle = 'X [pc]'
           endelse
           oplot, xpc, pressure, color=fsc_color("green")
        endif
     endif

    ;Temperature
    if (keyword_set(plottemp)) then begin
        !x.omargin=[0,0]
        !y.omargin=[0,0]
        if (i mod floor(nfiles/25.) eq 0) then begin
            plot, x, temp, title = t, /ylog, yra=[1,max(temp)], /nodata
            oplot, x, temp, color=fsc_color("red")
        endif
     endif

    ;Density and Temperature plot
    if (keyword_set(plotdens)) then begin
       plotsym, 0, .3, /fill
       ;; !x.omargin=[0,0]
       ;; !y.omargin=[0,0]

       if (i mod 40 eq 0 AND t lt 1.3d15) then begin
          ;; if (i eq 0) then begin
             plot, xpc, d, title = t, /ylog, yra=[1d-23,1d-20], /nodata, xsty=1, xtitle='X [pc]'
             axis, yaxis = 0, yrange=[1d-23, 1d-20], /ylog, ytitle = 'Density [g cm!E-3!N]', color=fsc_color("purple")
;             oplot, xpc, dneutral, color=fsc_color("blue");, thick = 1
             oplot, xpc, dneutral, color=fsc_color("blue");, psym = 8
             oplot, xpc, d, color=fsc_color("purple");, psym = 8

             oplot, [xpc[maxdloc], xpc[maxdloc]], [min(d),max(d)], linestyle = 2, thick=1
;             if (maxdloc gt 14) then begin
;                oplot, [xpc[xend[i]-1], xpc[xend[i]-1]], [min(d),max(d)], linestyle = 2
;                print, "Diff", xpc[xend[i]-1] - xpc[maxdloc]
;             endif 
             axis, yaxis = 1, yrange= [10, 10000], /save, ytitle = 'Temperature [K]', color=fsc_color("pink"), charsize = 1.3
             oplot, xpc, temp, color=fsc_color("pink");, psym = 8
         plot, xpc, ifrac, xsty = 1, ysty = 2, title = t, /nodata, xtitle = 'X [pc]', ytickformat="(A1)"
         axis, yaxis = 1, /save, ytitle = 'Ionization fraction'
            oplot, [-50,50], [0,0], linestyle =2, color=fsc_color("gray")
            oplot, [-50,50], [1,1], linestyle =2, color=fsc_color("gray")
            oplot, xpc, ifrac, color=fsc_color("navy"), psym = 8
         ;;  endif else begin
         ;;     plot, xpc, d, title = t, /ylog, yra=[1d-23,1d-19], /nodata, xsty=1, xtitle='X [pc]', ysty = 1
         ;;     oplot, xpc, dneutral, color=fsc_color("blue")
         ;;     oplot, xpc, d, color=fsc_color("purple")

         ;;     oplot, [xpc[maxdloc], xpc[maxdloc]], [min(d),max(d)], linestyle = 2
         ;;     ;axis, yaxis = 1, yrange= [10, 10000], ytickformat ='logticks_exp', /save
         ;;     oplot, xpc, temp*1d-24, color=fsc_color("pink")
         ;; plot, xpc, ifrac, xsty = 1, ysty = 2, title = t, /nodata, xtitle = 'X [pc]', ytitle= 'Ionization fraction'
         ;;    oplot, [-50,50], [0,0], linestyle =2, color=fsc_color("gray")
         ;;    oplot, [-50,50], [1,1], linestyle =2, color=fsc_color("gray")
         ;;    oplot, xpc, ifrac, color=fsc_color("navy")
         ;;  endelse
;          plot, x[i], dens[i], /ylog, yra=[1d-23,1d-20], /nodata, xsty=1
;          oplot, x[i], dens[i], color=fsc_color("magenta")
       endif
    endif


;     if (keyword_set(plotni2)) then begin
;         ni_analytic2 = sqrt(flux/alpha/pos)
;         plotsym, 0, .4, /fill
;          plot, xpc , d, ytitle = 'Density [g cm!E-3!N]', xtitle = 'Position [pc]', ysty = 2, xsty = 1, /ylog
;          oplot, pos, ni_analytic2, color=fsc_color("teal")
;          if (i mod 12 eq 0) then l = get_kbrd()
;      endif
endfor

;Close postscript files for timestep plots
if (keyword_set(plotvel)) then begin
   if(keyword_set(ps)) then begin
      device, /close
   endif
   stop
endif

if (keyword_set(plotifrac)) then begin
   if(keyword_set(ps)) then begin
      device, /close
   endif
   stop
endif

;     if (keyword_set(plotdens)) then begin
;      !x.omargin=[0,0]
;      !y.omargin=[0,0]
; !p.multi=0
;      plot, pos, dens, /ylog, yra=[1d-23,1d-20], /nodata, xsty=1
;      oplot, pos, dens, color=fsc_color("magenta")
;  endif

if (keyword_set(plotmass)) then begin
    if (keyword_set(ps)) then begin
        set_plot, 'ps'
        device, filen='coolfrac_mass.eps', xsize=10., ysize=10., /inches, /encapsulated
    endif

    !p.multi=[0,2,2]
    plot, times, massini, xsty = 1, ysty = 2, title='Swept up unionized mass/area ', xtitle='Time[s]', xra=[0, 8d13]
    plot, times, massnow, xsty = 1, ysty = 2, title='Mass in front/area', xtitle='Time[s]', xra=[0, 8d13]
    oplot, times, massnow2, linestyle =2

    plot, times, massini, xsty = 1, ysty = 2, /nodata, xtitle='Time[s]', xra=[0, 8d13]
    oplot, times, massini, color=fsc_color("blue"), linestyle =0
    oplot, times, massnow, color=fsc_color("cyan"), linestyle = 0
    ;; legend, ['Swept up unionized mass', 'Mass in front'], colors=[fsc_color("blue"), fsc_color("cyan")], psym =[8,8]
    plot, times, massnow/massini,  xsty =1, ysty = 2 , title='Front mass/Swept mass', xtitle='Time[s]', yra=[0, 2]
    if (keyword_set(ps)) then begin
        device, /close
    endif
    stop
endif



if (keyword_set(plotdens)) then begin
    if(keyword_set(ps)) then begin
        device, /close
    endif
stop
endif


if (keyword_set(plottemp)) then begin
    if(keyword_set(ps)) then begin
        device, /close
    endif
stop
endif

if (keyword_set(plotpres)) then begin
    if(keyword_set(ps)) then begin
        device, /close
    endif
stop
endif


 ;Analytic approximation
 cs2_t = kb * 1.d4 / mu_i
 print, min(temperature), max(temperature)
 coeff = ( 25./12. * sqrt(flux/alphab(temperature)) *cs2 /n0*(mu_i/mu_n))^.4
 pos_t = coeff * times^(4/5.) 


 coeff_ni = 25./12. *cs2 /n0*(1/1.7)
 print, sqrt(mean(cs2))
 ;pos_ni = sqrt(coeff_ni * niprod)*times
 ;pos_ni_t = flux/alpha/(niprod*niprod)

 if (keyword_set(plotni)) then begin
    if (keyword_set(ps)) then begin
       set_plot, 'ps'
       device, filen='numberdensityg_largeflux.eps', xsize=8., ysize=10.5, /inches, /encapsulated
    endif
    
    ni_analytic = sqrt(flux/alpha/pos_t)
    ni_analytic2 = sqrt(flux/alpha/pos)
    plotsym, 0, .4, /fill
    !p.multi=[0,1,2]
    plot, times, ni_analytic, ytitle = 'Mean number density [cm!E-3!N]', xtitle = 'Time [s]', ysty = 2, xsty = 1, /xlog,  xra=[7e12, 2d15], yra=[10, 100], /nodata, /ylog, title='Number density of Ionized Gas'
    oplot, times, ni_analytic, color=fsc_color("navy")
    oplot, times, ni_analytic2, color=fsc_color("teal")
    oplot, [1e12, 1e16], [n0, n0], linestyle = 2, color=fsc_color("gray")
    oplot, times, ni_sim
    legend, ['Simulated','Analytic using Simulated Sound Speed', 'Analytic using Simulated Position', 'Neutral number density'], $
            linestyle = [0,0,0,2], colors=[fsc_color("black"), fsc_color("navy"), fsc_color("teal"), fsc_color("gray")], /bottom
    plot, times, ni_analytic2/ni_sim,  /xlog,  xra=[7e12, 2d15], xtitle='Time[s]', ytitle = 'Analytic / Simulated number density', ysty=2, xsty=1
    if(keyword_set(ps)) then begin
       device, /close
    endif
stop
endif

 if (keyword_set(plotni2)) then begin
    if (keyword_set(ps)) then begin
       set_plot, 'ps'
       device, filen='numberdensity2_largeflux.eps', xsize=8., ysize=8., /inches, /encapsulated
    endif
    
    ni_analytic = sqrt(flux/alphab(temperature)/pos_t)
    ni_analytic2 = sqrt(flux/alphab(temperature)/pos)
    plotsym, 0, .4, /fill
    !p.multi=[0,1,1]
    tindex = where (times lt 2d15)
    plot, pos[1:*]/3.086d18 , ni_sim[1:*], ytitle = 'Mean number density [cm!E-3!N]', xtitle = 'Front Position [pc]', ysty = 2, xsty = 1, /xlog,  /nodata, /ylog, title='Number density of Ionized Gas', xra=[4,100]
    oplot, pos[tindex]/3.086d18, ni_analytic2[tindex], color=fsc_color("teal")
    oplot, pos[tindex]/3.086d18, ni_sim[tindex], psym =8
    oplot, [0, 10000], [63, 63], linestyle = 2, color=fsc_color("gray")
    ;; legend, ['Simulated', 'Analytic using Simulated Position'], $
    ;;              linestyle = [1,0], colors=[fsc_color("black"), fsc_color("teal")], /bottom
;        plot, times, ni_analytic2/ni_sim,  /xlog,  xra=[7e12, 2d15], xtitle='Time[s]', ytitle = 'Analytic / Simulated number density', ysty=2, xsty=1
    if(keyword_set(ps)) then begin
       device, /close
    endif
    stop
 endif
 

 ;Make plots
 if (keyword_set(ps)) then begin
    set_plot, 'ps'
    device, filen='reflecting_largeflux.eps', xsize=8., ysize=11., /inches, /encapsulated
 endif
 xerr = abs(x[2]-x[1])
 errors = dblarr(n_elements(pos)) + xerr/2

 !p.multi=[0,1,3]
 plotsym, 0, .4, /fill
 plot, times, pos_t, xstyle = 1, ystyle =1, xtitle = 'Time [s]', ytitle = 'X-Location of Max Density [cm]', title = flux, /nodata, /xlog, /ylog, xra=[1d13, 1.3d14]
 oploterror, times, pos, errors, psym = 8, /nohat
 oplot, times, pos_t, color=fsc_color("red"), psym =8
; legend, ['Simulated','Analytic'], psym = [8,8], colors=[fsc_color("black"), fsc_color("red")]

; plot, times, dens, xstyle = 1, ystyle =1, xtitle = 'Time [s]', ytitle = 'Max Density'
 plot, times, pos_t/pos, xsty = 1, ysty = 1, xtitle = 'Time [s]', ytitle='Analytic/Simulated' ,/xlog, psym=8,xra=[1d13, 1.3d14]
 
 plot, times, sqrt(cs2), xsty = 1, ysty = 2, xtitle = 'Time[s]', ytitle = 'Sound speed [cm/s]' ,/xlog, title = 'Ionized Sound Speed', psym=8, xra=[1d13, 1.3d14]


; plot, niprod, pos_ni, xtitle='Number density [cm^-3]', ytitle ='X-Location of Max Density [cm]'
; oplot, niprod, pos_ni_t, color=fsc_color("blue")
; oplot, niprod, pos, color=fsc_color("green")

if(keyword_set(ps)) then begin
   device, /close
endif

stop
END
