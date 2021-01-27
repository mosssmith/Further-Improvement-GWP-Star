pro itmo_test2
; test simple climate model used for metric calculations in AR5
; revise parameters to work entirely with tCO2, so NOT compatible with FaIR standard

data=read_csv('annual_rf_figure2.csv',header=header)
print,header
t_inp=data.(0)
nyr=n_elements(t_inp)
F_inp=fltarr(nyr,5)
F_inp(*,0)=data.(6)
F_inp(*,1)=data.(1)+data.(4)
F_inp(*,2)=data.(3)
F_inp(*,3)=data.(5)
F_inp(*,4)=data.(2)
F_tot=data.(7)

data=read_csv('fig2_rawems_plushist.csv',header=header)
print,header
t_obs=data.(0)
nob=n_elements(t_obs)
E_inp=fltarr(nob,3)
E_inp(*,0)=data.(1)*44./12.
E_inp(*,1)=data.(2)/1000.
E_inp(*,2)=data.(3)

;data=read_csv('California_cattle_emissions.csv',header=header)
;print,header
;t_cal=data.(0)
;E_cal=fltarr(nob)
;E_cal(where(t_obs ge min(t_cal) and t_obs le max(t_cal)))=data.(1)/25
;E_cal(where(t_obs le 1970))=11./25
;index=where(t_obs ge 1970 and t_obs le 2000,nin)
;print,index
;E_cal(index)=E_cal(index(0))+findgen(nin)*(E_cal(index(nin-1))-E_cal(index(0)))/nin
;index=where(t_obs ge max(t_cal),nin)
;drate=-0.003
;E_cal(index)=E_cal(index(0))*exp(drate*findgen(nin))

;p1=plot(t_obs,E_cal)

; first set up AR5 model parameters
print,nob
ny2=nob
tim2=findgen(ny2)+1

m_atm=5.1352d18 ; AR5 official mass of atmosphere in kg
m_air=28.97d-3  ; AR5 official molar mass of air
m_car=12.01d-3
m_co2=44.01d-3
m_ch4=16.043d-3

scl=1e3
a_ar5=fltarr(20)
a_ar5(0:3)=[0.21787,0.22896,0.28454,0.26863]
a_ar5(4)=1.e12*1.e6/m_co2/(m_atm/m_air) ; old value = 0.471 ; convert GtCO2 to ppm
a_ar5(5:8)=[1.e8,381.330,34.7850,4.12370]
a_ar5(9)=0. ; disable feedback
a_ar5(10:11)=[0.631,0.429]*0.7  ; 0.7 x AR5 sensitivity coeffs
a_ar5(13)=1.37e-2 ; rad efficiency in W/m2/ppm
a_ar5(14)=0. ; disable logarithmic forcing
a_ar5(15:16)=[8.400,409.5] ; AR5 thermal response parameters
a_ar5(17:19)=0.

a_ch4=a_ar5
a_ch4(0:3)=[0.,1.0,0.,0.]
;a_ch4(4)=1.e15*1.e9/16./1.81e20 ; convert GtCH4 to ppb
a_ch4(4)=1.e12*1.e9/m_ch4/(m_atm/m_air) ; convert GtCH4 to ppb
a_ch4(5:8)=[1.,12.4,1.,1.]
a_ch4(9)=0.
a_ch4(10:11)=[0.631,0.429]*0.7  ; 0.7 x AR5 sensitivity coeffs
a_ch4(13)=1.65*3.6324360e-4 ; RE W/m2/ppb
a_ch4(14)=0.
a_ch4(15:16)=[8.400,409.5]
a_ch4(17)=0.
a_ch4(18)=0.
a_ch4(19)=0.

AR5_AGWP=findgen(ny2,20)
AR5_AGWP(*,0)=a_ar5(4)*a_ar5(13)*a_ar5(0)*tim2
for j=1, 3 do AR5_AGWP(*,0)=AR5_AGWP(*,0)+a_ar5(4)*a_ar5(13)*a_ar5(j)*a_ar5(j+5)*(1.-exp(-tim2/a_ar5(j+5)))
AR5_AGWP(*,1)=a_ch4(4)*a_ch4(13)*a_ch4(1)*a_ch4(6)*(1.-exp(-tim2/a_ch4(6)))

H=100
s=[0.25,0.143,0.4]
g=(1.-exp(-s/(1-s)))/s
r=s/(1-s)/H
print,s,g,r
for i=0, 2 do AR5_AGWP(*,2+i)=AR5_AGWP(H,0)*(1-exp(-tim2*s(i)/(H*(1-s(i)))))/s(i)/g(i)

;stop

; compute operators to convert emissions to forcing
Fco2=EFmod(ny2,a_ar5)
Fch4=EFmod(ny2,a_ch4)

; and operator to convert forcing to temperature
Tar5=FTmod(ny2,a_ar5)

; operator to convert direct from emissions to temperature
Tco2=ETmod(ny2,a_ar5)
Tch4=ETmod(ny2,a_ch4)

; Calculate thermal pulse-response to methane and CO2
E4=fltarr(ny2)
E4(0)=1000. ; create 1000 GtCO2 pulse in year 0
; Compute forcing and temperature response
F4=Fco2#E4
T4r=Tar5#F4
T4=Tco2#E4

;Fco2cum=EFcum(ny2,a_ar5)
F4c=total(F4,/cumulative) ; Fco2cum#E4
;print,'Forcing matrix:'
;pm,Fco2(0:4,0:4)*1000.
;print,total(F4(0:10),/cumulative)
;print,F4c(0:10)

rhovals=0.005*(findgen(1000)+1.)/1000

r100=0.0047
a100=F4c(99)*r100/(1-exp(-r100*tim2(99)))
print,'rho,alpha(100,rho):',r100,a100
print,'Years, AGFP(CO2), AGWP(CO2), (alpha-AGFP)/AGWP, recalculated rho
for n=19, 100-1, 20 do begin
  ratvals=rhovals/(exp(rhovals*tim2(n))-1.)
  print,tim2(n),F4(n),F4c(n),(a100-F4(n))/F4c(n),interpol(rhovals,ratvals,F4(n)/F4c(n))
;  print,(1+r100*tim2(n))*exp(-r100*tim2(n))-1,(1+r100*tim2(n))*(1-r100*tim2(n))-1
endfor
;print,'Years, AGFP(CO2), AGWP(CO2), AGWP/H, d(AGWP/H)/dH, approx0, alpharho'
;for n=19, 100-1, 20 do begin
;  print,tim2(n),F4(n),F4c(n),F4c(n)/tim2(n),a100*(1-exp(-r100*tim2(n)))/r100/tim2(n),100*((F4(n)+F4(n+1))/2/tim2(n) - F4c(n)/tim2(n)^2),$
;    100*a100*((1+r100*tim2(n))*exp(-r100*tim2(n))-1)/r100/tim2(n)^2,100.*a100*(-r100/2.+r100^2*100./3.)
;  ;  print,(1+r100*tim2(n))*exp(-r100*tim2(n))-1,(1+r100*tim2(n))*(1-r100*tim2(n))-1
;endfor

stop

E5=fltarr(ny2)
E5(0)=1000. ; create 1000 GtCH4 pulse
; Compute forcing and temperature response
F5=Fch4#E5
T5r=Tar5#F5
T5=Tch4#E5

; Calculate CO2 warming-equivalent to E5, 1000 GtCH4 pulse of methane
E5s=fltarr(ny2)
GWP_H=AR5_AGWP(H,1)/AR5_AGWP(H,0)
print,'GWP_H:',GWP_H
c2=GWP_H*H*g(0)*(1-s(0))/20
c1=c2+GWP_H*g(0)*s(0)
print,'GWP* coefficients:',c1,c2,c1-c2
print,'GWP*/GWP coefficients:',c1/GWP_H,c2/GWP_H,c1/GWP_H-c2/GWP_H

E5s(0)=1000.*c1
; old GWP* definition
;E5s(20)=-1000.*c2
E5s(1:39)=-1000.*c2/39.

; Calculate methane warming-equivalent to E4, 1000 GtCO2 pulse of CO2
E4s=fltarr(ny2)
; old GWP* definition
;for n=0, ny2-1, 20 do E4s(n)=(1000./c1)*(c2/c1)^(n/20)
E4s(0)=1000./c1
for n=1, ny2-1 do for i=1, n<39 do E4s(n)=E4s(n)+E4s(n-i)*c2/c1/39.

; Compute forcing and temperature response
F4s=Fch4#E4s
T4s=Tch4#E4s

; Compute total methane equivalent to 1000 GtCO2 under GWP*:
print,'Total methane equivalent to 1000 GtCO2 under GWP*:'
print,(1000./c1)/(1-(c2/c1))

; Compute forcing and temperature response
F5s=Fco2#E5s
T5s=Tco2#E5s

; Compute old (2018) vintage GWP*
c0=GWP_H*H/20
E5s18=E5s*0.
E5s18(0)=1000.*c0
E5s18(20)=-1000.*c0
F5s18=Fco2#E5s18

p0=plot(tim2,AR5_AGWP(*,0),color='r',thick=3,xrange=[0,200],yrange=[0,0.2],$
  xtitle='Years of constant emissions',ytitle='Radiative forcing (W/m!U2!N per GtCO!D2!N/year)',font_size=12)
p0=plot(tim2,AR5_AGWP(*,2),thick=3,/overplot)
p0=plot(tim2,AR5_AGWP(*,3),line=1,/overplot)
p0=plot(tim2,AR5_AGWP(*,4),line=2,/overplot)
p0=plot(tim2,AR5_AGWP(*,1)/GWP_H,thick=3,color='b',/overplot)
p0=plot(tim2,total(F5s,/cumulative)/1000./GWP_H,thick=3,color='b',line=2,/overplot)
p0=plot(tim2,total(F5s,/cumulative)/1000./GWP_H*112/128,thick=1,color='b',line=2,/overplot)
p0=plot(tim2,total(F5s18,/cumulative)/1000./GWP_H,thick=1,color='b',line=1,/overplot)
;p0=plot(tim2,AR5_AGWP(*,3),line=1,/overplot)
;p0=plot(tim2,AR5_AGWP(*,4),line=2,/overplot)

;stop

p0.save,'ITMO_CO2_AGWP_2.png'

p0b=plot(tim2,AR5_AGWP(*,0)*1000./tim2,color='r',thick=4,xrange=[0,200],$
  xtitle='Years of constant emissions',ytitle='AGWP!DH!N/H (W/m!U2!N per 1000 GtCO!D2!N)',font_size=12)
p0b=plot(tim2,AR5_AGWP(*,2)*1000./tim2,thick=2,/overplot)
p0b.save,'ITMO_CO2_AGWP_3.png'

stop
; Compute methane emissions to give forcing F4 (response to 1000 GtCO2 pulse)
E4f=invert(Fch4)#F4
print,'Methane emissions equivalent to 1000 GtCO2 pulse (MtCH4)'
print,E4f(0:5),'50-yr avg:',mean(E4f(2:50))

;E4f=invert(Tch4)#T4
; Compute CO2 emissions to give forcing F5 (response to 1000 GtCH4 pulse)
E5f=invert(Fco2)#F5
print,'CO2 emissions equivalent to 1 GtCH4 pulse (GtCO2)'
print,E5f(0:5)/1000.,'50-yr avg:',mean(E5f(2:50))/1000.


;E5f=invert(Tco2)#T5
Finvco2=invert(Fco2)
;print,'Inverse CO2 GFP matrix:'
;pm,Finvco2(0:5,0:5)
Finvch4=invert(Fch4)
;print,'Inverse CH4 GFP matrix:'
;pm,Finvch4(0:5,0:5)

; calculate CO2-warming-equivalent methane emissions
E_ch4_we=Finvco2#Fch4#E_inp(*,1)
; calculate GWP* warming-equivalent methane emissions
E_ch4_st=E_inp(*,1)*c1
; old-style calculation
;E_ch4_st(20:*)=E_inp(20:*,1)*128.-E_inp(0:nob-21,1)*120.
for n=1, nob-1 do E_ch4_st(n)=E_inp(n,1)*c1-total(E_inp(0>(n-39):n-1,1))*c2/39.

;p1g=plot(t_obs,E_inp(*,0),yrange=[-10,50],color='r',thick=3,xtitle='Years',ytitle='Emissions (billion tCO!D2!N-"equivalent"/year)',font_size=12,axis_style=1)
;p1g=plot(t_obs,E_inp(*,1)*28.,color='b',thick=3,/overplot)
;p1g=plot(t_obs,t_obs*0.,thick=1,/overplot)
;p1g=plot([2016,2016],[-20,60],thick=1,line=1,/overplot)
;t1g=text(0.22,0.82,'Carbon dioxide emissions', color='red',font_size=12)
;t1g=text(0.22,0.77,'Methane using GWP!D100!N',color='blue',font_size=12)

;p1g.save,'CO2_CH4_GWP100.png'

;p1g=plot(t_obs,E_inp(*,0),yrange=[-10,50],color='r',thick=3,xtitle='Years',ytitle='Emissions (billion tCO!D2!N-"equivalent"/year)',font_size=12,axis_style=1)
;p1g=plot(t_obs,E_inp(*,1)*28.,color='b',line=1,/overplot)
;p1g=plot(t_obs,E_inp(*,1)*84.,color='b',thick=3,/overplot)
;p1g=plot(t_obs,t_obs*0.,thick=1,/overplot)
;p1g=plot([2016,2016],[-10,50],thick=1,line=1,/overplot)
;t1g=text(0.22,0.82,'Carbon dioxide emissions', color='red',font_size=12)
;t1g=text(0.22,0.77,'Methane using GWP!D20!N',color='blue',font_size=12)

;p1g.save,'CO2_CH4_GWP20.png'

p1g=plot(t_obs,E_inp(*,0),color='r',thick=3,yrange=[-20,50],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years',ytitle='Emissions (billion tCO!D2!N/year)',font_size=12)
p1g=plot(t_obs,E_inp(*,1)*28.4,color='b',thick=1,line=0,/overplot)
p1g=plot(t_obs,E_inp(*,1)*83.6,color='b',thick=1,line=3,/overplot)
p1g=plot(t_obs,E_inp(*,1)*4.25,color='b',thick=1,line=1,/overplot)
p1g=plot(t_obs,E_ch4_we,color='purple',thick=3,/overplot)
p1g=plot(t_obs,E_ch4_st,color='purple',line=2,thick=1,/overplot)
p1g=plot(t_obs,t_obs*0.,thick=1,line=0,/overplot)
p1g=plot([2015,2015],[0,60],thick=1,line=1,/overplot)
t1g=text(0.12,0.82,'a) Annual GHG emissions', font_size=12)
t1g=text(0.12,0.25,'Carbon dioxide',font_color='red', font_size=12)
t1g=text(0.12,0.2,'Methane metric-equivalent',color='blue', font_size=12)
t1g=text(0.12,0.15,'Methane warming-equivalent',color='purple', font_size=12)
t1g=text(2060,15,'GWP!D20!N',font_color='blue', font_size=10, /data)
t1g=text(2060,6,'GWP!D100!N',font_color='blue', font_size=10, /data)
t1g=text(2060,1.5,'GTP!D100!N',font_color='blue', font_size=10, /data)
t1g=text(2070,-5,'LWE',font_color='purple', font_size=10, /data)
t1g=text(2045,-12,'GWP*',font_color='purple', font_size=10, /data)

T_inp=E_inp
T_inp(*,0)=Tco2#E_inp(*,0)
T_inp(*,1)=Tch4#E_inp(*,1)

;p2g=plot(t_obs,T_inp(*,0),yrange=[0,2.],color='r',thick=3,xtitle='Years',ytitle='Warming (!Uo!NC since 1850)',font_size=12,axis_style=1)
;p2g=plot(t_obs,T_inp(*,1),color='b',thick=3,/overplot)
;p2g=plot(t_obs,t_obs*0.,thick=1,line=1,/overplot)
;p2g=plot([2016,2016],[-10,50],thick=1,line=1,/overplot)
;t2g=text(0.22,0.82,'Carbon-dioxide-induced warming',color='red',font_size=12)
;t2g=text(0.22,0.77,'Methane-induced warming',color='blue',font_size=12)

TCRE=max(T_inp(*,0))/max(total(E_inp(*,0),/cumulative))
print,'TCRE=',TCRE

p2g=plot(t_obs,T_inp(*,0),yrange=[0,2.],color='r',thick=6,xtitle='Years',ytitle='Warming (!Uo!NC since 1850)',dimensions=[1200,500],$
  position=[0.6,0.1,0.95,0.9],/current,font_size=12,/nodata)
p2g=plot(t_obs,T_inp(*,0),color='r',thick=6,/overplot,transparency=75)
p2g=plot(t_obs,T_inp(*,1),color='b',thick=6,/overplot,transparency=75)
p2g=plot(t_obs,T_inp(*,0)+T_inp(*,1),thick=6,/overplot,transparency=75)
p2g=plot(t_obs,t_obs*0.,thick=1,line=0,/overplot)
p2g=plot([2015,2015],[0,1.8],thick=1,line=1,/overplot)
p2g=plot(t_obs,total(E_inp(*,0),/cumulative)*TCRE,color='r',thick=3,/overplot)
p2g=plot(t_obs,total(E_ch4_we,/cumulative)*TCRE,color='purple',thick=3,/overplot)
p2g=plot(t_obs,total(E_ch4_st,/cumulative)*TCRE,color='purple',line=2,/overplot)
p2g=plot(t_obs,total(E_inp(*,1)*28.4,/cumulative)*TCRE,color='b',/overplot,line=0)
p2g=plot(t_obs,total(E_inp(*,0)+E_ch4_we,/cumulative)*TCRE,thick=3,/overplot)
p2g=plot(t_obs,total(E_inp(*,0)+E_ch4_st,/cumulative)*TCRE,/overplot,line=2)
p2g=plot(t_obs,total(E_inp(*,0)+E_inp(*,1)*4.25,/cumulative)*TCRE,/overplot,line=1)
p2g=plot(t_obs,total(E_inp(*,0)+E_inp(*,1)*28.4,/cumulative)*TCRE,/overplot,line=0)
p2g=plot(t_obs,total(E_inp(*,0)+E_inp(*,1)*83.6,/cumulative)*TCRE,/overplot,line=3)
ax=p2g.axes
ax[3].hide=1
yaxis = AXIS('Y', LOCATION='right', target=p2g, $
  TITLE='Cumulative emissions (trillion tCO!D2!N)',tickfont_size=12, $
  COORD_TRANSFORM=[0,1./1000/TCRE], $
  AXIS_RANGE=[0,2]/1000/TCRE)
t1g=text(928./1200,360./500,'GWP!D20!N', font_size=10)
t1g=text(1060./1200,380./500,'GWP!D100!N', font_size=10)
t1g=text(1050./1200,298./500,'GTP!D100!N', font_size=10)
t1g=text(1068./1200,352./500,'LWE', font_size=10)
t1g=text(1098./1200,335./500,'GWP*', font_size=10)
t1g=text(1060./1200,160./500,'GWP!D100!N', font_color='blue', font_size=10)
t1g=text(1044./1200,130./500,'LWE', font_color='purple', font_size=10)
t1g=text(1080./1200,122./500,'GWP*', font_color='purple', font_size=10)
t1g=text(0.62,0.75,'Carbon dioxide',font_color='red', font_size=12)
t1g=text(0.62,0.70,'Methane',color='blue', font_size=12)
t1g=text(0.62,0.65,'Combined', font_size=12)

t2g=text(0.62,0.82,'b) Warming and cumulative emissions',font_size=12)

p1g.save,'CO2_CH4_CO2we_temp.png'

;stop
;p2g.save,'CO2_CH4_temp.png'

; calculate CO2-warming-equivalent methane emissions
;E_cal_we=gauss_smooth(Finvco2#Fch4#E_cal,1)

;p1g=plot(t_obs,E_cal*28.,xrange=[1970,2070],yrange=[0,50],color='b',thick=3,xtitle='Years',ytitle='Emissions (million tCO!D2!N-"equivalent"/year)',font_size=12,axis_style=1)
;p1g=plot(t_obs,E_cal_we,color='purple',thick=3,/overplot)
;t1g=text(0.22,0.82,'California agricultural methane emissions using GWP!D100!N', color='blue',font_size=12)
;t1g=text(0.22,0.77,'CO!D2!N-warming-equivalent methane emissions',color='purple',font_size=12)

;p1g.save,'CO2_CH4_Calwe.png'

;T_cal=Tch4#E_cal
;T_cal_co2=Tco2#E_cal*28.

;p2g=plot(t_obs,T_cal,xrange=[1970,2070],color='b',thick=3,xtitle='Years',ytitle='Warming (mK since 1850)',font_size=12,axis_style=1)
;p2g=plot(t_obs,T_cal_co2,color='red',thick=3,/overplot)
;t2g=text(0.22,0.82,'Warming caused by Californian methane',color='blue',font_size=12)
;t2g=text(0.22,0.77,'Warming caused by "equivalent" CO!D2!N using GWP!D100!N',color='red',font_size=12)

;p2g.save,'CO2_CH4_Caltemp.png'

;stop

E4sc=total(E4s,/cumulative)
E5sc=total(E5s,/cumulative)
E4fc=total(E4f,/cumulative)
E5fc=total(E5f,/cumulative)

p=barplot(tim2,E4s,color='b',fill_color='b',xrange=[0,200],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years after implementation',ytitle='Emissions (MtCH!D4!N per GtCO!D2!N)',transparency=75,font_size=12)
p=barplot(tim2,E4f,color='r',fill_color='r',/overplot,transparency=75)
t=text(0.12,0.82,'a) Methane equivalent to pulse CO!D2!N',font_size=12)
p=barplot(tim2,E4sc,color='b',fill_color='b',xrange=[0,200],position=[0.6,0.1,0.95,0.9],/current,$
  xtitle='Years after start of implementation',ytitle='Emissions (MtCH!D4!N per GtCO!D2!N/yr)',transparency=75,font_size=12)
p=barplot(tim2,E4fc,color='r',fill_color='r',/overplot,transparency=75)
t=text(0.62,0.82,'b) Methane equivalent to sustained CO!D2!N',font_size=12)

p.save,'ITMO_CH4forCO2emm.png'

p=barplot(tim2,E5f/1000.,thick=2,color='r',fill_color='r',transparency=75,xrange=[0,200],yrange=[-10,130],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years after implementation',ytitle='Emissions (MtCO!D2!N per MtCH!D4!N)',font_size=12)
p=barplot(tim2,E5s/1000.,thick=2,color='b',fill_color='b',transparency=75,/overplot)
t=text(0.12,0.82,'a) CO!D2!N equivalent to pulse methane',font_size=12)
p=barplot(tim2,E5sc/1000.,thick=2,color='b',fill_color='b',transparency=75,xrange=[0,200],yrange=[0,150],position=[0.6,0.1,0.95,0.9],/current,$
  xtitle='Years after start of implementation',ytitle='Emissions (MtCO!D2!N per MtCH!D4!N/yr)',font_size=12)
p=barplot(tim2,E5fc/1000.,thick=2,color='r',fill_color='r',transparency=75,/overplot)
t=text(0.62,0.82,'b) CO!D2!N equivalent to sustained methane',font_size=12)

p.save,'ITMO_CO2forCH4emm.png'
;stop

print,'A_co2 in W/m2/kg:',a_ar5(4)*a_ar5(13)*(m_car/m_co2)/1.e12
print,'A_ch4 in W/m2/kg:',a_ch4(4)*a_ch4(13)/1.e12

F4C=total(F4,/cumulative)
F5C=total(F5,/cumulative)
T4C=total(T4,/cumulative)/1000.
T5C=total(T5,/cumulative)/1000.
T4Cs=total(T4s,/cumulative)/1000.
T5Cs=total(T5s,/cumulative)/1000.

index=[0,9,19,49,99]
print,'Years to check: ',tim2(index)
print,'AR5 formula values:'
print,'CO2 AGWP (x10e15):',AR5_AGWP(index,0)*1000.
print,'CH4 AGWP (x10e15):',AR5_AGWP(index,1)*1000.

print,'CO2 AGWP (x10e15):',F4C(index)
print,'CH4 AGWP (x10e15):',F5C(index)
print,'CH4 GWP          :',F5C(index)/F4C(index)
print,'CO2 AGTP (x10e15):',T4(index)
print,'CO2 AGTP (recalc):',T4r(index)
print,'CH4 AGTP (x10e15):',T5(index)
print,'CH4 AGTP (recalc):',T5r(index)
print,'CH4 GTP          :',T5(index)/T4(index)
print,'CO2 iAGTP(x10e12):',T4C(index)
print,'CH4 iAGTP(x10e12):',T5C(index)
print,'CH4 iGTP         :',T5C(index)/T4C(index)
print,'CH4 iGTP/GWP     :',T5C(index)/T4C(index)/(F5C(index)/F4C(index))

p6=plot(tim2,replicate(0.,ny2),linestyle=0,thick=1,xrange=[0,200],yrange=[-2,2],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years after implementation',ytitle='Impact on GMST (!Uo!NC per 1000 GtCO!D2!N offset)',font_size=12)
p6=plot(tim2,T4,thick=4,color='light grey',/overplot)
p6=plot(tim2,T4-T5/28,thick=3,/overplot)
p6=plot(tim2,T4-T5/84,thick=3,linestyle=3,/overplot)
p6=plot(tim2,T4-T5/4.2,thick=3,linestyle=1,/overplot)
p6=plot(tim2,T4-T4s,thick=3,linestyle=2,/overplot)
t6=text(0.12,0.82,'a) Methane offsetting pulse emission of CO!D2!N',font_size=12)
t6=text(70,-0.6,/data,'Metric used to calculate offset:',font_size=12)
p6=plot([100,125],[-0.8,-0.8],thick=4,color='light grey',/overplot)
t6=text(130,-0.85,/data,'No offset',font_size=12)
p6=plot([100,125],[-1.0,-1.0],thick=3,/overplot)
t6=text(130,-1.05,/data,'GWP!D100!N',font_size=12)
p6=plot([100,125],[-1.2,-1.2],thick=3,linestyle=3,/overplot)
t6=text(130,-1.25,/data,'GWP!D20!N',font_size=12)
p6=plot([100,125],[-1.4,-1.4],thick=3,linestyle=1,/overplot)
t6=text(130,-1.45,/data,'GTP!D100!N',font_size=12)
p6=plot([100,125],[-1.6,-1.6],thick=3,linestyle=2,/overplot)
t6=text(130,-1.65,/data,'GWP*',font_size=12)

p7=plot(tim2,replicate(0.,ny2),linestyle=0,thick=1,xrange=[0,200],yrange=[-0.1,0.1],position=[0.6,0.1,0.95,0.9],/current,$
  xtitle='Years after start of implementation',ytitle='Impact on GMST (!Uo!NC per GtCO!D2!N/yr offset)',font_size=12)
p7=plot(tim2,T4C,thick=4,color='light grey',/overplot)
p7=plot(tim2,T4C-T5C/28,thick=3,/overplot)
p7=plot(tim2,T4C-T5C/84,thick=3,linestyle=3,/overplot)
p7=plot(tim2,T4C-T5C/4.2,thick=3,linestyle=1,/overplot)
p7=plot(tim2,T4C-T4Cs,thick=3,linestyle=2,/overplot)
t7=text(0.62,0.82,'b) Methane offsetting sustained emission of CO!D2!N',font_size=12)

p6.save,'ITMO_CH4forCO2new.png'

p8=plot(tim2,replicate(0.,ny2),linestyle=0,thick=1,xrange=[0,200],yrange=[-0.1,0.1],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years after implementation',ytitle='Impact on GMST (!Uo!NC per GtCH!D4!N offset)',font_size=12)
p8=plot(tim2,T5/1000.,thick=4,color='light grey',/overplot)
p8=plot(tim2,T5/1000.-T4*28./1000.,thick=3,/overplot)
p8=plot(tim2,T5/1000.-T4*4.2/1000.,thick=3,linestyle=1,/overplot)
p8=plot(tim2,T5/1000.-T4*84/1000.,thick=3,linestyle=3,/overplot)
p8=plot(tim2,T5/1000.-T5s/1000.,thick=3,linestyle=2,/overplot)
t8=text(0.12,0.82,'a) CO!D2!N offsetting pulse emission of methane',font_size=12)
t8=text(70,0.06,/data,'Metric used to calculate offset:',font_size=12)
p8=plot([100,125],[0.05,0.05],thick=4,color='light grey',/overplot)
t8=text(130,0.0475,/data,'No offset',font_size=12)
p8=plot([100,125],[0.04,0.04],thick=3,/overplot)
t8=text(130,0.0375,/data,'GWP!D100!N',font_size=12)
p8=plot([100,125],[0.03,0.03],thick=3,linestyle=3,/overplot)
t8=text(130,0.0275,/data,'GWP!D20!N',font_size=12)
p8=plot([100,125],[0.02,0.02],thick=3,linestyle=1,/overplot)
t8=text(130,0.0175,/data,'GTP!D100!N',font_size=12)
p8=plot([100,125],[0.01,0.01],thick=3,linestyle=2,/overplot)
t8=text(130,0.0075,/data,'GWP*',font_size=12)
;t1=text(70,-0.08,/data,'Metric used to calculate offset:',font_size=12)

p9=plot(tim2,replicate(0.,ny2),linestyle=0,thick=1,xrange=[0,200],yrange=[-2.5,2.5],position=[0.6,0.1,0.95,0.9],/current,$
  xtitle='Years after start of implementation',ytitle='Impact on GMST (!Uo!NC per GtCH!D4!N/yr offset)',font_size=12)
p9=plot(tim2,T5C,thick=4,color='light grey',/overplot)
p9=plot(tim2,T5C-T4C*28,thick=3,/overplot)
p9=plot(tim2,T5C-T4C*84,thick=3,linestyle=3,/overplot)
p9=plot(tim2,T5C-T4C*4.2,thick=3,linestyle=1,/overplot)
p9=plot(tim2,T5C-T5Cs,thick=3,linestyle=2,/overplot)
t9=text(0.62,0.82,'b) CO!D2!N offsetting sustained emission of methane',font_size=12)

p8.save,'ITMO_CO2forCH4new.png'

; generate square-pulse emissions in MtCO2 for a 40-year lifetime 500MW gas plant with CO2 emissions of 0.5tCO2/MWh
E6=fltarr(ny2)
E6(10:49)=0.5*24*365*500/1e6
lrate=0.015
E7=E6*lrate*m_ch4/m_co2
E8=invert(Fco2)#Fch4#E7
p=plot(tim2-0.5,E6,color='r',fill_color='r',xrange=[0,100],dimensions=[1200,500],position=[0.1,0.1,0.45,0.9],$
  xtitle='Years',ytitle='Emissions (MtCO!D2!N-we/year)',thick=2,font_size=12)
p=plot(tim2-0.5,E8,color='b',fill_color='b',/overplot,thick=2)
p=plot(tim2-0.5,E7*28,color='b',fill_color='b',/overplot,line=2)
p=plot(tim2-0.5,E7*84,color='b',fill_color='b',/overplot,line=3)
t8=text(0.12,0.82,'a) Contributions to warming rate',font_size=12)
p=plot(tim2-0.5,total(E6,/cumulative),color='r',fill_color='r',xrange=[0,100],dimensions=[1200,500],position=[0.6,0.1,0.95,0.9],/current,$
  xtitle='Years',ytitle='Cumulative emissions (MtCO!D2!N-we)',thick=2,font_size=12)
p=plot(tim2-0.5,total(E8,/cumulative),color='b',fill_color='b',/overplot,thick=2)
p=plot(tim2-0.5,total(E7*28,/cumulative),color='b',fill_color='b',/overplot,line=2)
p=plot(tim2-0.5,total(E7*84,/cumulative),color='b',fill_color='b',/overplot,line=3)
t9=text(0.62,0.82,'b) Contributions to warming',font_size=12)

p.save,'CHLS_plot.png'

return
end

function EFmod,nyr,a
  Fcal=fltarr(nyr,nyr)
  time=findgen(nyr+1)
  F_0=a(4)*a(13)*a(0)*time
  for j=1, 3 do F_0=F_0+a(j)*a(4)*a(13)*a(j+5)*(1-exp(-time/a(j+5)))
  for i=0, nyr-1 do Fcal(i,0)=F_0(i+1)-F_0(i)
  for j=1, nyr-1 do Fcal(j:nyr-1,j)=Fcal(0:nyr-j-1,0)
  return,Fcal
end

function FTmod,nyr,a
  Tcal=fltarr(nyr,nyr)
  time=findgen(nyr)+0.5
  for j=0, 1 do Tcal(*,0)=Tcal(*,0)+(a(j+10)/a(j+15))*exp(-time/a(j+15))
  for j=0, nyr-1 do Tcal(j:nyr-1,j)=Tcal(0:nyr-j-1,0)
  return,Tcal
end

function ETmod,nyr,a
  Tcal=fltarr(nyr,nyr)
  time=findgen(nyr)+1
  for j=0, 1 do begin
    Tcal(*,0)=Tcal(*,0)+a(4)*a(13)*a(0)*a(j+10)*(1-exp(-time/a(j+15)))
    for i=1, 3 do begin
      Tcal(*,0)=Tcal(*,0)+a(4)*a(13)*a(i)*a(i+5)*a(j+10)*(exp(-time/a(i+5))-exp(-time/a(j+15)))/(a(i+5)-a(j+15))      
    endfor
  endfor
  for j=0, nyr-1 do Tcal(j:nyr-1,j)=Tcal(0:nyr-j-1,0)
  return,Tcal
end

