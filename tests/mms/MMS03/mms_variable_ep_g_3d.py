# Purpose: SymPy script to generate manufactured solution source terms
# by substituting the manufactured solution functions into the 
# governing equations.
#
# To run: $ python mms_variable_ep_g_3d.py > mms_variable_ep_g_3d.dat
# 
# Notes: 
#   - for incompressible and variable volume fractions, 
#     only gas continuity is being solved.
#
# Author: Aniruddha Choudhary (aniruddhac@gmail.com)
#
# Date: March 10, 2015
#######################################################################

# Import sympy. Sympy library should be in PYTHONPATH environment
# variable.
from sympy import *

# common symbols:
pi = symbols('pi')

# symbols for mms functions:
xt, yt, zt = symbols('xt yt zt')
ug, vg, wg = symbols('ug vg wg')
us, vs, ws = symbols('us vs ws')
pg = symbols('pg')
tg, ts, ths = symbols('tg ts ths')
eg, es, rops = symbols('eg es rops')

# symbols for mms source terms:
ugSrc, vgSrc, wgSrc = symbols('ugSrc vgSrc wgSrc')
usSrc, vsSrc, wsSrc = symbols('usSrc vsSrc wsSrc')
TgSrc, TsSrc, ThsSrc = symbols('TgSrc TsSrc ThsSrc')
ropgSrc, ropsSrc = symbols('ropgSrc ropsSrc')

# symbols for physical variables:
ros = symbols('ros')
rog = symbols('rog')
mug = symbols('mug')
mus = symbols('mus')
Cpg, kg = symbols('Cpg kg')
Cps, ks = symbols('Cps ks')

# symbols for MMS constants:
ug0,ugx,ugy,ugz,ugxy,ugyz,ugzx = symbols('ug0 ugx ugy ugz ugxy ugyz ugzx')
augx,augy,augz,augxy,augyz,augzx = symbols('augx augy augz augxy augyz augzx')

vg0,vgx,vgy,vgz,vgxy,vgyz,vgzx = symbols('vg0 vgx vgy vgz vgxy vgyz vgzx')
avgx,avgy,avgz,avgxy,avgyz,avgzx = symbols('avgx avgy avgz avgxy avgyz avgzx')

wg0,wgx,wgy,wgz,wgxy,wgyz,wgzx = symbols('wg0 wgx wgy wgz wgxy wgyz wgzx')
awgx,awgy,awgz,awgxy,awgyz,awgzx = symbols('awgx awgy awgz awgxy awgyz awgzx')

us0, vs0, ws0 = symbols('us0 vs0 ws0')

pg0,pgx,pgy,pgz,pgxy,pgyz,pgzx = symbols('pg0 pgx pgy pgz pgxy pgyz pgzx')
apgx,apgy,apgz,apgxy,apgyz,apgzx = symbols('apgx apgy apgz apgxy apgyz apgzx')

tg0,tgx,tgy,tgz,tgxy,tgyz,tgzx = symbols('tg0 tgx tgy tgz tgxy tgyz tgzx')
atgx,atgy,atgz,atgxy,atgyz,atgzx = symbols('atgx atgy atgz atgxy atgyz atgzx')

ts0,tsx,tsy,tsz,tsxy,tsyz,tszx = symbols('ts0 tsx tsy tsz tsxy tsyz tszx')
atsx,atsy,atsz,atsxy,atsyz,atszx = symbols('atsx atsy atsz atsxy atsyz atszx')

ths0,thsx,thsy,thsz,thsxy,thsyz,thszx = symbols('ths0 thsx thsy thsz thsxy thsyz thszx')
athsx,athsy,athsz,athsxy,athsyz,athszx = symbols('athsx athsy athsz athsxy athsyz athszx')

es0,esx,esy,esz,esxy,esyz,eszx = symbols('es0 esx esy esz esxy esyz eszx')
aesx,aesy,aesz,aesxy,aesyz,aeszx = symbols('aesx aesy aesz aesxy aesyz aeszx')

# symbols for divergence
divug, divus = symbols('divug divus')

# gas pressure variable
pg = ( pg0 + pgx*cos(apgx*pi*xt) + pgy*cos(apgy*pi*yt) + pgxy*cos(apgxy*pi*xt*yt) + 
       pgzx*cos(apgzx*pi*xt*zt) + pgz*sin(apgz*pi*zt) + pgyz*sin(apgyz*pi*yt*zt) )

# gas temperature variable
tg = ( tg0 + tgx*cos(atgx*pi*xt) + tgy*cos(atgy*pi*yt) + tgxy*cos(atgxy*pi*xt*yt) + 
       tgzx*cos(atgzx*pi*xt*zt) + tgz*sin(atgz*pi*zt) + tgyz*sin(atgyz*pi*yt*zt) )

# solid temperature variable
ts = ( ts0 + tsx*cos(atsx*pi*xt) + tsy*cos(atsy*pi*yt) + tsxy*cos(atsxy*pi*xt*yt) + 
       tszx*cos(atszx*pi*xt*zt) + tsz*sin(atsz*pi*zt) + tsyz*sin(atsyz*pi*yt*zt) )

# solid volume fraction
es = ( es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esz*sin(aesz*pi*zt) )
#es = ( es0 + esx*cos(aesx*pi*xt) + esy*cos(aesy*pi*yt) + esxy*cos(aesxy*pi*xt*yt) + 
#       eszx*cos(aeszx*pi*xt*zt) + esz*sin(aesz*pi*zt) + esyz*sin(aesyz*pi*yt*zt) )

# gas volume fraction
eg = 1.0 - es

# rops
rops = ros*es

# granular temperature
#ths = ( ths0 + thsx*cos(athsx*pi*xt) + thsy*cos(athsy*pi*yt) + 
#      thsxy*cos(athsxy*pi*xt*yt) + thszx*cos(athszx*pi*xt*zt) + 
#      thsz*sin(athsz*pi*zt) + thsyz*sin(athsyz*pi*yt*zt) )

# vector field of general sinusoidal functions
Fvecx, Fvecy, Fvecz = symbols('Fvecx Fvecy Fvecz')

Fvecx = ( ug0 + ugx*sin(augx*pi*xt) + ugy*cos(augy*pi*yt) + 
  ugz*cos(augz*pi*zt) + ugxy*cos(augxy*pi*xt*yt) + 
  ugyz*sin(augyz*pi*yt*zt) + ugzx*cos(augzx*pi*xt*zt) )

Fvecy = ( vg0 + vgx*sin(avgx*pi*xt) + vgy*cos(avgy*pi*yt) + 
  vgz*cos(avgz*pi*zt) + vgxy*cos(avgxy*pi*xt*yt) + 
  vgyz*sin(avgyz*pi*yt*zt) + vgzx*cos(avgzx*pi*xt*zt) )

Fvecz = ( wg0 + wgx*cos(awgx*pi*xt) + wgy*sin(awgy*pi*yt) + 
  wgz*cos(awgz*pi*zt) + wgxy*sin(awgxy*pi*xt*yt) + 
  wgyz*sin(awgyz*pi*yt*zt) + wgzx*cos(awgzx*pi*xt*zt) )

# gas velocity variables (i.e., curl of F)
ug = ( diff(Fvecz, yt) - diff(Fvecy, zt) )/eg

vg = ( diff(Fvecx, zt) - diff(Fvecz, xt) )/eg

wg = ( diff(Fvecy, xt) - diff(Fvecx, yt) )/eg

# check gas continuity
divug = diff(ug,xt)+diff(vg,yt)+diff(wg,zt)
print('\n Gas velocity field meets continuity criteria:')
print(simplify(diff(eg*ug,xt)+diff(eg*vg,yt)+diff(eg*wg,zt))==0)

# solid velocity variable
us = ( us0*sin((pi*(xt + yt + zt))/2.0)**2 )/es

vs = ( vs0*cos((pi*(xt + yt + zt))/2.0)**2 )/es

ws = ( ws0 )/es

# check solid continuity
divus = (diff(us,xt)+diff(vs,yt)+diff(ws,zt)).subs(us0,vs0)
print('\n Solid velocity field meets continuity criteria:')
print(simplify((diff(es*us,xt)+diff(es*vs,yt)+diff(es*ws,zt)).subs(us0,vs0))==0)

# print mms functions
print('\n ****** MMS Functions ******')

print('\n pg =')
print(pg)

print('\n ug =')
print(ug)

print('\n vg =')
print(vg)

print('\n wg =')
print(wg)

print('\n us =')
print(us)

print('\n vs =')
print(vs)

print('\n ws =')
print(ws)

print('\n tg =')
print(tg)

print('\n ts =')
print(ts)

print('\n ths =')
print(ths)

print('\n es =')
print(es)

print('\n rops =')
print(rops)

# mms source terms
#ugSrc = ( rog*eg*(diff(ug**2, xt) + diff(ug*vg, yt) + diff(ug*wg, zt)) - 
#   mug*(diff(ug, xt, 2) + diff(ug, yt, 2) + diff(ug, zt, 2)) - (diff(
#      mug*diff(ug, xt), xt) + diff(mug*diff(vg, xt), yt) + 
#     diff(mug*diff(wg, xt), zt) + diff(-Rational(2,3)*mug*divug, xt)) + eg*diff(pg, xt) )
#
#vgSrc = ( rog*eg*(diff(ug*vg, xt) + diff(vg**2, yt) + diff(vg*wg, zt)) - 
#   mug*(diff(vg, xt, 2) + diff(vg, yt, 2) + diff(vg, zt, 2)) - (diff(
#      mug*diff(ug, yt), xt) + diff(mug*diff(vg, yt), yt) + 
#     diff(mug*diff(wg, yt), zt) + diff(-Rational(2,3)*mug*divug, yt)) + eg*diff(pg, yt) )
#
#wgSrc = ( rog*eg*(diff(ug*wg, xt) + diff(vg*wg, yt) + diff(wg**2, zt)) - 
#  mug*(diff(wg, xt, 2) + diff(wg, yt, 2) + diff(wg, zt, 2)) - (diff(
#     mug*diff(ug, zt), xt) + diff(mug*diff(vg, zt), yt) + 
#    diff(mug*diff(wg, zt), zt) + diff(-Rational(2,3)*mug*divug, zt)) + eg*diff(pg, zt) )
#
#usSrc = ( rops*(diff(us**2, xt) + diff(us*vs, yt) + diff(us*ws, zt)) - 
#  mus*(diff(us, xt, 2) + diff(us, yt, 2) + diff(us, zt, 2)) - (diff(
#     mus*diff(us, xt), xt) + diff(mus*diff(vs, xt), yt) + 
#    diff(mus*diff(ws, xt), zt) + diff(-Rational(2,3)*mus*divus, xt)) + es*diff(pg, xt) )
#
#vsSrc = ( rops*(diff(us*vs, xt) + diff(vs**2, yt) + diff(vs*ws, zt)) - 
#  mus*(diff(vs, xt, 2) + diff(vs, yt, 2) + diff(vs, zt, 2)) - (diff(
#     mus*diff(us, yt), xt) + diff(mus*diff(vs, yt), yt) + 
#    diff(mus*diff(ws, yt), zt) + diff(-Rational(2,3)*mus*divus, yt)) + es*diff(pg, yt) )
#
#wsSrc = ( rops*(diff(us*ws, xt) + diff(vs*ws, yt) + diff(ws**2, zt)) - 
#  mus*(diff(ws, xt, 2) + diff(ws, yt, 2) + diff(ws, zt, 2)) - (diff(
#     mus*diff(us, zt), xt) + diff(mus*diff(vs, zt), yt) + 
#    diff(mus*diff(ws, zt), zt) + diff(-Rational(2,3)*mus*divus, zt)) + es*diff(pg, zt) )

ugSrc = ( rog*(diff(eg*ug**2, xt) + diff(eg*ug*vg, yt) + diff(eg*ug*wg, zt)) - 
   mug*(diff(ug, xt, 2) + diff(ug, yt, 2) + diff(ug, zt, 2)) - (diff(
      mug*diff(ug, xt), xt) + diff(mug*diff(vg, xt), yt) + 
     diff(mug*diff(wg, xt), zt) + diff(-Rational(2,3)*mug*divug, xt)) + eg*diff(pg, xt) )

vgSrc = ( rog*(diff(eg*ug*vg, xt) + diff(eg*vg**2, yt) + diff(eg*vg*wg, zt)) - 
   mug*(diff(vg, xt, 2) + diff(vg, yt, 2) + diff(vg, zt, 2)) - (diff(
      mug*diff(ug, yt), xt) + diff(mug*diff(vg, yt), yt) + 
     diff(mug*diff(wg, yt), zt) + diff(-Rational(2,3)*mug*divug, yt)) + eg*diff(pg, yt) )

wgSrc = ( rog*(diff(eg*ug*wg, xt) + diff(eg*vg*wg, yt) + diff(eg*wg**2, zt)) - 
  mug*(diff(wg, xt, 2) + diff(wg, yt, 2) + diff(wg, zt, 2)) - (diff(
     mug*diff(ug, zt), xt) + diff(mug*diff(vg, zt), yt) + 
    diff(mug*diff(wg, zt), zt) + diff(-Rational(2,3)*mug*divug, zt)) + eg*diff(pg, zt) )

usSrc = ( ros*(diff(es*us**2, xt) + diff(es*us*vs, yt) + diff(es*us*ws, zt)) - 
  mus*(diff(us, xt, 2) + diff(us, yt, 2) + diff(us, zt, 2)) - (diff(
     mus*diff(us, xt), xt) + diff(mus*diff(vs, xt), yt) + 
    diff(mus*diff(ws, xt), zt) + diff(-Rational(2,3)*mus*divus, xt)) + es*diff(pg, xt) )

vsSrc = ( ros*(diff(es*us*vs, xt) + diff(es*vs**2, yt) + diff(es*vs*ws, zt)) - 
  mus*(diff(vs, xt, 2) + diff(vs, yt, 2) + diff(vs, zt, 2)) - (diff(
     mus*diff(us, yt), xt) + diff(mus*diff(vs, yt), yt) + 
    diff(mus*diff(ws, yt), zt) + diff(-Rational(2,3)*mus*divus, yt)) + es*diff(pg, yt) )

wsSrc = ( ros*(diff(es*us*ws, xt) + diff(es*vs*ws, yt) + diff(es*ws**2, zt)) - 
  mus*(diff(ws, xt, 2) + diff(ws, yt, 2) + diff(ws, zt, 2)) - (diff(
     mus*diff(us, zt), xt) + diff(mus*diff(vs, zt), yt) + 
    diff(mus*diff(ws, zt), zt) + diff(-Rational(2,3)*mus*divus, zt)) + es*diff(pg, zt) )

TgSrc = ( eg*rog*Cpg*(ug*diff(tg, xt) + vg*diff(tg, yt) + wg*diff(tg, zt)) + 
  diff((-kg*diff(tg, xt)), xt) + diff((-kg*diff(tg, yt)), yt) + 
  diff((-kg*diff(tg, zt)), zt) )

TsSrc = ( rops*Cps*(us*diff(ts, xt) + vs*diff(ts, yt) + ws*diff(ts, zt)) + 
  diff((-ks*diff(ts, xt)), xt) + diff((-ks*diff(ts, yt)), yt) + 
  diff((-ks*diff(ts, zt)), zt) )

ThsSrc = ( Rational(3,2)*rops*(diff(us*ths, xt) + diff(vs*ths, yt) + diff(ws*ths, zt)) - 
  diff(ks*diff(ths, xt), xt) - diff(ks*diff(ths, yt), yt) - diff(ks*diff(ths, zt), zt) )

ropgSrc = ( diff(eg*rog*ug, xt) + diff(eg*rog*vg, yt) + diff(eg*rog*wg, zt) )

ropsSrc = ( diff(es*ros*us, xt) + diff(es*ros*vs, yt) + diff(es*ros*ws, zt) )

# print mms source terms
print('\n ****** MMS Source terms ******')

# print('\n pg =')
# print(pg)

print('\n ugSrc =')
print(ugSrc)

print('\n vgSrc =')
print(vgSrc)

print('\n wgSrc =')
print(wgSrc)

print('\n usSrc =')
print(usSrc)

print('\n vsSrc =')
print(vsSrc)

print('\n wsSrc =')
print(wsSrc)

print('\n TgSrc =')
print(TgSrc)

print('\n TsSrc =')
print(TsSrc)

print('\n ThsSrc =')
print(ThsSrc)

print('\n ropgSrc =')
print(simplify(ropgSrc))

print('\n ropsSrc =')
print(simplify(ropsSrc))
