# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 13:26:57 2020

@author: adbp
"""

import math as m
import numpy as np
import matplotlib.pyplot as plt

#Ka*m.exp(-Ka*t)-Ke*m.exp(-Ke*t)      Ke =  0.154(h^-1)   t = 1.5(h)
Ke = 0.154/60
t = 1.5*60
def Ka(x):
    return x*m.exp(-x*t)-Ke*m.exp(-Ke*t)

def Ka_(x):
    return m.exp(-t*x)-t*x*m.exp(-t*x)

def bissKa(a,b):
    counter = 0
    while abs(Ka(a)-Ka(b)) > 0.0000001:
        counter += 1
        m = (a+b)/2
        if Ka(a)*Ka(m) < 0:
            b = m
        else:
            a = m
    return (a,b,counter)


def ropeKa(a,b):
    counter = 0
    w_ = (a*Ka(b) - b*Ka(a)) / (Ka(b)-Ka(a))
    w = 100
    while abs(w - w_) > 0.0000001:
        counter += 1
        w = w_
        if Ka(a)*Ka(w) < 0 :
            b = w
        else:
            a = w
        w_ = (a*Ka(b) - b*Ka(a)) / (Ka(b)-Ka(a))
    return (a,b,counter)

def newtonKa(x):
    counter = 0
    x_ = x - Ka(x)/Ka_(x)
    while abs(x-x_) > 0.0000001:
        counter += 1
        x = x_
        x_ = x - Ka(x)/Ka_(x)
    return x,counter

def zero(x):
    return 0

def drawKa(a,b):
    plt.title("Grafico Ka")
    x = np.linspace(a,b,100)
    f2 = np.vectorize(Ka)
    f3 = np.vectorize(zero)
    plt.xlabel("Tempo (minutos)")
    plt.plot(x,f2(x),label = "Ka")
    plt.plot(x,f3(x),label = "y=0")
    plt.legend()
    plt.show()
    
    
drawKa(0,0.05)
print("Newton Ka =",newtonKa(0.04))
Ka_val = newtonKa(0.04)[0]
Ke_val = 0.154/60
print("Bad Newton Ka = ",newtonKa(0)) # equal to Ke
print("Rope Ka=",ropeKa(0.01,0.04))
print("Bad Rope Ka= ", ropeKa(0,0.01))
print("bissection Ka =",bissKa(0.01,0.04))
print("Bad Bissection Ka =",bissKa(0,0.01))

#diferencial equations solution


def D(x): # função de administração
    if (x > 365*24*60):
        return 0
    elif (x % 480 >= 0 and x % 480 <= 90):
        return (x % 480)*40/4050
    else:
        return 0

def drawD(a,b,n):
    x_list = [] 
    y_list = []
    for i in range(a,b,n):
        x_list.append(i)
        y_list.append(D(i))
        
    plt.title("Fução administrção")
    plt.xlabel("Tempo (minutos)")
    plt.ylabel("Dose administrada (mg)")
    plt.plot(x_list,y_list)
    plt.show()

#VAR_____________VAR_____________VAR_____________VAR_____________VAR_____________VAR_____________
tf = 24*60*7
h = 30
Vap = 3250

print("Função de admnistração 0-24 horas")
drawD(0,tf,h) # tempo é em minutos

def f1(x,y,z): #equação diferencial primeira
    return D(x) - Ka_val*y

def f2(x,y,z): #equação diferencial segunda
    return Ka_val*y-Ke_val*z

x_list = []
y_list = []
z_list = []

def euler(x,y,z,h,xf,makeList = True):
    counter = 0
    while (abs(x - xf) > 0.00001):
        counter += 1
        if makeList:
            x_list.append(x)
            y_list.append(y/Vap)
            z_list.append(z/Vap)
        var1 = f1(x,y,z)
        var2 = f2(x,y,z)
        x += h
        y = y + var1*h
        z = z + var2*h
        
    return y,z,counter

print("grafico Euler compartimento centrar e plasmatico")
def drawE(t,y,z,h,tf): 
    print(euler(t,y,z,h,tf))
    plt.title("Grafico Euler")
    plt.xlabel("Tempo (minutos)")
    plt.ylabel("Concentração farmaco (mg mL^-1)")
    plt.plot(x_list,y_list,label = "C.Central")
    plt.plot(x_list,z_list,label = "C.Plasmático")
    plt.legend()
    plt.show()
h =  3.75
drawE(0,0,0,h,tf)

x_list = []
y_list = []
z_list = []

def rk2(x,y,z,h,xf,makeList = True):
    counter = 0
    while (abs(x - xf) > 0.00001):
        counter += 1
        if makeList:
            x_list.append(x)
            y_list.append(y/Vap)
            z_list.append(z/Vap)
        var1 = f1(x + h/2, y + h/2*f1(x,y,z) , z + h/2*f2(x,y,z))
        var2 = f2(x + h/2, y + h/2*f1(x,y,z) , z + h/2*f2(x,y,z))
        x += h
        y = y + var1*h
        z = z + var2*h
    return y,z,counter
        
print("grafico rk2 compartimento centrar e plasmatico")
def drawRK2(t,y,z,h,tf): 
    print(rk2(t,y,z,h,tf))
    plt.title("Grafico RK2")
    plt.xlabel("Tempo (minutos)")
    plt.ylabel("Concentração farmaco (mg mL^-1)")
    plt.plot(x_list,y_list,label = "C.Central")
    plt.plot(x_list,z_list,label = "C.Plasmático")
    plt.legend()
    plt.show()

h = 7.5
drawRK2(0,0,0,h,tf)


x_list = []
y_list = []
z_list = []
 
def rk4(x,y,z,h,xf,makeList = True):
    while (abs(x - xf) > 0.00001):
        if makeList:
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
        d1f1 = h * f1(x,y,z)
        d1f2 = h * f2(x,y,z)
        d2f1 = h * f1(x + h/2, y + d1f1/2, z + d1f2/2)
        d2f2 = h * f2(x + h/2, y + d1f1/2, z + d1f2/2)
        d3f1 = h * f1(x + h/2, y + d2f1/2, z + d2f2/2)
        d3f2 = h * f2(x + h/2, y + d2f1/2, z + d2f2/2)
        d4f1 = h * f1(x + h, y + d3f1, z + d3f2) 
        d4f2 = h * f2(x + h, y + d3f1, z + d3f2)
        x += h
        y = y + d1f1/6 + d2f1/3 + d3f1/3 + d4f1/6
        z = z + d1f2/6 + d2f2/3 + d3f2/3 + d4f2/6
    return (y,z)
        
#print("grafico rk4 compartimento centrar e plasmatico")      
def drawRK4(t,y,z,h,tf): 
    rk4(t,y,z,h,tf)
    plt.plot(x_list,y_list)
    plt.plot(x_list,z_list)
    plt.show()
h = 60 
#drawRK4(0,0,0,h,tf)    


def Qc(s1,s2,s3):
    qc0 = (s2[0]-s1[0])/(s3[0]-s2[0])
    qc1 = (s2[1]-s1[1])/(s3[1]-s2[1])
    return qc0,qc1

def Erro(s2,s3,ordem):
    e0 = (s3[0]-s2[0])/(2**ordem-1)
    e1 = (s3[1]-s2[1])/(2**ordem-1)
    return e0,e1
  
#QC-VARS____________QC-VARS____________QC-VARS____________QC-VARS____________QC-VARS_____________
h = 15
tf = 24*60*3
s1E = euler(0,0,0,h,tf,False) 
s2E = euler(0,0,0,h/2,tf,False)
s3E = euler(0,0,0,h/4,tf,False)
s4E = euler(0,0,0,h/8,tf,False)
s5E = euler(0,0,0,h/16,tf,False)     

#print("Qc Euler =", Qc(s1E,s2E,s3E))
#print("Qc Euler =", Qc(s2E,s3E,s4E))     
print("Qc Euler =", Qc(s3E,s4E,s5E),s3E[2]) 
print("Erro Euler = ", Erro(s4E,s5E,1))

s1Rk2 = rk2(0,0,0,h,tf,False) 
s2Rk2 = rk2(0,0,0,h/2,tf,False)
s3Rk2 = rk2(0,0,0,h/4,tf,False) 
s4Rk2 = rk2(0,0,0,h/8,tf,False)  

#print("Qc y RK2 =", Qc(s1Rk2,s2Rk2,s3Rk2))
print("Qc y RK2 =", Qc(s2Rk2,s3Rk2,s4Rk2),s2Rk2[2], "h=3.5",s3Rk2[2])
print("Erro RK2 = ", Erro(s3Rk2,s4Rk2,2))

s1Rk4 = rk4(0,0,0,h,tf,False) 
s2Rk4 = rk4(0,0,0,h/2,tf,False)
s3Rk4 = rk4(0,0,0,h/4,tf,False)
s4Rk4 = rk4(0,0,0,h/8,tf,False)  
#s5Rk4 = rk4(0,0,0,h/16,tf,False)
#s6Rk4 = rk4(0,0,0,h/32,tf,False)
#s7Rk4 = rk4(0,0,0,h/64,tf,False)
#s8Rk4 = rk4(0,0,0,h/128,tf,False)
#s9Rk4 = rk4(0,0,0,h/256,tf,False)

print("Qc y RK4 =", Qc(s1Rk4,s2Rk4,s3Rk4))
print("Qc y RK4 =", Qc(s2Rk4,s3Rk4,s4Rk4))
#print("Qc y RK4 =", Qc(s3Rk4,s4Rk4,s5Rk4))
#print("Qc y RK4 =", Qc(s4Rk4,s5Rk4,s6Rk4))
#print("Qc y RK4 =", Qc(s7Rk4,s8Rk4,s9Rk4))



print("Mono compartimental")






   

    