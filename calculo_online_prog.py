import numpy as np # linear algebra
from IPython.display import Latex
from ipywidgets import interact, interact_manual,interactive,widgets,Layout
import ht


# Datos Gasoil
P1,P2 = (np.array([2302.55849485,   -6.61174115]),
 np.array([-6.94642857e-01,  1.02629170e+03]))
viscosidad_f = np.poly1d(P1)
densidad_f = np.poly1d(P2)
Tmedia = 13+273
rho = densidad_f(Tmedia)
mu = viscosidad_f(1/Tmedia)/1000
nu = mu/rho
#covnersion
psipa = 6894.75
#Datos sensor DP. volt a psi
P3 = [0.579,-0.03]
f_ajuste = np.poly1d(P3)
psivolt  = np.poly1d([1/P3[0],-P3[1]/P3[0]])
#datos loop
D_interno_1 = 0.8366*25.4e-3
Longitudes = np.array([6532.7,6532.7+290+1590.7+290,6532.7+(290+1590.7+290)*2,
                       6532.7+(290+1590.7+290)*2+6518.4])/1000
Li = Longitudes[-1]
def f_calcula(b=1):
    mu = viscosidad.value/1000
    #mu = 14e-4
    delta_w = espesor_w_i.value*1e-3
    G = caudal_i.value/1000/60
    Area = D_interno_1**2*np.pi/4
    Area2 = (D_interno_1-2*delta_w)**2*np.pi/4

    U = G/Area
    U2 = G/Area2

    Re = U * D_interno_1/(mu/rho)
    Re2  = (U2 * D_interno_1-2*delta_w)/(mu/rho)
    f = ht.conv_internal.friction_factor(Re)
    f2 = ht.conv_internal.friction_factor(Re2)
    #f = 64/Re
    #f2 = 64/Re2
    Reynolds.value = f'Re={Re:.0f}'
    deltap1 = Li/D_interno_1 * rho/2*U**2*f/psipa
    deltap2 = Li/(D_interno_1-2*delta_w) * rho/2*U2**2*f/psipa



    c = 64; n = 1;
    aux1  = 2*c*rho*Li*(mu/rho)**n*(4*G/np.pi)**(5-n)
    aux2 = (D_interno_1-2*delta_w)**(5-n)
    deltap_f1 = aux1/aux2
    display(Latex(f'$\Delta P_f = $ {deltap1:.3g}'))
    salida_presion.value = f'$\Delta p_0 = $ {deltap1:.2f} psi'
    salida_presion_w.value = f'$\Delta p_w = $ {deltap2:.2f} psi'
    V0,Vw = [psivolt(deltap1), psivolt(deltap2)]
    salida_tension.value = f'$V_0 = $ {V0:.2f} volt'
    salida_tension_w.value = f'$V_w = $ {Vw:.2f} volt'
    return deltap1,deltap2

output_deltap = widgets.HBox(layout={'border': '4px solid black'})
caudal_i = widgets.FloatSlider( value=5,     min=1,     max=25.0,     step=0.1,     
    disabled=False,      orientation='horizontal',     readout=True,
    readout_format='.1f',description='Caudal [lt/min]',style={'description_width': '100px'})


espesor_w_i = widgets.FloatSlider( value=1,     min=1,     max=4,     step=0.1, 
      description='$\delta_w$ [mm]', disabled=False,     orientation='horizontal',     readout=True,
    readout_format='.1f' ,style={'description_width': '100px'})
viscosidad = widgets.FloatText(    value=np.around(mu*1000,3),     description='$\mu[cP]$',
                               disabled=False,step=0.05,readout_format='.2f')

Reynolds = widgets.Label('Re=')
boton_calcula = widgets.Button(description='calcula')
salida_presion = widgets.Label(value='$\Delta p_0=$')
salida_presion.layout=widgets.Layout(width='120px')
#caudal_i.layout = widgets.Layout(description_width='300px')
#espesor_w_i.layout = widgets.Layout(width='120px')
salida_presion_w = widgets.Label(value='$\Delta p_w=$')
salida_presion_w.layout=widgets.Layout(width='120px')

salida_tension = widgets.Label(value='$V_0=$')
salida_tension.layout=widgets.Layout(width='110px')

salida_tension_w = widgets.Label(value='$V_w=$')
salida_tension_w.layout=widgets.Layout(width='110px')

Reynolds.layout=widgets.Layout(width='150px')
boton_calcula.on_click(f_calcula)
parametros_input = widgets.VBox() 
parametros_input.children = ([caudal_i,espesor_w_i,viscosidad])
panel_control = widgets.VBox(layout={'border': '1px solid black'})
panel_salida_Re = widgets.VBox(layout={'border': '1px solid black'})
panel_salida_p = widgets.VBox(layout={'border': '1px solid black'})
panel_salida_v = widgets.VBox(layout={'border': '1px solid black'})
panel_control.children = ([ parametros_input,boton_calcula])
panel_salida_Re.children = ([Reynolds])
panel_salida_p.children = ([ salida_presion,salida_presion_w])
panel_salida_v.children = ([ salida_tension,salida_tension_w])
output_deltap.children = ([panel_control,panel_salida_Re,panel_salida_p,panel_salida_v]) 