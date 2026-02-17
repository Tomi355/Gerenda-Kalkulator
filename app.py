import streamlit as st
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from sympy import symbols, Eq, solve
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

fajlnev = "VasbetonAnyagtulajdonsagok.xlsx"
df_beton = pd.read_excel(fajlnev, sheet_name='Beton')
df_acel = pd.read_excel(fajlnev, sheet_name='Acel')
df_acel.columns = df_acel.columns.str.strip()
df_beton.columns = df_beton.columns.str.strip()

"------------------------------------------------"
st.title("🏗️ Vasbeton Gerenda Méretező")
st.sidebar.header("Geometria és Anyag")
b = st.sidebar.number_input("Szélesség (b) [mm]", value=300)
h = st.sidebar.number_input("Magasság (h) [mm]", value=500)


valasztott_beton = st.sidebar.selectbox("Beton",["C12/15","C16/20","C20/25","C25/30","C30/37","C35/45","C40/50","C45/55","C50/60"])
Beton = valasztott_beton
st.info(f"A használt betonminőségünk: {Beton}")
beton_sor = df_beton[df_beton['Betonminőség'] == Beton].iloc[0]
    
#------------------------ECM--------------------------------------------------------

modulus_tipus = st.sidebar.selectbox("Rugalmassági Modulus (E)", ["Ecm","Eceff","Ec_custom"])
if "Ecm" in modulus_tipus:
    E_hasznalt = float(str(beton_sor['Ecm']).replace(',', '.'))
    st.info(f"Rövid idejű rugalmassági modulussal számolunk: {E_hasznalt} GPa")
elif "Eceff" in modulus_tipus:
    E_hasznalt = float(str(beton_sor['Eceff']).replace(',', '.'))
    st.info(f"Hosszú idejű (kúszott) rugalmassági modulussal számolunk: {E_hasznalt} GPa")
elif "Ec_custom" in modulus_tipus:
    E_exel = float(str(beton_sor['Ecm']).replace(',', '.'))
    phi = st.sidebar.slider("Kúszási tényező ($\phi$)", min_value=1.5, max_value=3.5, value=2.5, step=0.1)
    E_hasznalt = E_exel/(1+phi)
    st.info(f"Egyedi modulussal számolunk ($\phi$={phi}): **{E_hasznalt:.2f}")
    

#---------------------------------------------------------------------------
valasztott_acel = st.sidebar.selectbox("Acél",["B500(MH)","B400(MH)","B240(MH)","B500(HH)","B500(HH)"])
Acel = valasztott_acel
st.info(f"A használt acélminőségünk: {Acel} GPa")
acel_sor = df_acel[df_acel['Eurocode'] == Acel].iloc[0]
    


fck = beton_sor['fck']     # N/mm2
fctm = beton_sor['fctm']    # N/mm2
Eceff = beton_sor['Eceff']
Ecm = beton_sor['Ecm']
fyk = acel_sor['fyk']      # N/mm2
Es = acel_sor['Es'] 

fcd = fck / 1.5             # Betonszilárdság tervezési értéke
fyd = fyk / 1.15            # Acélszilárdság tervezési értéke

ae = Es/E_hasznalt

fi = st.sidebar.number_input("Vasak átmárője (fi)", value=20)
n = st.sidebar.number_input("Vasak száma (n)", value=4)
kengyel_fi=st.sidebar.number_input("kengyel atmero(fi_kengyel)", value=10)
cc=st.sidebar.number_input("betontakaras", value=20)

d = cc+kengyel_fi+fi/2 #mm

b_min= (2*cc+2*kengyel_fi+(n*fi)+((n-1)*max(20,fi)))
if b_min > b:
    print("❌ HIBA: A vasak NEM férnek el egy sorban!")
else:
    print("✅ Rendben: A vasak kényelmesen elférnek egy sorban.")

#"-------I. fesz. állapot.---------"

As=round(n*(fi**2*3.1415)/4,0)
Ai = (h*b)+(As*ae)-As
Sx = (h*b*h/2)+(As*(ae-1)*d)
xi1 = round(Sx/Ai,0)
I1=round((((b*xi1**3)/3)+((b*(h-xi1)**3)/3)+As*(ae-1)*((d-xi1)**2)),2)
Mcr=round(((fctm*I1)/(h-xi1))/10**6,0) #kNm

#"-------II. fesz. állapot.---------"

xi2 = symbols('xi2')
egyenlet = Eq(xi2, (b*xi2*(xi2/2)+As*ae*d)/((xi2*b)+(As*ae)))
xi2 = solve(egyenlet,xi2)
xi2=round(max(xi2),0) #mm
I2 = ((b*xi2**3)/3)+(As*(ae)*(d-xi2)**2)
epsilon_s = fyd / (Es*1000) # Es-t MPa-ban add meg, pl. 200000
epsilon_c = fcd/ (E_hasznalt*1000)
kappa_c = epsilon_c/ xi2
kappa_s = epsilon_s / (d - xi2)
kappa = min(kappa_c,kappa_s)
M2 = round((E_hasznalt*10**3 * I2 * kappa)/10**6,0)

#"-------III. fesz. állapot.---------"

from sympy import symbols, Eq, solve
xc = symbols('xc')
vetegyenlet = Eq(xc*b*fcd,As*fyd)
xc= solve(vetegyenlet,xc)
xc=round(xc[0],0)
Mrd = xc*b*fcd*(d-(xc/2)) #mm*N
Mrd = round(Mrd/10**6,0) #kNm

# Logikai ellenőrzés látványosan
if b_min > b:
    st.error(f"❌ NEM FÉR EL! Szükséges szélesség: {b_min} mm")
else:
    st.success("✅ A vasalás elhelyezhető.")

Kszi_c = xc/d
Kszi_co = 560/(fyd+700)
if Kszi_c > Kszi_co:
    st.error(f"❌ Az acélbetetétek nem folynak: $\\xi_c$ = {Kszi_c:.3f}")
    st.write(f"Mivel $\\xi_c$ > $\\xi_{{co}}$ ({Kszi_c:.3f} > {Kszi_co}), a keresztmetszet túlvasalt.")
else:
    st.success(f"✅ Az acélbetétek folynak: $\\xi_c$ = {Kszi_c:.3f}")
    st.write(f"Mivel $\\xi_c$ < $\\xi_{{co}}$ ({Kszi_c:.3f} < {Kszi_co}), a szerkezet duktilis.")


st.header("I. Feszültségi állapot")
col1,= st.columns(1)
with col1:
    st.metric("I. fesz állapothoz tartozó inercia", f"{int(I1/10000)} cm4")
    st.metric("I. fesz állapothoz tartózó nyomott zóna magassága (xi1) ", f"{int(xi1)} mm")
    st.metric("Repesztő nyomaték (Mcr)", f"{int(Mcr)} kNm")
#-------------------------
st.header("II. Feszültségi állapot")
col1,= st.columns(1)
with col1:
    st.metric("II. fesz állapothoz tartozó inercia", f"{int(I2/10000)} cm4")
    st.metric("II. fesz állapothoz tartózó nyomott zóna magassága (xi2) ", f"{int(xi2)} mm")
    st.metric("II. fesz állapothoz tartózó nyomaték (M2)", f"{int(M2)} kNm")
#--------------------------
st.header("III. Feszültségi állapot")
col1,= st.columns(1)
with col1:
    st.metric("III. fesz állapothoz tartózó nyomott zóna magassága (xc) ", f"{int(xc)} mm")
    st.metric("Teherbírási nyomaték (Mrd)", f"{int(Mrd)} kNm")

#---------------------Rajzok----------------------------
def rajzol_profi_metszet(b, h, d, x_aktualis, As, allapot_nev):
    fig, ax = plt.subplots(figsize=(6, 9))

# Színek beállítása
    color_beton = '#E0E0E0'  # Középszürke beton
    color_nyomott_keret = 'red'
    color_huzott_keret = 'blue'

# 1. ALAP: A teljes beton szürke kitöltéssel
    ax.add_patch(patches.Rectangle((0, 0), b, h, facecolor=color_beton, edgecolor='none', zorder=1))

# 2. NYOMOTT ZÓNA (xi) - Felső rész jelölése kerettel
# A nyomott zóna a tetejétől indul lefelé
    ax.add_patch(patches.Rectangle((0, h - x_aktualis), b, x_aktualis, 
                                   linewidth=3, edgecolor=color_nyomott_keret, 
                                   facecolor='none', zorder=5, label='Nyomott zóna'))
    
# 3. HÚZOTT ZÓNA - Alsó rész jelölése kék kerettel
    ax.add_patch(patches.Rectangle((0, 0), b, h - x_aktualis, 
                                   linewidth=2, edgecolor=color_huzott_keret, 
                                   facecolor='none', zorder=4,))
    
# 4. ACÉL JELÖLÉSE (As területű téglalap)
    # A téglalap magassága (vastagsága), hogy a területe As legyen a b szélesség mellett:
    b_acel = b - 100
    v_acel = As / b_acel
    y_acel_kozep = h - d # d távolság fentről, tehát h-d alulról

# A fekete téglalap d távolságra fentről (középpontosan illesztve a 'd' vonalra)
    ax.add_patch(patches.Rectangle((50, y_acel_kozep - v_acel/2), b_acel, v_acel, 
                                   facecolor='black', edgecolor='black', zorder=10))

# --- KÓTÁZÁS ---
    # d jelölése fentről
    ax.annotate('', xy=(b + 40, h), xytext=(b + 40, h - d),
                arrowprops=dict(arrowstyle='<|-|>', color='black', lw=1.5))
    ax.text(b + 50, h - d/2, f"d={d}", va='center', fontweight='bold')

    # x jelölése fentről (pirossal)
    ax.annotate('', xy=(b + 100, h), xytext=(b + 100, h - x_aktualis),
                arrowprops=dict(arrowstyle='<|-|>', color='red', lw=2))
    ax.text(b + 110, h - x_aktualis/2, f"x={x_aktualis:.1f}", color='red', va='center', fontweight='bold')

    # h és b jelölése
    ax.annotate('', xy=(-40, 0), xytext=(-40, h), arrowprops=dict(arrowstyle='<|-|>', color='black'))
    ax.text(-50, h/2, f"h={h}", rotation=90, va='center', ha='right')
    
    ax.annotate('', xy=(0, h + 40), xytext=(b, h + 40), arrowprops=dict(arrowstyle='<|-|>', color='black'))
    ax.text(b/2, h + 50, f"b={b}", ha='center', va='bottom')

    # Tengelyvonal jelzése az acélnál
    ax.axhline(y=h-d, color='black', linestyle=':', alpha=0.5, zorder=3)

    # Megjelenítési beállítások
    ax.set_title(allapot_nev, fontsize=14, pad=30)
    ax.set_xlim(-150, b + 200)
    ax.set_ylim(-50, h + 150)
    ax.set_aspect('equal')
    ax.axis('off')
    
    return fig        
st.divider()
st.subheader("📐 Keresztmetszeti ábrák")

# Fülek létrehozása az áttekinthetőségért
tab1, tab2, tab3 = st.tabs(["I. Állapot", "II. Állapot", "III. Állapot"])

with tab1:
    st.write("Repedésmentes keresztmetszet")
    st.pyplot(rajzol_profi_metszet(b, h, d, xi1, As, "I. Feszültségi állapot"))

with tab2:
    st.write("Berepedezett keresztmetszet")
    st.pyplot(rajzol_profi_metszet(b, h, d, xi2, As, "II. Feszültségi állapot"))

with tab3:
    st.write("Teherbírási hatrállapot")
    st.pyplot(rajzol_profi_metszet(b, h, d, xc, As, "III. Feszültségi állapot"))

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8))

#---------------------Diagram-----------------------------------
st.header("📈 Teherbírás optimalizálása")
As_range = np.linspace(100,10000,50) 
col_info1, col_info2 = st.columns(2)
with col_info1:
    st.info(f"**Aktuális vasalás ($A_s$):**\n\n {As:.0f} mm²")
with col_info2:
    if b_min <= b:
        st.success(f"**Elhelyezhetőség:**\n\n ✅ A vasak elférnek 1 sorban\n(Szükséges: {b_min:.0f} mm)")
    else:
        st.error(f"**Elhelyezhetőség:**\n\n ❌ NEM férnek el 1 sorban!\n(Szükséges: {b_min:.0f} mm)")

Mrd_list = []
xi_list = []

for As_temp in As_range:
    xc_temp = (As_temp * fyd) / (b * fcd)
    Mrd_temp = xc_temp * b * fcd * (d - (xc_temp / 2)) / 10**6
    xi_temp = xc_temp / d
    
    Mrd_list.append(Mrd_temp)
    xi_list.append(xi_temp)

fig2, ax2 = plt.subplots(figsize=(8, 5))

for i in range(len(As_range)-1):
    # A Kszi_co a te kódban a határérték (pl. 0.493)
    color = 'green' if xi_list[i] <= Kszi_co else 'red'
    ax2.plot(As_range[i:i+2], Mrd_list[i:i+2], color=color, lw=3)

ax2.scatter([As], [Mrd], color='blue', s=100, zorder=5, label='Jelenlegi vasalás')

# Formázás és feliratozás (LaTeX görög betűkkel)
ax2.set_xlabel("Vasalás területe ($A_s$) [mm²]")
ax2.set_ylabel("Teherbírás ($M_{Rd}$) [kNm]")
ax2.set_title("Teherbírás változása a vasmennyiség függvényében (100-10000 mm²)")
ax2.grid(True, linestyle='--', alpha=0.7)

# Segédvonal a jelenlegi vasmennyiségnél
ax2.axvline(x=As, color='blue', linestyle=':', alpha=0.5)

# Egyedi jelmagyarázat összeállítása
custom_lines = [Line2D([0], [0], color='green', lw=3),
                Line2D([0], [0], color='red', lw=3),
                Line2D([0], [0], color='blue', marker='o', linestyle='None')]
ax2.legend(custom_lines, ['Alulvasalt ($\\xi < \\xi_{co}$)', 'Túlvasalt ($\\xi > \\xi_{co}$)', 'Aktuális pont'])

st.pyplot(fig2)

st.divider() # Egy kis elválasztó vonal a profi kinézetért
st.metric(
    label="Aktuális Teherbírás ($M_{Rd}$)", 
    value=f"{Mrd:.2f} kNm",
    help="Ez a jelenleg beállított vasmennyiséghez tartozó maximális nyomatéki teherbírás."
)