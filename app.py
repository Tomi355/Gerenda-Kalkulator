import streamlit as st
import os
import pandas as pd

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
d = h-50 #mm

valasztott_beton = st.sidebar.selectbox("Beton",["C12/15","C16/20","C20/25","C25/30","C30/37","C35/45","C40/50","C45/55","C50/60"])
Beton = valasztott_beton
beton_sor = df_beton[df_beton['Betonminőség'] == Beton].iloc[0]

modulus_tipus = st.sidebar.selectbox("Számítás alapja (E)", ["Ecm","Eceff"])
if "Ecm" in modulus_tipus:
    E_hasznalt = float(str(beton_sor['Ecm']).replace(',', '.'))
    st.info(f"Rövid idejű rugalmassági modulussal számolunk: {E_hasznalt} GPa")
else:
    E_hasznalt = float(str(beton_sor['Eceff']).replace(',', '.'))
    st.info(f"Hosszú idejű (kúszott) rugalmassági modulussal számolunk: {E_hasznalt} GPa")

valasztott_acel = st.sidebar.selectbox("Acél",["B500(MH)","B400(MH)","B240(MH)","B500(HH)","B500(HH)"])
Acel = valasztott_acel
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

fi = st.sidebar.number_input("Vasak átmárője (fi)", value=24)
n = st.sidebar.number_input("Vasak száma (n)", value=7)
kengyel_fi=st.sidebar.number_input("kengyel atmero(fi_kengyel)", value=10)
cc=st.sidebar.number_input("betontakaras", value=20)

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

from sympy import symbols, Eq, solve
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

st.header("Eredmények")
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Repesztő nyomaték (Mcr)", f"{int(Mcr)} kNm")
with col2:
    st.metric("II.F.á. nyomaték (M2)", f"{int(M2)} kNm")
with col3:
    st.metric("Teherbírási nyomaték (Mrd)", f"{int(Mrd)} kNm")

# Logikai ellenőrzés látványosan
if b_min > b:

    st.error(f"❌ NEM FÉR EL! Szükséges szélesség: {b_min} mm")
else:
    st.success("✅ A vasalás elhelyezhető.")

print(Mcr)

print(f"A választott beton: {Beton}")
print(f"Számított Mcr érték: {Mcr}")







