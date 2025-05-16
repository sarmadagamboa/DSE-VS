import math
import customtkinter as ck

ck.set_appearance_mode("system")
ck.set_default_color_theme("dark-blue")

root = ck.CTk()
root.geometry("1600x500")

def main():
        
    P = float(P1.get())
    P_up = float(P_up1.get())
    L_TX = float(L_TX1.get())
    L_TX_up = float(L_TX_up1.get())
    f = float(f1.get())
    f_up_down_ratio = float(f_up_down_ratio1.get())
    d = float(d1.get())
    D = float(D1.get())
    orbit_alt = float(orbit_alt1.get())
    destination = destination1.get()
    theta_ES = float(theta_ES1.get())
    e_TX = float(e_TX1.get())
    B_R_up = float(B_R_up1.get())
    P_SW = float(P_SW1.get())
    P_PXS = float(P_PXS1.get())
    P_BPPX = float(P_BPPX1.get())
    D_C = float(D_C1.get())
    T_DL = float(T_DL1.get())
 


    
    if destination != "Earth" and destination != "Earths Moon":
        theta_ES = math.radians(theta_ES)
    
    f = f*10**9

    f_up = f*f_up_down_ratio

    c = 300000000
    pi = math.pi

    ##Downlink:

    ##Transmitter performances
    P_TX = 10*math.log10(P) ##transforming the power from watts to decibels
    L_TX = 10*math.log10(L_TX) ##transforming the loss factor to decibels
    G_TX = 20*math.log10(d) + 20*math.log10(f/10**9) + 17.8 ##finding the gain  of the antenna  in dB
    EIRP = P_TX + L_TX + G_TX ##finding EIRP

    ##Transmitter antenna pointing loss
    G_t = 10**(G_TX/10) ##Finding gain of spacecraft antenna
    Lambda = c/f ##Finding downlink wavelength
    D_t = math.sqrt((G_t*Lambda**2)/(0.55*pi**2))
    Alpha_half = 21/(f*D_t/10**9)
    L_TP = 12*((e_TX/Alpha_half)**2) ##Transmitter antenna pointing loss

    ##Free space loss
    d_E = 149597871000 ##Earth-Sun distance
    if destination == "Earth": ##Earth
        S = 6371*(math.sqrt(((orbit_alt/1000 + 6371)/6371)**2 - (math.cos(math.radians(10)))**2) - math.sin(math.radians(10)))
        S = S*1000 ##Distance in m
        L_s = ((4*pi*S)/Lambda)**2 ##Free space loss
        R_planet = 6371 ##Defining radius for later use
        g = 3.986*10**(14) ##Defining gravitational parameter for later use
    elif destination == "Earths Moon": ##Earths moon
        S = 384400000
        S = 6371*(math.sqrt(((S/1000 + 6371)/6371)**2 - (math.cos(math.radians(10)))**2) - math.sin(math.radians(10)))
        S = S*1000
        L_s = ((4*pi*S)/Lambda)**2
        R_planet = 1737.4
        g = 4.904*10**12
    elif destination == "Mars": ##Mars
        d_S = 1.52*d_E
        S = math.sqrt(d_E**2+d_S**2-2*d_E*d_S*math.cos(theta_ES))
        L_s = ((4*pi*S)/Lambda)**2
        R_planet = 3389.5
        g = 4.282*10**13
    elif destination == "Venus":  ##venus
        d_S = 0.72*d_E
        S = math.sqrt(d_E**2+d_S**2-2*d_E*d_S*math.cos(theta_ES))   
        L_s = ((4*pi*S)/Lambda)**2
        R_planet = 6051.8
        g = 3.248*10**14
    L_s = 10*math.log10(L_s)  ##transforming free space loss to dB

    ##Atmospheric loss
    L_A = 0.035/math.sin(math.radians(10))

    ##Receiver G/T
    G_RX = 20*math.log10(D) + 20*math.log10(f/10**9) + 17.8
    G_r = 10**(G_RX/10)
    GoverT = G_r/135
    GoverT = 10*math.log10(GoverT)

    ##Required data rate
    V_ground = (R_planet/(R_planet + orbit_alt/1000))*math.sqrt(g/(R_planet*1000+orbit_alt)) ##ground velocity
    W_line = 2*orbit_alt*(math.sin(math.radians(P_PXS/120))/math.cos(math.radians(P_PXS/120))) ##width per line 
    LinePerS = V_ground/W_line ##Lines captured per second
    R_G = P_BPPX*LinePerS*(P_SW*60/P_PXS) ##Data generated per second
    T_DL = T_DL/24 ##Downlink time
    D_C = D_C/100 ##Payload duty cycle
    B_R = R_G*(D_C/T_DL) ##Required data rate
    B_R = 10*math.log10(B_R) ##Required data rate in dB

    ##Received SNR
    SNR = EIRP - L_TP - L_s - L_A + GoverT + 228.6 - B_R

    ##Required SNR (provided a req BER = 10^-6) and BPSK modulation
    SNR_req = 10.4

    ##Downlinkg margin
    margin_down = SNR-SNR_req


    ##uplink


    ##Transmitter performances
    P_TX_up = 10*math.log10(P_up)
    L_TX_up = 10*math.log10(L_TX_up)
    EIRP_up = P_TX_up + L_TX_up + G_RX

    ##Transmitter antenna pointing loss
    L_TP_up = L_TP

    ##Free space loss
    d_E = 149597871000
    Lambda_up = c/f_up
    if destination == "Earth":
        S_up = 6371*(math.sqrt(((orbit_alt/1000 + 6371)/6371)**2 - (math.cos(math.radians(10)))**2) - math.sin(math.radians(10)))
        S_up = S_up*1000
        L_s_up = ((4*pi*S_up)/Lambda_up)**2
    elif destination == "Earths Moon":
        S_up = 384400000
        S_up = 6371*(math.sqrt(((S_up/1000 + 6371)/6371)**2 - (math.cos(math.radians(10)))**2) - math.sin(math.radians(10)))
        S_up = S_up*1000
        L_s_up = ((4*pi*S_up)/Lambda_up)**2
    elif destination == "Mars":
        d_S_up = 1.52*d_E
        S_up = math.sqrt(d_E**2+d_S_up**2-2*d_E*d_S_up*math.cos(theta_ES))
        L_s_up = ((4*pi*S_up)/Lambda_up)**2
    elif destination == "Venus":
        d_S_up = 0.72*d_E
        S_up = math.sqrt(d_E**2+d_S_up**2-2*d_E*d_S_up*math.cos(theta_ES))   
        L_s_up = ((4*pi*S_up)/Lambda_up)**2
    L_s_up = 10*math.log10(L_s_up)

    ##Atmospheric loss
    L_A_up = 0.035/math.sin(math.radians(10))

    ##Receiver G/T
    GoverT_up = G_t/135
    GoverT_up = 10*math.log10(GoverT_up)

    ##Required data rate
    B_R_up = 10*math.log10(B_R_up)

    ##Received SNR
    SNR_up = EIRP_up - L_s_up - L_A_up - L_TP_up + GoverT_up + 228.6 - B_R_up

    ##Required SNR (provided a req BER = 10^-6) and BPSK modulation
    SNR_req_up = 10.4

    ##Uplink margin
    margin_up = SNR_up-SNR_req_up




    print('Your EIRP is', EIRP,'for downlink and', EIRP_up, 'for uplink')
    print('Your Pointing loss is', L_TP, 'for downlink and', L_TP_up, 'for uplink')
    print('Your Free space loss is', L_s, 'for downlink and', L_s_up, 'for uplink')
    print('Your atmospheric Loss is', L_A, 'for downlink and', L_A_up, 'for uplink')
    print('Your G/T is', GoverT, 'for downlink and', GoverT_up, 'for uplink')
    print('Your SNR for downlink is',  SNR, 'and for uplink it is',  SNR_up)
    print('Your required SNR is',  SNR_req, 'for downlink and', SNR_req_up, 'for uplink')
    print('Your link margin is', margin_down, 'for downlink and', margin_up, 'for uplink')

import tkinter as tk

# Scrollable container
container = ck.CTkFrame(master=root)
container.pack(pady=0, padx=60, fill="both", expand=True)

canvas = tk.Canvas(container, bg="#2b2b2b", highlightthickness=0)  # match CTk dark background
scrollbar = ck.CTkScrollbar(master=container, orientation="vertical", command=canvas.yview)
scrollbar.pack(side="right", fill="y")

scrollable_frame = ck.CTkFrame(master=canvas)
scrollable_frame.bind(
    "<Configure>",
    lambda e: canvas.configure(
        scrollregion=canvas.bbox("all")
    )
)

canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
canvas.configure(yscrollcommand=scrollbar.set)
canvas.pack(side="left", fill="both", expand=True)

label = ck.CTkLabel(master=scrollable_frame, text="What is your spacecraft transmitter power (in Watts)?")
label.pack(pady=0, padx=0, side="top")
P1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Watts")
P1.pack(pady=0, padx=0)
#print('What is your spacecraft transmitter power (in Watts)?')
#P = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your groundstation transmitter power (in Watts)?")
label.pack(pady=0, padx=0, side="top")
P_up1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Watts")
P_up1.pack(pady=0, padx=0)
#print('What is your groundstation transmitter power (in Watts)?')
#P_up = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your spacecraft transmitter loss factor?")
label.pack(pady=0, padx=0, side="top")
L_TX1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Factor")
L_TX1.pack(pady=0, padx=0)
#print('What is your spacecraft transmitter loss factor?')
#L_TX = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your ground station transmitter loss factor?")
label.pack(pady=0, padx=0, side="top")
L_TX_up1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Factor")
L_TX_up1.pack(pady=0, padx=0)
#print('What is your ground station transmitter loss factor?')
#L_TX_up = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is the frequency of your downlink signal (in GHz)?")
label.pack(pady=0, padx=0, side="top")
f1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="GHz")
f1.pack(pady=0, padx=0)
#print('What is the frequency of your downlink signal (in GHz)?')
#f =  float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your turnaround ratio (in decimal form)?")
label.pack(pady=0, padx=0, side="top")
f_up_down_ratio1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Decimal")
f_up_down_ratio1.pack(pady=0, padx=0)
#print('What is your turnaround ratio (in decimal form)?')
#f_up_down_ratio = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is the diameter of your spacecraft antenna (in meters)?")
label.pack(pady=0, padx=0, side="top")
d1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Meters")
d1.pack(pady=0, padx=0)
#print('What is the diameter of your spacecraft antenna?')
#d = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is the diameter of your receiving antenna (in meters)?")
label.pack(pady=0, padx=0, side="top")
D1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Meters")
D1.pack(pady=0, padx=0)
#print('What is the diameter of your receiving antenna (in m)?')
#D = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="At what altitude above the celestial body does your spacecraft orbit (in meters)?")
label.pack(pady=0, padx=0, side="top")
orbit_alt1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Meters")
orbit_alt1.pack(pady=0, padx=0)
#print('At what altitude above the celestial body does your spacecraft orbit (in m)?')
#orbit_alt = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="Which celestial body is your spacecraft orbiting?")
label.pack(pady=0, padx=0, side="top")
destination1 = ck.CTkComboBox(master = scrollable_frame, values=["Earth", "Earths Moon", "Mars", "Venus"])
destination1.pack(pady=0, padx=0)
#print('Which celestial body is your spacecraft orbiting?')
#print('1. Earth')
#print('2. Earths Moon')
#print('3. Mars')
#print('4. Venus')
#print('Type the number associated with your number')
#destination = int(input('>'))


label = ck.CTkLabel(master=scrollable_frame, text="What is your elongation  angle (in degrees)?")
label.pack(pady=0, padx=0, side="top")
theta_ES1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Degrees")
theta_ES1.pack(pady=0, padx=0)


label = ck.CTkLabel(master=scrollable_frame, text="What is your pointing offset angle (in degrees)?")
label.pack(pady=0, padx=0, side="top")
e_TX1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Degrees")
e_TX1.pack(pady=0, padx=0)
#print('What is your pointing offset angle (in deg)?')
#e_TX = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your required uplink data rate (in bit/s)?")
label.pack(pady=0, padx=0, side="top")
B_R_up1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="bit/s")
B_R_up1.pack(pady=0, padx=0)
# print('What is your required uplink data rate (in bit/s)?')
# B_R_up = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your Payload swath width angle (in degrees)?")
label.pack(pady=0, padx=0, side="top")
P_SW1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Degrees")
P_SW1.pack(pady=0, padx=0)
# print('What is your Payload swath width angle (in deg)?')
# P_SW = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your payload pixel size (in arcmin)?")
label.pack(pady=0, padx=0, side="top")
P_PXS1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Arcmin")
P_PXS1.pack(pady=0, padx=0)
# print('What is your payload pixel size (in arcmin)?')
# P_PXS = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your payload bits per pixel?")
label.pack(pady=0, padx=0, side="top")
P_BPPX1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="bits/pixel")
P_BPPX1.pack(pady=0, padx=0)
# print('What is your payload bits per pixel?')
# P_BPPX = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="What is your payload duty cycle (in percentage)?")
label.pack(pady=0, padx=0, side="top")
D_C1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Percentage")
D_C1.pack(pady=0, padx=0)
# print('What is your payload duty cycle (in percent form)?')
# D_C = float(input('>'))

label = ck.CTkLabel(master=scrollable_frame, text="How many hours per day of downlink time does your payload have?")
label.pack(pady=0, padx=0, side="top")
T_DL1 = ck.CTkEntry(master=scrollable_frame, placeholder_text="Hours")
T_DL1.pack(pady=0, padx=0)
# print('How many hours per day of downlink time does your payload have?')
# T_DL = float(input('>'))

button = ck.CTkButton(master=scrollable_frame, text="Run", command=main)
button.pack(pady=0, padx=0, side="top")

root.mainloop()