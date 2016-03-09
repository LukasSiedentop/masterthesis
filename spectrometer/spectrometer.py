#! ~/Scripts/spectrometer python
# coding: utf8
# -*- coding: utf-8 -*-

# @author: pknappe, dropers, lsiedentop

# COM-Port
import serial
# National Instruments
import niusb6211
# Sleep
import time
# stdout flush
import sys
# gnuplot
import subprocess

# COM-Ports identifizieren
import serial.tools.list_ports

monoPort = ''
motionPort = ''
ports = list(serial.tools.list_ports.comports())
for p in ports:
    if 'ATEN' in p[1]:
        motionPort = p[0]
    if 'Princeton' in p[1]:
        monoPort = p[0]

# gnuplot öffnen
gp = subprocess.Popen(['gnuplot'],
                    stdin=subprocess.PIPE,
                    )

# kopnvertiert die gegebenen Listen, sodass Gnuplot sie in seinem stdin interpretieren kann
def listsToString(wavelengths, data, ref=list()):
    string = ""

    i = 0

    if (len(ref) == 0):
        while (i<len(data)):
            string += str(wavelengths[i]) + " " + str(data[i]) + "\n"
            i+=1
    else :    
        while (i<len(data)):
            string += str(wavelengths[i]) + " " + str(data[i]) + " " + str(ref[i]) + " " + str(data[i]/ref[i]) + "\n"
            i+=1
            
    return string + "e\n"

# Liest für jeden Wert in der Liste wavelengths [nm] die Diode aus und gibt sie in einer Liste gleicher Länge zurück.
def getSpectrum(ser, wavelengths):
    # Den Monochromator aktivieren und die Erste Wellenlänger erreichen lassen
    if (not ser.isOpen()): ser.open()
    ser.write(str('{:10.2f}'.format(wavelengths[0])) + ' GOTO\r')
    time.sleep(0.5)
    
    # gnuplot range anpassen
    gp.stdin.write('set xrange ['+str(wavelengths[0])+':' + str(wavelengths[len(wavelengths)-1]) + ']\n')
    
    # Liste initialisieren
    values = list()
     
    for wavelength in wavelengths:
	   # Monochromator die Wellenlänge anfahren lassen (der Monochromator nimmt das Format '1.23 GOTO\r')
        ser.write(str('{:10.2f}'.format(wavelength)) + ' GOTO\r')
        
        # Den Monochromator das Ziel erreichen lassen...
        time.sleep(.2) #0.5
        
        # Werte speichern
        value = niusb6211.main()
        values.append(value)
        
        # plotten
        gp.stdin.write("plot '-' w l t ''\n")
        gp.stdin.write(listsToString(wavelengths, values))
 

    # Alles Licht durchlassen (wavelength = 0), Monochromator deaktivieren
    ser.write(str('0.00 GOTO\r'))
    ser.close()
    
    return values

# wartet n [s] Sekunden und zählt runter (http://stackoverflow.com/questions/17220128/display-a-countdown-for-the-python-sleep-function)
def wait(n=5, message='Meassurement starts in '):
    i=n
    sys.stdout.write(message)
    sys.stdout.flush()
    while (i>0):
        sys.stdout.write(str(i)+' ')
        sys.stdout.flush()
        time.sleep(1)
        i-=1
    print('')

# bewegt die probe relativ
def move(ser, direction = '+', distance = 4):
    # Motor aktivieren + Probe bewegen
    if (not ser.isOpen()): ser.open()
    # Motor anschalten
    ser.write('1MO;1MO?\r')
    # Geschwindigkeit setzen VU: maximale Geschwindigkeit setzen, VA: Geschwindigkeit setzen
    ser.write('1VU6;1VA6\r')
    motor = ser.readline().rstrip()
    if (motor != '1'):
        print('ERROR: Motor cannot be turned on. Check power cords. Exiting.')
        ser.close()
        return
        #monoSer.close()
        #exit()
    # position relativ um +5 ändern
    ser.write('1PR' + direction + str(distance) + '\r')
    ser.close()

# gibt das minimum von zwei Listen zurück
def findMin(listA, listB):
    minimum = listA[0] 

    for valA in listA:
        if(valA < minimum):
            minimum = valA
    
    for valB in listB:
        if(valB < minimum):
            minimum = valB
    
    return minimum
   
# Parameter bekommen
try:
    sample = str(raw_input('Name of sample (without space): '))
except ValueError:
    sample = "sample"
    print "no string input; sample is named sample"
    
print('Select scanningrange (0nm - 2800nm, min 0.1nm (?)) and resolution (default: 800nm).')
try:
    start__wl = float(raw_input('Start Wavelength [nm]: '))
except ValueError:
    start__wl = 800
    print "800 nm"
    
try:
    end__wl = float(raw_input('End Wavelength [nm] (default: 1700nm): '))
except ValueError:
    end__wl = 1700
    print "1700 nm"
    
try:
    stepsize = float(raw_input('Stepsize [nm] (default: 10 nm): '))
except ValueError:
    stepsize = 10
    print "10 nm"

# Wellenlängen berechnen
wavelengths = list()
# den ersten Messpunkt später wegschmeißen
#wavelengths.append(start__wl-stepsize)
while (start__wl <= end__wl):
	wavelengths.append(start__wl)
	start__wl += stepsize

# Serial-Connection des Monochromators konfigurieren
monoSer = serial.Serial(
	port=monoPort,
	baudrate=9600,
	parity=serial.PARITY_NONE,
	stopbits=serial.STOPBITS_ONE,
	bytesize=serial.EIGHTBITS
)

# Serial-Connection des Motioncontrolers konfigurieren
motionSer = serial.Serial(
	port=motionPort,
	baudrate=19200,
	parity=serial.PARITY_NONE,
	stopbits=serial.STOPBITS_ONE,
	bytesize=serial.EIGHTBITS,
     timeout=1
)

str(raw_input('Please put in the transmission-sample and press Enter to start measuring.'))

# kurz warten
wait()

# Fahre Wellenlängen für Messung des Transmissionsspektrums durch
valuesTrans = getSpectrum(monoSer, wavelengths)

# Probe wegbewegen
move(motionSer, '+', 3)

# kurz warten
wait(message="Referencemeassurement starts in ")

# Fahre Wellenlängen für Messung des Referenzspektrums durch
valuesRef = getSpectrum(monoSer, wavelengths)

# minimum von min(ref) und min(trans) finden und zu allen werten dazuaddieren
minimum = findMin(valuesTrans, valuesRef)

if (minimum < 0): 
    i = 0
    while(i<len(valuesTrans)):
        valuesTrans[i] += abs(minimum)
        valuesRef[i] += abs(minimum)
        i+=1
    
    print('|' + str(minimum) + '| added to all values.')

# Probe wieder hinbewegen
move(motionSer, '-', 3)

# Datei speichern C:\Users\Hyperion\Desktop\Lukas\Data
s = 'Wavelength\tTransmission\tReference\tSpectrum\n'
filename = 'C:\Users\Hyperion\Desktop\Lukas\Data\\' + sample + '_' + str(wavelengths[0]) + 'nm-' + str(stepsize) + 'nm-' + str(wavelengths[len(wavelengths)-1]) + 'nm_' + time.strftime('%Y-%m-%d_%H-%M-%S')
f = open(filename + '.sp', 'w')
f.write(s)
i = 0
while (i < len(wavelengths)):
    #valuesRef[i]==0?'nan', .../...
    s = str(wavelengths[i]) + '\t' + str(valuesTrans[i]) + '\t' + str(valuesRef[i]) + '\t' + str(valuesTrans[i]/valuesRef[i]) + '\n'
    f.write(s)
    i+=1
f.close()

print("Data saved as " + filename + ". Enjoy!")

gp.stdin.write("set terminal push\n")
gp.stdin.write("set terminal svg size 640,480 fname 'Verdana' fsize 14\n")
gp.stdin.write("set output '" + filename + ".svg'\n")

gp.stdin.write('set xrange ['+str(wavelengths[0])+':' + str(wavelengths[len(wavelengths)-1]) + ']; set yrange [0:1]\n')
gp.stdin.write('set ytics 0,0.1,1\n')
gp.stdin.write("set key center bottom\n")

gp.stdin.write("set title '" + sample + "' noenhanced\n")
gp.stdin.write("set label '" + filename + "' at screen 0.0, 0.9 font 'Verdana,10' noenhanced\n")

gp.stdin.write("plot '-' u 1:2 w l t 'sample', '-' u 1:3 w l t 'reference', '-' u 1:4 w l t 'transmission'\n")
gp.stdin.write(listsToString(wavelengths, valuesTrans, valuesRef))
gp.stdin.write(listsToString(wavelengths, valuesTrans, valuesRef))
gp.stdin.write(listsToString(wavelengths, valuesTrans, valuesRef))

gp.stdin.write("set terminal pop\n")
gp.stdin.write("set output\n")
gp.stdin.write("replot\n")

str(raw_input('Press enter to exit.'))
# gnupolot schließen
gp.stdin.write("q\n")

# falls etwas schief gelaufen ist
motionSer.close()
monoSer.close()