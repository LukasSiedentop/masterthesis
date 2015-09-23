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

# Liest für jeden Wert in der Liste wavelengths [nm] die Diode aus und gibt sie in einer Liste gleicher Länge zurück.
def getSpectrum(ser, wavelengths):
    # Den Monochromator aktivieren und die Erste Wellenlänger erreichen lassen
    if (not ser.isOpen()): ser.open()
    ser.write(str('{:10.2f}'.format(wavelengths[0])) + ' GOTO\r')
    time.sleep(2)
    
    # Liste & Fortschrittsbalken initialisieren
    values = list()
     
    for wavelength in wavelengths:
	   # Monochromator die Wellenlänge anfahren lassen (der Monochromator nimmt das Format '1.23 GOTO\r')
        ser.write(str('{:10.2f}'.format(wavelength)) + ' GOTO\r')
        # Den Monochromator das Ziel erreichen lassen...
        time.sleep(.5)
        # Werte speichern
        value = niusb6211.main()
        values.append(value)
        
        # Fortschrittsbalken
        #print('Meassure Wavelength ' + str(wavelength) + 'nm.')
        #print('wl: ' + str(wavelength) + 'nm, val: ' + str(value) + 'V.\r')
						
    # Alles Licht durchlassen (wavelength = 0), Monochromator deaktivieren
    ser.write(str('0.00 GOTO\r'))
    ser.close()
    #print('Aufnahme bei lambda=' + str(wavelengths[len(wavelengths)-1]) + 'nm beendet.')
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

# Parameter bekommen

try:
    sample = str(raw_input('Name of sample: '))
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

# Motor aktivieren + Probe bewegen
if (not motionSer.isOpen()): motionSer.open()
# Geschwindigkeit setzen
motionSer.write('1VU10;1VA10\r')
# Motor anschalten
motionSer.write('1MO;1MO?\r')
motor = motionSer.readline().rstrip()
if (motor != '1'):
    print('ERROR: Motor cannot be turned on. Check power cords. Exiting.')
    motionSer.close()    
    monoSer.close()
    exit()
# position relativ um +5 ändern
motionSer.write('1PR+5\r')
motionSer.close()

# kurz warten
wait(message="Referencemeassurement starts in ")

# Fahre Wellenlängen für Messung des Referenzspektrums durch
valuesRef = getSpectrum(monoSer, wavelengths)

# Motor aktivieren + Probe zurückbewegen
if (not motionSer.isOpen()): motionSer.open()
motionSer.write('1MO;1MO?\r')
motor = motionSer.readline().rstrip()
if (motor != '1'):
    print('ERROR: Motor cannot be turned on. Check power cords. Exiting.')
    motionSer.close()    
    monoSer.close()
    exit()
# position relativ um -5 ändern
motionSer.write('1PR-5\r')
motionSer.close()

# Datei speichern C:\Users\Hyperion\Desktop\Lukas\Data
s = 'Wavelength\tTransmission\tReference\tSpectrum\n'
filename = 'C:\Users\Hyperion\Desktop\Lukas\Data\\' + sample + '_' + str(wavelengths[0]) + 'nm-' + str(stepsize) + 'nm-' + str(wavelengths[len(wavelengths)-1]) + 'nm_' + time.strftime('%Y-%m-%d_%H-%M-%S') + '.sp'
f = open(filename, 'w')
f.write(s)
i = 0
while (i < len(wavelengths)):
    s = str(wavelengths[i]) + '\t' + str(valuesTrans[i]) + '\t' + str(valuesRef[i]) + '\t' + str(valuesTrans[i]/valuesRef[i]) + '\n'
    f.write(s)
    i+=1
f.close()

print("Data saved as " + filename + ". Enjoy!")
# TODO : gnuplot ausführen und anzeigen lassen