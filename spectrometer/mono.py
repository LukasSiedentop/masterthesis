import time
import serial
#import niusb6211
import numpy
#import pylab

#print 'Enter your commands below.\r\nInsert "exit" to leave the application.'


input = 1
condition = 1
#while (condition):
start__wl = 0 # wavelength to start with
end__wl = 0
# get keyboard input
print "wavelength?"
wavelength = raw_input(">> ")

wavelength = numpy.float64(wavelength)

print "wavelength:", wavelength

# configure the serial connections (the parameters differs on the device you are connecting to)
ser = serial.Serial(
	port='COM6',
	baudrate=9600,
	parity=serial.PARITY_NONE,
	stopbits=serial.STOPBITS_ONE,
	bytesize=serial.EIGHTBITS
)
ser.isOpen()
print(not ser.isOpen())
#ser.open()

# Python 3 users
# input = input(">> ")
#if input == 'exit':
#    ser.close()
#    print 'serial-port successfully closed'
#    condition = 0
#else:


# send the character to the device
# (note that I happend a \r\n carriage return and line feed to the characters - this is requested by my device)
ser.write(str("{:10.2f}".format(wavelength)) + ' GOTO\r')
out = ''
# let's wait one second before reading output (let's give device time to answer)
time.sleep(.5)
#ser.write('?nm\r')
while ser.inWaiting() > 0:
    out+= ser.read(1)
print out
ser.close()