#! ~/Scripts/shapes python
# coding: utf8
# -*- coding: utf8 -*-

# Sinus, Cosinus, pipapo
import numpy as np
import math

# Um zu animieren
import subprocess
import sys

# Fortschrittsbalken
from progressbar import Bar, ProgressBar, Percentage, ETA

########################
##### Berechnungen #####
########################

# Skalarprodukt zweier Vektoren vector1 [µm^2] und vector2 [µm^2]
def scalarProduct(vector1, vector2):
	scalarProduct = 0
	for i in range(3):
		scalarProduct += vector1[i]*vector2[i]
	return scalarProduct

# Bestimmt den Euklidschen Abstand zweier Punkte point1 [µm^3] und point2 [µm^3].
def euklidNorm(point1, point2=np.array([0,0,0])):
	vector = point1 - point2
	return math.sqrt(scalarProduct(vector, vector))

# Gibt den Winkel zwischen der Geraden zwischen point1 [µm^3] und point2 [µm^3] und der z-Achse.
def cuttingAngleZ(point1, point2):
	return math.asin(math.fabs(scalarProduct((point1-point2), np.array([0,0,1]) ))/euklidNorm(point1, point2))

# Gibt die Fläche zurück, die ein um alpha [] um die y-Achse gedrehtes Ellipsoid auf die y-z-Ebene Projiziert.
# Das gedrehts Ellipsoid hat die x-y-z-Halbachsen a [µm], a [µm], gamma*a [µm].
def areaProjectedEllipsoid(alpha, a, gamma):
	# zwecks weniger Rechenaufwand vorrechnen
	sin2 = math.sin(alpha)**2
	cos2 = math.cos(alpha)**2
	ammag2 = 1/gamma**2
	return math.pi * a**2 * 1/math.sqrt(sin2 + ammag2*cos2 - ((1-ammag2)**2 * sin2 * cos2)/(cos2 + ammag2*sin2))

# Gibt das Volumen zurück, das zwischen Punkten die mit einem Ellipsoid verbunden werden auftritt.
# a [µm] - Länge der kurzen Halbachsen in x- und y-Richtung
# gamma [] - Faktor, um den die Halbachse in z-Richtung gestreckt ist, im Vergleich zu a
def pathVolume(points, a, gamma):
	# Punkteliste in Pfadliste umwandeln (an den writes trennen)
	paths = splitList(points)
	
	# Volumen des Ellipsoids
	ellipsoidVolume = math.pi * a**3 * gamma
	
	# Volumensammler
	volume = 0
	
	# Über alle Pfade...
	for path in paths:
		# Über alle Punktindexe außer dem ersten eines Pfades...
		for pointIdx in range(1,len(path)):
			
			# Winkel der Geraden durch die Punkte mit der z-Achse
			alpha = cuttingAngleZ(path[pointIdx-1], path[pointIdx])
			
			# Ellipsenfläche * Pfadlänge - Volumen zweier halb Ellipsoide = Volumen zwischen den Punkten
			# Die zwei Halbellipsoide werden abgezogen, da die gerade das Volumen wäre, das mehrfach gezählt würde
			volume += euklidNorm(path[pointIdx-1], path[pointIdx]) * areaProjectedEllipsoid(alpha, a, gamma) - ellipsoidVolume
	
	# Anzahl der Knoten bestimmen und für jeden das Volumen eines Ellipsoids hinzufügen.
	volume += len(uniquify(points)) * 4/3 * math.pi * a**3 * gamma
	
	return volume

# Gibt das Volumen zurück, das von zwei Eckpunkten einer Box mit Kanten parallel zu den Koordinatenachsen aufgespannt wird.
def boxVolume(box):
	volume = 1
	for i in range(len(box[0])):
		volume *= box[1][i]-box[0][i]
	return volume
	
# Gibt die Pfadlänge einer Punkteliste points [Punkte] zurück.
def pathLength(points):
	# Punkteliste in Pfadliste umwandeln (an den writes trennen)
	paths = splitList(points)
	
	# Distanzensammler
	length = 0
	
	# Über alle Pfade...
	for path in paths:
		# Über alle Punktindexe außer dem ersten eines Pfades...
		for pointIdx in range(1,len(path)):
			# Distanz zum Vorgängerpunkt bestimmen und der Länge hinzufügen
			length += euklidNorm(path[pointIdx-1], path[pointIdx])	
				
	return length

# Schätzt den Füllfaktor einer gegebenen Punktelinie ab.
# a [µm] - Länge der kurzen Halbachsen in x- und y-Richtung
# gamma [] - Faktor, um den die Halbachse in z-Richtung gestreckt ist, im Vergleich zu a
def ff(points, a=0.045, gamma=3):

	pv = pathVolume(points, a=0.09, gamma=3)
	
	ff='Volumen Struktur: ' + str(pv) + 'µm^3'

	bb = getBoundingBox(points)
	bb[0] -= a
	bb[0][2] -= 2*a
	bb[1] += a
	bb[1][2] += 2*a	
	ff+='\nVolumen BB: ' + str(boxVolume(bb)) + 'µm^3'
	
	ff+='\nFüllfaktor: ' + str(100 * pv/boxVolume(bb)) + '%'
	

	return ff

######################
##### Nützliches #####
######################

# Trennt eine Punkteliste an den write-Zeilen auf und gibt eine Liste von Listen (Pfaden) zurück.
# Nach http://stackoverflow.com/questions/4322705/split-a-list-into-nested-lists-on-a-value, ssplit
def splitList(points, splitters=['write']):
	
	seq=list(points)
	
	# wenn splitters in seq überhaupt vorkommt...
	if splitters and seq:
		
		# Ergebnissammler
		result=[]
		
		# Anfangsposition der nächste Unterliste
		begin=0
		
		# gehe über alle Elementpositionen in seq
		for end in range(len(seq)):
			
			# wenn ein Element sowohl in splitters als auch seq vorkommt...
			if seq[end] in splitters:
				
				# wenn die Position dieses Elements nach dem Letzten doppelten ist, Element hinzufügen
				if end > begin:
					result.append(seq[begin:end])
				# neuen Anfang setzen
				begin=end+1
		# wenn kein splitter in seq vorgekommen ist		
		if begin<len(seq):
			result.append(seq[begin:])
		return result
		
	return [seq]

# Gibt die Anzahl der Einzigartigen Punkte zurück
def getNodes(points):
	
	return numberUniques
	
# Gibt die Länge der Liste points [Punkte] ohne doppelte Einträge und ohne Strings zurück.
def uniquify(points):
	# Punkte als Strings ohne writes etc
	pts = list()
	
	# Fortschrittsbalken
	bar = ProgressBar(widgets=['Doppelte Einträge löschen... ', Percentage(), Bar(), ' - ' , ETA()], maxval=len(points)).start()
	i=0
	
	# alle 'write's etc in points löschen und die arrays in strings umwandeln
	for point in points:
		if (type(point) != str):
			pts.append(str(point))
			
		# Fortschrittsbalken
		bar.update(i)
		i += 1
		
	bar.finish()
	
	return set(pts)
	
# Verschiebt die gegebenen Punkte points um den Verschiebungsvektor shifter
def shifting(points, shifter):
	shifted = list()
	for point in points:
		if (type(point)==str): shifted.append(point)
		else: shifted.append(point + shifter)
	
	return shifted

# Rotiert die gegebenen Punkte points um (x,0,0) um den Winkel alpha [].
def rotateX(points, alpha):
	rotatedPoints = list()
	
	for point in points:
		if (type(point)==str):
			rotatedPoints.append(point)
			continue
		
		# wende 2d Drehmatrix an
		y = math.cos(alpha)*point[1] - math.sin(alpha)*point[2]
		z = math.sin(alpha)*point[1] + math.cos(alpha)*point[2]
			
		# schreibe neue Punkte
		rotatedPoints.append(np.array([point[0],y,z]))

	return rotatedPoints	

# Rotiert die gegebenen Punkte points um (0,0,z) um den Winkel alpha [].
def rotateZ(points, alpha):
	rotatedPoints = list()
	
	for point in points:
		if (type(point)==str):
			rotatedPoints.append(point)
			continue
		
		# wende 2d Drehmatrix an
		x = math.cos(alpha)*point[0] - math.sin(alpha)*point[1]
		y = math.sin(alpha)*point[0] + math.cos(alpha)*point[1]
		
		# schreibe neue Punkte
		rotatedPoints.append(np.array([x,y,point[2]]))

	return rotatedPoints

# Gibt die Punkteliste points rückwärts geschrieben wieder aus.
def backward(points):
	rev = list()
	
	for point in reversed(points):
		rev.append(point)
		
	return rev

# Gibt Punkte auf einem Sinus der Amplitude amplitude [µm], Wellenlänge wavelength [µm] und Phase phase [] in der x-z-Ebene zurück.
# Der Wellenzug auf denen die n [Punkte/µm] Punkte verteilt liegen ist length [µm] lang.
def sin(amplitude, wavelength, phase, n, length):
	
	# Abstand (auf x-Achse)
	dl = length/n
	
	# Punkteliste und Laufparameter initialisieren
	sin = list()
	i=0
	
	while (i<n):
		x = i*dl
		sin.append( np.array([x, 0, amplitude*math.sin(phase - x*2*math.pi/wavelength)]) )
		i += 1
	
	return sin

# Gibt eine Punkteliste points [Punkte] gespiegelt an der x-y-Ebene zurück.
def flipXY(points):
	
	for point in points:
		if (type(point)!=str): point[2] *= -1
	return points

# Gibt die Punkteliste points [Punkte] in den Positiven x-y-z-Quadranten verschoben zurück.
def positivate(points):
	# Initialisiere geringste Werte
	xMin = 0
	yMin = 0
	zMin = 0
	for point in points:
		if (type(point)!=str):
			xMin = min(xMin, point[0])
			yMin = min(yMin, point[1])
			zMin = min(zMin, point[2])
	
	return shifting(points, np.array([-xMin, -yMin, -zMin]))

# Skaliert die gegebnenen Punkte points [Punkte] mit dem Faktor a.
def scaleStructure(points, a):
	if (a == 1): return points
	
	scaledPoints = list()
	
	for point in points:
		if (type(point)==str):
			scaledPoints.append(point)
			continue
		
		scaledPoints.append(a * point)
		
	return scaledPoints
	
# Gibt die zwei Extremsten Eckpunkte der Struktur points [Punkte] zurück.
def getBoundingBox(points):
	# Suche ersten Punkt in der Struktur für die Startwerte
	i=0
	while (type(points[i]) == str): i+=1
	
	# Initialisiere Startwerte
	xMin = points[i][0]
	yMin = points[i][1]
	zMin = points[i][2]
	
	xMax = points[i][0]
	yMax = points[i][1]
	zMax = points[i][2]
	
	for point in reversed(points):
		if (type(point)!=str):
			xMin = min(xMin, point[0])
			yMin = min(yMin, point[1])
			zMin = min(zMin, point[2])
			xMax = max(xMax, point[0])
			yMax = max(yMax, point[1])
			zMax = max(zMax, point[2])
	return [np.array([xMin, yMin, zMin]), np.array([xMax, yMax, zMax])]

# Gibt den geometrischen Mittelpunkt einer Punkteliste points [Punkte] zurück.
def getCenter(points):
	bb = getBoundingBox(points)
	return (bb[0]+bb[1])/2

# Bewegt eine Punkteliste points [Punkte] so, dass ihr geometrischer Mittelpunkt bei center [Punkt] liegt
def setCenter(points, center):
	c = getCenter(points)
	return shifting(points, (center - c))
	
# Schreibt den text [] and die Position x, y, z [µm^3] in den Lack.
def writeTextAt(x, y, z, writing):
	# Punkteliste initialisieren
	text = list()
	
	# an die Position gehen
	text.append('TextPositionX ' + str(x))
	text.append('TextPositionY ' + str(y))
	text.append('TextPositionZ ' + str(z))
	
	# text schreiben, chr(34) = "
	text.append('WriteText ' + chr(34) + writing + chr(34))
	
	# zurückgeben
	return text	

# Gibt den header einer .gwl Datei als Liste zurück.
def header():
	header = list()
	# When TimeStampOn each message in the Message Log is marked with the current date and time.
	header.append('TimeStampOn')

	# Sets the current ScanMode to either piezowriting (0) or stagewriting (1).
	# In ScanMode 0 all coordinates are addressed via the piezo movement in x-, y- and z-directions.
	# In ScanMode 1 the x- and y coordinates are addressed via a stage movement whereas the z-coordinate is still addressed via the piezo.
	# Note, that in ScanMode 1 all coordinates are absolute with the origin at the sample center, whereas in ScanMode 0 all coordinates are relative to the current position of the xy-stage on the sample with the origin at the corner of the piezo.
	header.append('ScanMode 0')

	# Each z-coordinate is multiplied with this factor.
	# The DefocusFactor compensates for shifts of the z-coordinate due to the index mismatch between the immersion medium, the substrate and the photoresist.
	header.append('DefocusFactor 1.0')
	
	# Triggers the inversion of the z-axis. To conserve a right-handed coordinate system the x-axis is inverted at the same time.
	# 0 deactivates InvertZAxis
	# 1 activates InvertZAxis
	header.append('InvertZAxis 1')
	
	# Outputs the current position of the microscope z-drive to the Message Log.
	header.append('ZDrivePosition')

	# Displays important parameters in the Message Log.
	header.append('ShowParameter')	

	# Recalibrates the laser power.
	header.append('Recalibrate')
	
	# Switches PerfectShape on.
	#print('ATTENTION! Velocity/PS not set in header.')
	header.append('PerfectShapeQuality')
	
	
	return header

# Setzt die Schreibgeschwindigkeit auf v [µm/s] bei der updateRate [punkte/s]. Die pointDistance wird dann berechnet.
def setVelocity(v=20, updateRate=2500):
	
	pointDistance = int(v/updateRate * 1000) #nm
	
	velocity = list()
	velocity.append('MessageOut ' + chr(34) + 'Schreibgeschwindigkeit v=' + str(v) + 'µm/s (Updaterate ' + str(updateRate) + 'Punkten/s, Punktabstand ' + str(pointDistance) + 'nm).' + chr(34))
		
	# Deactivates PerfectShape
	velocity.append('PerfectShapeOff')	

	# In this mode programmed line segments are exposed continuously.
	# The exposure is paused after each write command.
	# The writing speed is either set via ScanSpeed or via PointDistance and UpdateRate .
	velocity.append('ContinuousMode')
	
	# This command is applicable in ScanMode 0. It interpolates additional points between the programmed coordinates.
	# The distance of the interpolated points is specified via PointDistance.
	# These points are then sent to the piezo with the given UpdateRate.
	velocity.append('ConnectPointsOn')
	
	# Sets the distance between the interpolated points when ConnectpointsOn.
	# Command is only valid in ScanMode 0!
	velocity.append('PointDistance ' + str(pointDistance)) #nm
	
	# Sets the rate with which the programmed or interpolated points are sent to the piezo.
	# Note: This command is only valid in ScanMode 0!
	# The frequency of resonance is at around 100Hz, so stay away from that!
	velocity.append('UpdateRate ' + str(updateRate)) #Hz
	
	return velocity

# Gibt den footer einer .gwl Datei als Liste zurück.
def footer():
	footer = list()
	
	# Displays the typed text on the Message Log.
	footer.append('MessageOut ' + chr(34) + 'Die Struktur ist fertig. So richtig fertig.' + chr(34) )
	
	# Saves the current content of the Message Log in the file specified.	
	footer.append('SaveMessages')
	
	return footer

# Kopiert-Rotiert die gegebene Punkteliste points [Punkte] in alle drei Ebenen
def threePlanes(points):
	threePlanes = list()
	
	# x-y-Ebene
	threePlanes.extend( setCenter(points, np.array([0,0,0])) )
	# x-z-Ebene
	threePlanes.extend( rotateX( setCenter(points, np.array([0,0,0])), (math.pi/2) ) )
	# y-z-Ebene
	threePlanes.extend( rotateZ( rotateX( setCenter(points, np.array([0,0,0])), (math.pi/2) ), (math.pi/2) ) )
	
	return threePlanes
	
###########################
##### Einfache Formen #####
###########################

# Gibt eine Liste von numberPoints [] Punkten pro Umdrehung auf einer Spirale mit Anfangsradius radius [µm].
# Mit jedem Umlauf verändert sich der Radius um dr [µm] und die Höhe um dz [µm].
# Der Parameter n [Umläufe] bestimmt die Anzahl der Umläfue.
# Beispiel: dr = dz = 0µm, n = 1 ergibt einen Kreis.
def spiral(numberPoints, dz, radius, dr, n, leftie):
	# dr pro Winkel berechnen
	dr = dr/numberPoints
	
	# dz pro Winkel
	dz = dz/numberPoints
	
	# Winkel berechnen, Koordinaten & Punkte initialisieren
	alpha = 0
	dAlpha = 2*math.pi/numberPoints
	x = 0
	y = 0
	z = 0
	spiral = list()

	# Letzter Winkel, bisschen weiter wegen Rundungsfehlern. Nötig?
	alphaFinish = n*2*math.pi + dAlpha/4
	
	# links- bzw rechtsrum
	if(leftie):
		dAlpha *= -1
	
	# Punkte berechnen	
	while(math.fabs(alpha) < alphaFinish):
			# neue Koordinaten berechnen
			x = radius*math.cos(alpha)
			y = radius*math.sin(alpha)
			
			# Punkt speichern
			spiral.append(np.array([x,y,z]))
			
			# Parameter weiterlaufen lassen
			alpha += dAlpha
			radius += dr
			z += dz
		
	return spiral

# Gibt einen Bügel der Höhe h [µm] und der Länge l [µm] zurück.
# style 0: =]
# style 1: x]
# style 2: ~~]
# style 3: //	-	funktioniert ohne wand!
def hanger(l, h, style):
	# =]
	if (style == 0):
		hanger = [np.array([0,0,0]), np.array([l,0,0]), np.array([l,0,h]), np.array([0,0,h])]
			
	# x]
	if (style == 1):
		hanger = [np.array([0,0,0]), np.array([l/2, 0, h*1/3]), np.array([l,0,0]), np.array([l,0,h]), np.array([l/2, 0, h*2/3]), np.array([0,0,h])]
			
	# ~~]
	if (style == 2):
		hanger = list()
		# unterer Sinus
		hanger.extend( shifting( sin((h/3), (l/3), math.pi/2, 50, l), np.array([0,0,(-h/3)]) ) )
		# Verbindung
		hanger.extend([np.array([l,0,0]), np.array([l,0,h])])
		# oberer Sinus
		hanger.extend( shifting( backward(sin((h/3), (l/3), -math.pi/2, 50, l)), np.array([0,0,(h*4/3)]) ) )
	
	# //
	if (style == 3):
		d = 25
		hanger = [np.array([0,(-d/2),0]), np.array([l,(-d/2),h]), np.array([l,(d/2),h]), np.array([0,(d/2),0])]
	
	hanger.append('write')
		
	return hanger

# Verbindet zwei Bügel mit zwei Streben. Gibt die zwei Streben zurück.
# style 0: =]
# style 1: x]
# style 2: ~~]
# style 3: //	-	funktioniert ohne wand und ohne Streben!
def bars(hanger1, hanger2, style):
	# Index unterer Eckpunkt
	# =]
	if (style==0):
		lowerIndex = 1
	# x]
	if (style==1):
		lowerIndex = 2
	# ~~]	
	if (style==2):
		lowerIndex = 50
	# //
	if (style == 3):
		return []
		
	# Index oberer Eckpunkt
	upperIndex = lowerIndex + 1
	
	# Punkteliste initialisieren
	bars = list()
	
	# Strebe oben
	bars.append(hanger1[upperIndex])
	bars.append(hanger2[upperIndex])
	bars.append('write')
	
	# Strebe unten
	bars.append(hanger2[lowerIndex])
	bars.append(hanger1[lowerIndex])
	bars.append('write')
	
	return bars

# Gibt eine Liste von Punkten zurück, die ein Gitter der Gitterkonstante (=Strichabstand) g [µm/Strich] der Breite w [µm] bilden.
# Die Striche sind l [µm] lang, auf der Höhe z [µm] und parallel zur x-Achse. Sie werden atniparallel geschrieben.
def grating(w, l, g, z):
	# Anzahl der Striche
	n = math.floor(w/g)
	
	# Punktliste und Laufparameter initialisieren
	grating = list()
	i=0
	
	while(i<n):
		grating.append(np.array([0,i*g,z]))
		grating.append(np.array([l,i*g,z]))
		grating.append('write')
		i += 1
		grating.append(np.array([l,i*g,z]))
		grating.append(np.array([0,i*g,z]))
		grating.append('write')
		i += 1
	
	return grating

# Gibt einen Pfahlwald der Größe w x l x h [µm^3] zurück auf dem Strukturen geschrieben werden können um die Schrumpfung zu messen.
# Die Pfähle haben einen Abstand d [µm] zu den nächsten vier Nachbarn.
def stakes(w, l, h, d):
	# Anzahl Pfosten in x- bzw y-Richtung, abgerundet
	nx = math.floor(w/d)
	ny = math.floor(l/d)
	
	# Punkteliste und Laufparameter initialisieren
	stakes = list()
	stakes.append('%stakes with w x l x h = ' + str(w) + ' x ' + str(l) + ' x ' + str(h) + ' um^3 and a distance of ' + str(d) + 'um.')	
	ix = 0
	iy = 0
	
	# x-Reihe
	xRow = list()
	while (ix <= nx):
		xRow.extend(shifting([np.array([0,0,0]), np.array([0,0,h]), 'write'], np.array([(d*ix),0,0])))
		ix += 1
		
	# y-Reihe
	while (iy <= ny):
		stakes.extend(shifting(xRow, np.array([0,(d*iy),0])))
		iy += 1
	
	return stakes

# Genau wie Stakes, nur das diese oben mit einem Gitter verbunden sind das mit der Intensität power [%] geschrieben ist.
def connectedStakes(w, l, h, d, powerStakes, powerPile):
	# Die Stelzen schauen ue [µm] ins Woodpile hinein
	ue = 0.5
	connectedStakes = list()
	# Stelzen schreiben
	connectedStakes.append('LaserPower ' + str(powerStakes))
	connectedStakes.extend(stakes(w, l, h, d))
	# Verbindungsgitter mit power schreiben
	connectedStakes.append('LaserPower ' + str(powerPile))
	connectedStakes.extend( grating(w, l, d, (h-ue)) )
	connectedStakes.extend( shifting( rotateZ( grating(l, w, d, (h-ue)), (-math.pi/2) ),  np.array([0,l,0])) )
	
	return connectedStakes
	
###########################
##### Komplexe Formen #####
###########################

# Gibt eine Liste von nC [] Eckpunkten pro Umlauf in einem runden Wall.
# der Höhe h [µm] mit Innenradius innerRadius [µm] und Außenradius outerRadius [µm] zurück.
def roundWall(nC, h, innerRadius, outerRadius):
	### optionale Parameter ###
	# Der z-Versatz sollte dem axialen Fokus entsprechen um eine solide Wand zu bekommen.
	dz = 0.6 #µm. 
	# Die Radiuserhöhung sollte dem lateralen Fokus entsprechen um eine solide Wand zu bekommen.
	dr = 0.2 #µm.
	###########################
	
	# Berechne die Anzahl der Umdrehungen der Spirale um den Außenradius zu bekommen.
	n = math.floor((outerRadius - innerRadius)/dr)
	
	# Erstelle erste Disk
	disk = spiral(nC, 0, innerRadius, dr, n, 0)
		
	# Zweite Disk ist die erste rückwärts, nach oben verschoben
	disk.extend( shifting( backward(spiral(nC, 0, innerRadius, dr, n, 1)), np.array([0,0,dz]) ) )
	
	# initialisiere Punkteliste
	roundWall = list()
	roundWall.append('%roundWall with ' + str(nC) + ' Points, ' + str(h) + 'um high, ' + str(innerRadius) + 'um inner radius, ' + str(outerRadius) + 'um outer radius, dr = ' + str(dr) + 'um, dz = ' + str(dz) + 'um.')
	roundWall.extend(disk)
	
	# Starthöhe
	z = 2*dz
	
	while (z<h):
		roundWall.extend(shifting(disk, np.array([0,0,z])))
		z+=(dz*2)
		
	# Füge 'write' Befehl alle nWrite Schritte hinzu
	nWrite = 100
	i = nWrite
	
	while( i < len(roundWall) ):
		roundWall.insert(i, 'write')
		# damit nicht eine Linie übersprungen wird den letzten Punkt nach das write kopieren
		roundWall.insert(i+1, roundWall[i-1])
		i += nWrite
		
	# Abschließendes 'write'
	roundWall.insert(i, 'write')
	
	return roundWall

# Genau wie roundWall, gibt aber eine Hülle zurück statt einer soliden runden Wand.
def hull(nC, h, innerRadius, outerRadius):
	### optionale Parameter ###
	# Der z-Versatz sollte dem axialen Fokus entsprechen um eine solide Wand zu bekommen.
	dz = 0.4 #µm. 
	# Die Radiuserhöhung sollte dem lateralen Fokus entsprechen um eine solide Wand zu bekommen.
	dr = 0.2 #µm.
	# Überlapp um sicher zu sein das die Teile verbunden sind.
	ue = 1 #µm
	###########################
	
	# Anzahl Umdrehungen bis Höhe h [µm] erreicht ist
	nz = math.ceil((h+ue)/dz)
	
	# Anzahl Umdrehungen bis Disk geschrieben ist
	nr = math.ceil(((outerRadius+ue) - (innerRadius-ue))/dr)
	
	# Punkteliste initialisieren
	hull = list()
	hull.append('%hull with ' + str(nC) + ' Points, ' + str(h) + 'um high, ' + str(innerRadius) + 'um inner radius, ' + str(outerRadius) + 'um outer radius, dr = ' + str(dr) + 'um, dz = ' + str(dz) + 'um, ue = ' + str(ue) + 'um overlapp.')
	
	# Innenwand
	hull.extend(spiral(nC, dz, innerRadius, 0, nz, 0))
	hull.append('write')
	
	# Außenwand
	hull.extend(spiral(nC, dz, outerRadius, 0, nz, 0))
	hull.append('write')
	
	# Deckel
	hull.extend(shifting(spiral(nC, 0, (innerRadius-ue), dr, nr, 0), np.array([0,0,h])))
	hull.append('write')
	
	return hull

# Baut das Gerüst mit n [] Haltern der Größe l [µm] x h [µm] um die Struktur in der Mitte zu halten.
# Das Gerüst hat den Innenradius innerRadius [µm] auf der Höhe z [µm] gemessen am unteren Rand des Gerüsts.
# Für den Parameter style siehe Methode hanger.
def frameWork(l, h, n, innerRadius, z, hangerStyle):
	# Berechne Drehwinkel für die Bügel
	alpha = 2*math.pi/n
	
	# Initialisiere Gerüst, baue ersten Bügel und verschiebe ihn
	hanger1 = shifting(hanger(l,h,hangerStyle), np.array([(-innerRadius-l),0,z]))
	framework = list()
	framework.append('%frameWork with l = ' + str(l) + 'um, h = ' + str(h) + 'um, n = ' + str(n) + ' hangers, ' + str(innerRadius) + 'um inner radius, ' + str(z) + 'um above x-y-plane, hangerStyle ' + str(hangerStyle) + '.')
	framework.extend(hanger1)
	
	hangeri = list()
	
	i=1
	while(i<n):
		# Bügel um (0,0,z) drehen
		hangeri = rotateZ(hanger1, alpha*i)

		# Hinzufügen
		framework.extend(hangeri)
		
		# Verbinden	
		hangerOld = rotateZ(hanger1, alpha*(i-1))
		framework.extend(bars(hangerOld, hangeri, hangerStyle))
		
		# weiter gehts
		i += 1
		
	# letzte Verbindung
	framework.extend(bars(hangeri, hanger1, hangerStyle))
	
	return framework

# Gibt eine Woodpilestruktur der Höhe h [µm], Länge l [µm] und Breite w [µm] mit der Gitterkonstanten g [µm/Strich] zurück.
# Der Abstand zweier Ebenen dh [µm] sollte maximal dem axialen Fokus entsprechen um die Ebenen verbinden zu können.
def woodpile(w, l, h, dh, g):
	# Erstes Gitter bauen und den Mittelpunkt nach [0,0,0] setzen zwecks leichterem drehen
	grat = shifting(grating(w, l, g, 0), np.array([-w/2, -l/2, 0]))
	
	# Punkteliste und Laufparameter initialisieren
	woodpile = list()
	woodpile.append('%woodpile g = ' + str(g) + 'um per line, w = ' + str(w) + 'um wide, l = ' + str(l) + 'um long, h = ' + str(h) + 'um high.')
	z=0
	
	while(z<h):
		woodpile.extend( shifting(grat, np.array([0, 0, z])) )
		z += dh
		# negative Rotation damit Anfangs- und Endpunkte übereinanderliegen
		woodpile.extend( shifting( rotateZ(grat, -math.pi/2), np.array([0, 0, z]) ) )
		z += dh
		# Rotation damit Anfangs- und Endpunkte übereinanderliegen
		woodpile.extend( shifting( rotateZ(grat, math.pi), np.array([0, g/2, z]) ) )
		z += dh
		
		woodpile.extend( shifting( rotateZ(grat, math.pi/2), np.array([g/2, 0, z]) ) )
		z += dh
		
	# wieder in den positiven Bereich schieben, damit NanoScribe drauf klar kommt
	woodpile = shifting(woodpile, np.array([w/2, l/2, 0]))
	
	return woodpile

# Gibt ein Woodpile auf Stelzen zurück. Die Stelzen sind oben mit einem Maschengitter verbunden, das mit der gleichen Intensität wie das woodpile geschrieben ist.
# powerStakes [%] - Laserintensität Stelzen
# powerPile [%] - Laserintensität woodpile
def woodpileOnStakes(w, l, hStakes, dStakes, hPile, dhPile, gPile, powerStakes, powerPile):
	wpos = list()
	
	# Stelzen
	wpos.extend( connectedStakes(w, l, hStakes, dStakes, powerStakes, powerPile) )
	
	# Woodpile
	wpos.extend( shifting( woodpile(w, l, hPile, dhPile, gPile), np.array([0,0,(hStakes-(1-0.39))])) )
	
	return wpos

# Gibt einen Trampolin-halter für eine Scheibenstruktur der Höhe h [µm] x Durchmesser d [µm]. Die Struktur befindet sich z [µm] über dem Boden.
# Für den Parameter hangerStyle siehe Methode hanger.
# wallStyle 0 -> solide Wand
# wallStyle 1 -> hohle Wand
def spongeHolder(h, d, z, hangerStyle, wallStyle, structureFile):
	### optionale Parameter ###
	t = 10#µm Dicke der Wand
	l = 20#µm Länge der Seile
	ue = 1#µm Überlapp
	trampolinPower = 40#% Laserpower für den Trampolin
	structurePower = 20.5#% Laserpower für die Struktur
	###########################
	
	# Radius
	r = d/2
	
	# Punkteliste initialisieren
	trampolin = list()
	trampolin.extend(header())

	print('Halterstruktur erstellen...')
	trampolin.append('%Trampolin for a Structure with diameter d=' + str(d) + 'um and height h=' + str(h) + 'um, z=' +str(z)+'um above the interface.')
	trampolin.append('LaserPower ' + str(trampolinPower))
	# Wall
	if (wallStyle == 1):
		trampolin.extend(hull(25, (z + h + 2*ue), (r + l), (r + l + t)))
	else: 
		trampolin.extend(roundWall(25, (z + h + 2*ue), (r + l), (r + l + t)))
	# Gerüst
	trampolin.extend(frameWork((l + (2*ue)),(h - 2*ue), 25, (r-ue), (z + ue), hangerStyle))
	trampolin.append('MessageOut ' + chr(34) + 'Trampolin mit P=' + str(trampolinPower) + 'pr geschrieben.' + chr(34))
	
	# Ins positive verschieben
	trampolin = positivate(trampolin)
	print('Halterstruktur erstellt!')
	
	# Mittelpunkt bestimmen
	centerTrampolin = getCenter(trampolin)
	centerTrampolin[2] = z+h/2
		
	# Die Scheibe in den Mittelpunkt legen
	print('Struktur verschieben...')
	trampolin.append('%HPU with diameter d=240um and height h=16um.')
	trampolin.append('LaserPower ' + str(structurePower))
	trampolin.extend(setVelocity(1))
	trampolin.extend( setCenter( readPoints(structureFile), centerTrampolin ) )
	trampolin.append('MessageOut ' + chr(34) + 'Struktur mit P=' + str(structurePower) + 'pr geschrieben.' + chr(34))
	print('Struktur verschoben!')
	
	# Abschluss
	trampolin.extend(footer())
	
	return trampolin

# Ergibt das Zurr-Zeug um eine Scheibenstruktur der Höhe h [µm] x Durchmesser d [µm] am Boden zu halten.
# Die Seile zeigen Radial l [µm] von der Scheibe weg . Eine gegebene Datei am Pfad structureFile wird eingebaut.
def spongeLasher(h, d, l, structureFile):
	### optionale Parameter ###
	shrinkage = 15 #% Schrumpfung der Struktur für hangerStyle 3
	ue = 1 #µm Überlapp
	n = 20 #lashes Anzahl der Leinen
	###########################
	
	# Skalierungsfaktor
	a = 1
	
	h *= a
	d *= a
	l *= a
	
	# Radius
	r = d/2	
				
	# Höhe z zum festzurren bestimmen
	z = math.sqrt( h**2 + l*d*shrinkage/100 + (shrinkage/100)**2 ) + 5.0/a
	
	# Punkteliste initialisieren
	lashes = list()
	lashes.extend(header())

	# Schutzwall
	print('Schutzwall erstellen...')
	lashes.append('LaserPower 40')
	lashes.append('PerfectShapeFast')
	lashes.extend(roundWall(20, (z+2*h), (r+l), (r+l+10)))
	lashes.append('MessageOut ' + chr(34) + 'Schutzwall mit P=40pr geschrieben.' + chr(34))
	print('Schutzwall erstellt!')

	# Leinen hinzufügen
	print('Halterstruktur erstellen...')
	lashes.append('%Lashes of length l=' + str(l) + 'um for a Structure with diameter d=' + str(d) + 'um and height h=' + str(h) + 'um, scale = ' + str(a) + '.')
	lashes.append('%The Structure is z=' + str('%.1f' % z) + 'um above the interface.')
	lashes.extend(frameWork(l, (z+h-ue), n, (r-ue), 0, 3))
	lashes.append('MessageOut ' + chr(34) + 'Schnuere mit P=40pr geschrieben.' + chr(34))
	print('Halterstruktur erstellt!')
	
	# Mittelpunkt bestimmen
	centerLashes = getCenter(lashes)
	centerLashes[2] = z+h/2
		
	# Die Scheibe in den Mittelpunkt legen
	print('Struktur erstellen...')
	lashes.append('%HPU with diameter d=' + str(d) + 'um and height h=' + str(h) + 'um, scale = ' + str(a) + '.')
	##########################################################################
	lashes.append('LaserPower 19.75')		##### LaserPower			######
	lashes.append('PerfectShapeQuality')	##### PerfectShape/Velocity	######
	##########################################################################
	#lashes.extend(setVelocity(15))
	lashes.extend( setCenter( scaleStructure(readPoints(structureFile), a), centerLashes ) )
	lashes.append('MessageOut ' + chr(34) + 'HPU mit P=19.75pr geschrieben.' + chr(34))
	print('Struktur erstellt!')
	
	# Abschluss
	lashes.extend(footer())
	
	return positivate(lashes)

# Gibt viele Strukturen zurück die in den Ofen gesteckt werden können.
# quadrat voll, quadrat leer, quadrat voll mauer, quadrat leer mauer, woodpile 45° gedreht
# Kreis voll, Kreis leer, Kreis voll mauer, Kreis leer mauer, woodpile
def shrinker():
	
	shapes = list()
	shapes.extend( shifting( roundWall(4, 40, 0, 20), np.array([20, 100, 0]) ) )
	shapes.extend( shifting( hull(4, 40, 0, 20), np.array([70, 100, 0]) ) )
	shapes.extend( shifting( roundWall(4, 40, 10, 20), np.array([120, 100, 0]) ) )
	shapes.extend( shifting( hull(4, 40, 10, 20), np.array([170, 100, 0]) ) )
	shapes.extend( shifting( rotateZ(woodpile(40,40,40, 0.39,1), math.pi/4), np.array([230, 70, 0]) ) )
	
	shapes.extend( shifting( roundWall(25, 40, 0, 20), np.array([20, 20, 0]) ) )
	shapes.extend( shifting( hull(25, 40, 0, 20), np.array([70, 20, 0]) ) )
	shapes.extend( shifting( roundWall(25, 40, 10, 20), np.array([120, 20, 0]) ) )
	shapes.extend( shifting( hull(25, 40, 10, 20), np.array([170, 20, 0]) ) )
	shapes.extend( shifting( woodpile(40,40,40, 0.39, 1), np.array([200, 0, 0]) ) )
	
	return shapes

# Gibt Woodpiles auf Stelzen zurück, geschrieben mit den Laser-Intensitäten von Pmin [%] bis Pmax [%] mit dem Unterschied dP [%].
# Schreibt die Strukturen mit dem Piezo, versetzt diese mit der Stage.
def shrinkerOnStakes(Pmin, dP, Pmax):
	
	# Intensitäten
	P = Pmin
	Pstakes = 40
	
	# Breite und Länge einer Struktur
	w = 30
	# Abstand untereinander
	d = 30
	# Länge der Stakes
	l = 15
	
	# Anzahl Strukturen
	ny = math.floor((Pmax-Pmin)/dP)
	iy = 0
	
	# Punkteliste initialisieren
	shrinker2 = list()
	shrinker2.extend(header())
	
	
	while (iy <= ny):
		shrinker2.append('MessageOut ' + chr(34) + 'Schreibe mit P=' + str(P) + 'pr...' + chr(34))
		# scale1
		#scaleStart = 0.4
		#scaleEnd = 0.9
		#dScale = 0.1
		#scale = scaleStart
		nx = 3#math.floor((scaleEnd - scaleStart)/dScale)
		ix = 0
				
		while (ix <= nx):	
			# Stelzen und Struktur richtig verschoben
			conStakes = connectedStakes(w,w,l,2, Pstakes, P)
			center = getCenter(conStakes)
			center[2] = 4.5+l
			struktur = setCenter(readPoints('./Hyperuniformstrukturen/hpu_30x30x10/hpu_30x30x10_scal0.6'), center)
			
			# Hinzufügen (Laserpower wird in conStakes geschrieben!)
			shrinker2.extend( conStakes )
			shrinker2.append('%HPU-Struktur hpu_30x30x10_scal0.6')
			if (ix != 0):
				shrinker2.extend(setVelocity(v=(ix*5), updateRate=2500))
			shrinker2.extend( struktur )
			shrinker2.append('MessageOut ' + chr(34) + 'Struktur hpu_30x30x10_scal0.6 geschrieben.' + chr(34))
			
			# Die Stage in x-Richtung verschieben.
			shrinker2.append('MoveStageX ' + str(w+d))
			
			# nächster Scale
			#scale += dScale
			ix += 1
			
		shrinker2.append('MessageOut ' + chr(34) + 'Strukturen mit P=' + str(P) + 'pr geschrieben.' + chr(34))
		
		# Stage zurück zu x=0 fahren
		shrinker2.append('MoveStageX ' + str(- ix * (w+d)))
		# Die Stage in y-Richtung verschieben.
		shrinker2.append('MoveStageY ' + str(w+d))
		
		# weiter gehts
		P += dP
		iy += 1
		
	shrinker2.extend(footer())
	
	return shrinker2
	
# miniscript zum testen von Offsetgeschichten
def offsetter(Pmin, dP, Pmax):
	offsetter = list()
	
	# Intensitäten
	P = Pmin
	Plines = 40
	
	# Länge einer Linie
	l = 10
	# Abstand untereinander
	d = 20
	
	# Anzahl Linien in y-Richtung
	n = math.floor((Pmax-Pmin)/dP)
	i = 0
	
	while (i <= n):
		# Woodpile auf Stelzen
		offsetter.append('XOffset 0')
		offsetter.extend( [np.array([0,0,0]), np.array([0,l,0]), 'write'] )
		
		offsetter.append('ZOffset 1')
		offsetter.extend( [np.array([1,0,0]), np.array([1,l,0]), 'write'] )
		offsetter.append('ZOffset 0')
		
		# Schrift neben Strukturen
		#offsetter.extend( writeTextAt(((l+d)*3), (l/2), 3, '<- LP ' + str(P) + 'pr') )
				
		offsetter.append('MessageOut ' + chr(34) + 'Struktur mit P=' + str(P) + 'pr geschrieben bei' + chr(34))
		offsetter.append('ZDrivePosition')
		
		# Die Stage in y-Richtung verschieben. ACHTUNG: nichts in der while-schleife muss also in y-Richtung verschoben werden!
		offsetter.append('MoveStageY ' + str(l))
		
		# weiter gehts
		P += dP
		i += 1
		
	offsetter.append('SaveMessages')
	
	return offsetter

# Gibt drei Spiralen in jeder Ebene zurück. Soll den Piezo testen.
def spiralLogger():
	# LogMode
	# This mode is equivalent to ContinuousMode . Additionally one log-file per line segment is created displaying the following 8 columns:
	# programmed coordinates [μm]: x, y, z, LaserPower [%], effective coordinates [μm]: x, y, z, AOM voltage [V].
	
	'''
	### Erstelle eine Spirale ###
	# Spiralenparameter
	innerRadius = 0.05#µm
	outerRadius = 0.8#µm
	nE = 4# Ecken
	n = 2# Umläufe
	testStructure = list()
	testStructure.extend(spiral(nE, 0, innerRadius, ((outerRadius-innerRadius)/n), n, 0))
	testStructure.append('write')
	#############################
	'''
	
	spiralLogger = list()
	spiralLogger.extend(header())
		
	# normal
	spiralLogger.append('LaserPower 10')
	spiralLogger.append('MessageOut ' + chr(34) + 'Normale Struktur' + chr(34))
	spiralLogger.append('LogMode')
	spiralLogger.extend(threePlanes(readPoints('hexagonLP')))
	
	# erweitert
	spiralLogger.append('MessageOut ' + chr(34) + 'Erweiterte Struktur' + chr(34))
	spiralLogger.append('LogMode')
	spiralLogger.extend(threePlanes(readPoints('hexagonLP')))
	
	spiralLogger.extend(footer())
	
	return shifting(positivate(spiralLogger), np.array([1,1,1]))
	
###########################
##### Nachbearbeitung #####
###########################
	
# Schreibt die gegebenen Punkte points in die Datei filename.
# mode:	'w': file neu erstellen/überschreiben,
#		'a': points anfügen,
#		'x': gibt einen Fehler wenn die Datei schon existiert.
def writePoints(points, filename='./shapesOut', mode='w'):
	# Fortschrittsbalken
	bar = ProgressBar(widgets=['Schreibe Datei ' + filename + ' ', Percentage(), Bar(), ' - ' , ETA()], maxval=len(points)).start()
	i=0
	
	#print('Datei schreiben...')
	f = open(filename, mode)

	# Absatz wenn Punkte angefügt werden sollen
	if(mode == 'a'):f.write('write\n')
	
	for point in points:
		if (type(point)==str): 
			f.write(point + '\n')
		else: 
			f.write('\t'.join(str('%.3f' % coord) for coord in point))
			f.write('\n')
			
		# Fortschrittsbalken
		bar.update(i)
		i += 1

	f.close()
	bar.finish()

# Liest eine .gwl Datei ein. Strings bleiben Strings und leerzeichengetrennte Zahlentripletts werden zu np.array([x,y,z]); Spaltenweise gesehen.
def readPoints(filename):
	# Punkteliste, lineList und Line initialisieren
	points = list()
	lineList = list()
		
	# Zeilen in Liste einlesen
	f = open(filename, 'r')
	lineList = f.read().splitlines()
	f.close()
	
	# Fortschrittsbalken
	bar = ProgressBar(widgets=['Datei ' + filename + ' einlesen ', Percentage(), Bar(), ' - ' , ETA()], maxval=len(lineList)).start()
	i=0
	
	# Strings zu Arrays machen #TODO: geht schöner...
	for line in lineList:
		point = line.split('\t')
		# \n ausschließen
		if(len(line)==0):
			continue
		# writes und andere Strings übernehmen
		elif (len(point)==1):
			points.append(point[0])
		# rest in Arrays umwandeln
		else:
			points.append(np.fromstring(line, dtype=float, count=-1, sep='\t'))
		# Fortschrittsbalken
		bar.update(i)
		i += 1
		
	bar.finish()
	return points
	
# Stellt den Schreibprozess als Video dar.
# TODO: dateinamen an gnuplot weitergeben, evtl über systemvariablen?
def animate(filename='writeprocess'):
	'''
	process0 = subprocess.Popen(
		[
		'mkdir',
		'./animation'],
		stdout=sys.stdout,
		stderr=sys.stderr
		)
	output = process0.communicate()[0]
	'''
	
	print("Plots erstellen...")
	DatatoPNG = "gnuplot animation/animate.gp"
	process1 = subprocess.Popen(DatatoPNG.split(), stdout=subprocess.PIPE)
	output = process1.communicate()[0]
	
	# umbenennen
	#rename frame frame0 frame?????
	
	print("Plots erstellt!")
	
	# chr(39) = '
	# chr(34) = "
	'''
	print('Video erstellen...')
	process2 = subprocess.Popen(
		[
		'mencoder',
		'mf:///home/lukas/Scripts/shapes/animation/*.png',
		'-mf','fps=25',
		'-ovc','lavc',
		'-lavcopts','vcodec=mpeg4',
		'-o',filename + '.avi',
		],
		stdout=sys.stdout,
		stderr=sys.stderr
		)
	output = process2.communicate()[0]
	print('Video erstellt!')
	
	
	print("Aufräumen...")
	process3 = subprocess.Popen(
		[
		'rm', '-R',
		'./animation'],
		stdout=sys.stdout,
		stderr=sys.stderr
		)
	output = process3.communicate()[0]
	print("Aufgeräumt!")
	'''
	
	print('Ich bin animiert!')

#####################
##### Ausführen #####
#####################

#writePoints(spongeLasher(16, 240, 20, 'Hyperuniformstrukturen/270mu_hpu_unscaled_del'), '../../Proben/probe07_2015_06_01/hpu_lashed_240um.gwl', 'w')
#writePoints(spongeLasher(16, 100, 15, 'Hyperuniformstrukturen/hpu_r50_h16_scal0.6'), '/home/lukas/Proben/probe09_2015_06_15/hpu_wall-lashed_100um_scale=0.6.gwl', 'w')
#writePoints(shrinkerOnStakes(19.5,0.25,20), '/home/lukas/Proben/probe08_2015_06_/shrinkerOnStakes_v=15.gwl', 'w')
#writePoints(shrinkerOnStakes(18.75,0.25,19.25), '/home/lukas/Proben/probe09_2015_06_15/shrinkerOnStakes_scal0.6.gwl', 'w')
animate()
#writePoints(scaleStructure(readPoints('./hexagon'), 0.5), './hexagon_scal0.5')

#writePoints(readPoints('./Hyperuniformstrukturen/3x3_hpu_unscaled'))
#p = readPoints('./shapesOut')
#print(ff(p))
