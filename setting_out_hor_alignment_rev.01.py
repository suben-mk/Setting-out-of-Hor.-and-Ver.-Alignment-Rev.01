# Setting Out of Horizontal Alignment Program Rev.01

import numpy as np 

################################## Function List ##################################

# Convert Angle : แปลงมุม (Credit Prajuab Riabroy's Blog)
PI = np.pi
DEG2RAD = PI / 180.0
RAD2DEG= 180.0 / PI

# แปลง degree > deg,min,sec (Credit Prajuab Riabroy's Blog)
def deg2dms(dd):
  sign=1
  if (dd < 0):
    sign = -1
  dd=abs(dd)
  minutes, seconds = divmod(dd*3600,60)
  degrees, minutes = divmod(minutes,60)
  return (sign, degrees, minutes, seconds)

# แปลง deg,min,sec > deg - min - sec (Credit Prajuab Riabroy's Blog)
def DMS2str(degree, minute, second, numdec):
  degree = abs(degree); minute = abs(minute); second = abs(second)
  s ='{:.%df}' % numdec
  ss = s.format(second)
  smin = s.format(60.0)
  mm ='{:.0f}'.format(minute)
  if (ss == smin):
    minute += 1
    ss = s.format(0.0)
    if (minute >= 60):
      mm = '0'
      degree += 1
    else:
      mm ='{:.0f}'.format(minute)
  return '{:.0f}'.format(degree)+"-"+mm+"-"+ss

#Direction Function (Distance and Azimuth)
def Direction(Estart, Nstart, Eend, Nend):
	dE = Eend - Estart
	dN = Nend - Nstart

	Dist = np.sqrt(dE**2 + dN**2)
	A = np.arctan(dE / dN) * RAD2DEG

	if dN < 0:
		Azi = 180 + A
	elif dE < 0:
		Azi = 360 + A
	else:
		Azi = A		
	return Dist, Azi

# Position Calculation 
def Position_Cal(Ecenter, Ncenter, Azimuth, Chainage, Offset):
	Ni = Ncenter + Chainage * np.cos(Azimuth * DEG2RAD) + Offset * np.cos((Azimuth + 90) * DEG2RAD)
	Ei = Ecenter + Chainage * np.sin(Azimuth * DEG2RAD) + Offset * np.sin((Azimuth + 90) * DEG2RAD)
	return Ni, Ei

#---------------------------------- Circular Curve ----------------------------------#
def CirCurve(D, Rc):
	# 1.Length of Curve from PC to PT (m.)
	Lc = Rc * (D * DEG2RAD)

	# 2.Tangent Line from PC to PIc or PT to PIc (m.)
	Tc = Rc * np.tan((D / 2) * DEG2RAD)
	return Lc, Tc

def HorPC(CPI, EPI, NPI, Az1, Tc):
	COOR = Position_Cal(EPI, NPI, Az1, -Tc, 0)
	CPC = CPI - Tc
	EPC = COOR[1]
	NPC = COOR[0]
	AzPC = Az1
	return CPC, EPC, NPC, AzPC

def HorPT(CPC, EPI, NPI, Az2, Tc, Lc):
	COOR = Position_Cal(EPI, NPI, Az2, Tc, 0)
	CPT = CPC + Lc
	EPT = COOR[1]
	NPT = COOR[0]
	AzPT = Az2
	return CPT, EPT, NPT, AzPT

#---------------------------------- Transition Curve ----------------------------------#
def TraCurve(D, Rc, Ls):		
	# 1.Spiral Angle and Circular Angle (rad.)
	Qs = Ls / (2 * Rc)
	Qc = (D * DEG2RAD) - (2 * Qs)

	# 2.Length of Circular (m.)
	Lc = Qc * Rc

	# 3.Offset Xs, Ys
	C1 = 1 / 3 ; C2 = -1 / 10 ; C3 = -1 / 42 ; C4 = 1 / 216 ; C5 = 1 / 1320 ; C6 = -1 / 9360 ; C7 = -1 / 75600 ; C8 = 1 / 685440
	Xs = Ls * (1 + (C2 * Qs ** 2) + (C4 * Qs ** 4) + (C6 * Qs ** 6) + (C8 * Qs ** 8))
	Ys = Ls * ((C1 * Qs) + (C3 * Qs ** 3) + (C5 * Qs ** 5) + (C7 * Qs ** 7))

	# 4.Offset from PCO Tangent to New Curve (m.)
	P = Ys - Rc * (1 - np.cos(Qs))
	  
	# 5.Distance from PCO Tangent to New Curve (m.)
	k = Xs - Rc * np.sin(Qs)
	  
	# 6.Tangent Line from TS to PI or ST to PI (m.)
	Ts = (Rc + P) * np.tan((D / 2) * DEG2RAD) + k
	return Lc, Ts, Xs, Ys

def HorTS(CPI, EPI, NPI, Az1, Ts):
	COOR = Position_Cal(EPI, NPI, Az1, -Ts, 0)
	CTS = CPI - Ts
	ETS = COOR[1]
	NTS = COOR[0]
	AzTS = Az1
	return CTS, ETS, NTS, AzTS  

def HorSC(CTS, ETS, NTS, LS1, Rc, Az1, Az2, Xs, Ys):
	CSC = CTS + LS1

	if Az2 - Az1 < 0:
		COOR = Position_Cal(ETS, NTS, Az1, Xs, -Ys)
		ESC = COOR[1]
		NSC = COOR[0]
	else:
		COOR = Position_Cal(ETS, NTS, Az1, Xs, Ys)
		ESC = COOR[1]
		NSC = COOR[0]

	if Az2 - Az1 < 0:
	  AzSC = Az1 - (LS1 / (2 * Rc)) * RAD2DEG
	else:
	  AzSC = Az1 + (LS1 / (2 * Rc)) * RAD2DEG
	return CSC, ESC, NSC, AzSC

def HorST(CSC, EPI, NPI, Az2, LS2, Lc, Ts):
	COOR = Position_Cal(EPI, NPI, Az2, Ts, 0)
	CST = CSC + Lc + LS2
	EST = COOR[1]
	NST = COOR[0]
	AzST = Az2
	return CST, EST, NST, AzST  

def HorCS(CST, EST, NST, LS2, Rc, Az1, Az2, Xs, Ys):
	CCS = CST - LS2

	if Az2 - Az1 < 0:
		COOR = Position_Cal(EST, NST, Az2, -Xs, -Ys)
		ECS = COOR[1]
		NCS = COOR[0]
	else:
		COOR = Position_Cal(EST, NST, Az2, -Xs, Ys)
		ECS = COOR[1]
		NCS = COOR[0]

	if Az2 - Az1 < 0:
	  AzCS = Az2 + (LS2 / (2 * Rc)) * RAD2DEG
	else:
	  AzCS = Az2 - (LS2 / (2 * Rc)) * RAD2DEG
	return CCS, ECS, NCS, AzCS


################################## Main Program ##################################

# นำเข้าข้อมูล HPI Data.csv (ประเภทข้อความ)                 
HPI = np.loadtxt('01_IMPORT_HSR-HPI.csv', delimiter=',', skiprows=1, dtype=str)

# input point เริ่มต้น [Chainage, Easting, Northing]
BEGIN_POINT = [26787.988, 688703.7479, 1518287.3616]

# นับจำนวณ HPI ใน List
Count = len(HPI)

# หัวตาราง SETTING OUT OF HORIZONTAL ALIGNMENT DATA
print('SETTING OUT OF HORIZONTAL ALIGNMENT DATA')
print('HPI NO., POINT, CHAINAGE, EASTING, NORTHING, AZIMUTH, DEFLECTION ANGLE, RADIUS, LS IN, LS OUT, LENGTH OF CURVE')

# กำหนดชื่อข้อมูลตาม Column ใน HPI Array
for i in range(Count):
	HPI_NO = HPI[i][0]		# Column 1 is HPI NO. (ประเภทข้อความ)
	E = float(HPI[i][1])	# Column 2 is Easting (ประเภทตัวเลข)
	N = float(HPI[i][2])	# Column 3 is Northing (ประเภทตัวเลข)
	R = float(HPI[i][3])	# Column 4 is Radius (ประเภทตัวเลข)
	LS1 = float(HPI[i][4])	# Column 5 is Spital Length in (ประเภทตัวเลข)
	LS2 = float(HPI[i][5])	# Column 6 is Spital Length out (ประเภทตัวเลข)

	# กรณี PI no curve (PI)
		# A คือ HPI Point แรก, B คือ HPI Point ตรงกลาง, C คือ HPI Point สุดท้าย
	if R == 0:

		# Begin Point
		if i == 0:
			DireAB = Direction(E, N, float(HPI[i+1][1]), float(HPI[i+1][2]))
			AzAB = DireAB[1] 	# Azimuth Point.A to Point.B
			CBP = BEGIN_POINT[0]		# CBP คือ Chainage begining point

			# เปลี่ยนค่า Beginning point จะใช้ค่า HPI เป็นจุดเริ่มในการคำนวณโค้งถัดไป
			BEGIN_POINT = [CBP, E, N]

			print('{}, PI, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(HPI_NO, CBP, E, N, AzAB))

		# End point หรือ PI no curve ใดๆ
		else:
			DireBC = Direction(float(HPI[i-1][1]), float(HPI[i-1][2]), E, N)
			AzBC = DireBC[1] 	# Azimuth Point.B to Point.C
			DireBPPI = Direction(BEGIN_POINT[1], BEGIN_POINT[2], E, N)
			DistBPPI = DireBPPI[0] 	# ระยะจาก Begining point to HPI 
			CPI = BEGIN_POINT[0] + DistBPPI 	# CPI คือ Chainage of HPI

			# เปลี่ยนค่า Beginning point จะใช้ค่า HPI เป็นจุดเริ่มในการคำนวณโค้งถัดไป
			BEGIN_POINT = [CPI, E, N]

			print('{}, PI, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(HPI_NO, CPI, E, N, AzBC))

	# กรณี Simple curve (PC, PI, PT)
	  # A คือ HPI Point แรก, B คือ HPI Point ตรงกลาง, C คือ HPI Point สุดท้าย
	elif R != 0 and LS1 == 0:
		DireAB = Direction(float(HPI[i-1][1]), float(HPI[i-1][2]), E, N)
		AzAB = DireAB[1] 	# Azimuth Point.A to Point.B
		DireBPPI = Direction(BEGIN_POINT[1], BEGIN_POINT[2], E, N)
		DistBPPI = DireBPPI[0] 	# ระยะจาก Begining point to HPI 

		DireBC = Direction(E, N, float(HPI[i+1][1]), float(HPI[i+1][2]))
		AzBC = DireBC[1] 	# Azimuth Point.B to Point.C

		# Deflection angle
		D = np.abs(AzBC - AzAB)
		sn, d, m, s = deg2dms(D) # แปลง Deflection angle Deg>DMS
		D_DMS = DMS2str(d, m, s, 2) # (ประเภทข้อความ)

		if AzBC - AzAB < 0: 
			D_String = D_DMS + " " + 'LT'
		else:
		  D_String = D_DMS + " " + 'RT'	

		# Simple Curve Data
		Curve_data = CirCurve(D, R)
		Lc = Curve_data[0]
		Tc = Curve_data[1]

		# Calculation of Simple Curve (HPI, PC, PT)
		CPI = BEGIN_POINT[0] + DistBPPI
		PC = HorPC(CPI, E, N, AzAB, Tc)
		PT = HorPT(PC[0], E, N, AzBC, Tc, Lc)

		# เปลี่ยนค่า Beginning point จะใช้ค่า PT เป็นจุดเริ่มในการคำนวณโค้งถัดไป
		BEGIN_POINT = [PT[0], PT[1], PT[2]]

		print('-, PC, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(PC[0], PC[1], PC[2], PC[3]))
		print('{}, PI, {:.3f}, {:.5f}, {:.5f}, -, {}, {:.3f}, -, -, {:.3f}'.format(HPI_NO, CPI, E, N, D_String, R, Lc))
		print('-, PT, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(PT[0], PT[1], PT[2], PT[3]))


	# กรณี Trasition curve (TS, SC, PI, CS, ST)
	  # A คือ HPI Point แรก, B คือ HPI Point ตรงกลาง, C คือ HPI Point สุดท้าย
	else:
		DireAB = Direction(float(HPI[i-1][1]), float(HPI[i-1][2]), E, N)
		AzAB = DireAB[1] 	# Azimuth Point.A to Point.B
		DireBPPI = Direction(BEGIN_POINT[1], BEGIN_POINT[2], E, N)
		DistBPPI = DireBPPI[0] 	# ระยะจาก Begining point to HPI 

		DireBC = Direction(E, N, float(HPI[i+1][1]), float(HPI[i+1][2]))
		AzBC = DireBC[1] 	# Azimuth Point.B to Point.C

		# Deflection angle
		D = np.abs(AzBC - AzAB)
		sn, d, m, s = deg2dms(D) # แปลง Deflection angle Deg>DMS
		D_DMS = DMS2str(d, m, s, 2) # (ประเภทข้อความ)

		if AzBC - AzAB < 0: 
			D_String = D_DMS + " " + 'LT'
		else:
		  D_String = D_DMS + " " + 'RT'	

		# Transition Curve Data
		Curve_data = TraCurve(D, R, LS1)
		Lc = Curve_data[0]
		Ts = Curve_data[1]
		Xs = Curve_data[2]
		Ys = Curve_data[3]

		# Calculation of Transition Curve (HPI, TS, SC, CS, ST)
		CPI = BEGIN_POINT[0] + DistBPPI
		TS = HorTS(CPI, E, N, AzAB, Ts)
		SC = HorSC(TS[0], TS[1], TS[2], LS1, R, AzAB, AzBC, Xs, Ys)
		ST = HorST(SC[0], E, N, AzBC, LS2, Lc, Ts)
		CS = HorCS(ST[0], ST[1], ST[2], LS2, R, AzAB, AzBC, Xs, Ys)

		# เปลี่ยนค่า Beginning point จะใช้ค่า ST เป็นจุดเริ่มในการคำนวณโค้งถัดไป
		BEGIN_POINT = [ST[0], ST[1], ST[2]]

		print('-, TS, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(TS[0], TS[1], TS[2], TS[3]))
		print('-, SC, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(SC[0], SC[1], SC[2], SC[3]))
		print('{}, PI, {:.3f}, {:.5f}, {:.5f}, -, {}, {:.3f}, {:.3f}, {:.3f}, {:.3f}'.format(HPI_NO, CPI, E, N, D_String, R, LS1, LS2, Lc))
		print('-, CS, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(CS[0], CS[1], CS[2], CS[3]))
		print('-, ST, {:.3f}, {:.5f}, {:.5f}, {:.5f}'.format(ST[0], ST[1], ST[2], ST[3]))