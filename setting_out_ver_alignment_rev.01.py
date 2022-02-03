# Setting Out of Vertical Alignment Program Rev.01

import numpy as np 

################################## Vertical Curve Function List ##################################

def Gradient(Ch1, El1, Ch2, El2):
	Gradient = ((El2 - El1) / (Ch2 - Ch1)) * 100
	return Gradient

def VerBVC(CPI, EPI, g1, L):
	Chainage = CPI - L / 2
	Elev = EPI - (g1 / 100) * (L / 2)
	return Chainage, Elev

def VerEVC(CPI, EPI, g2, L):
	Chainage = CPI + L / 2
	Elev = EPI + (g2 / 100) * (L / 2)
	return Chainage, Elev

def RCurve(VLC, g1, g2):
	R = VLC * 100 / (g2 - g1)
	return R

################################## Main Program ##################################

# นำเข้าข้อมูล PVI.csv (ประเภทข้อความ)                 
PVI = np.loadtxt('02_IMPORT_HSR-PVI.csv', delimiter=',', skiprows=1, dtype=str)

# นับจำนวณ PVI ใน List
Count = len(PVI)

# หัวตาราง SETTING OUT OF VERTICAL ALIGNMENT DATA
print('SETTING OUT OF VERTICAL ALIGNMENT DATA')
print('PVI NO., POINT, CHAINAGE, ELEVATION, GRADIENT IN, GRADIENT OUT, LENGTH OF CURVE, RADIUS')

# กำหนดชื่อข้อมูลตาม Column ใน PVI Array
for i in range(Count):
	PVI_NO = PVI[i][0]		# Column 1 is PVI NO. (ประเภทข้อความ)
	CH = float(PVI[i][1])	# Column 2 is Chainage (ประเภทตัวเลข)
	EL = float(PVI[i][2])	# Column 3 is Elevation (ประเภทตัวเลข)
	LVC = float(PVI[i][3])	# Column 4 is Length of Curve (ประเภทตัวเลข)

	# กรณี PVI no curve
	if LVC == 0:
		# Beginning Point
		if i == 0:
			g_out = Gradient(CH, EL, float(PVI[i+1][1]), float(PVI[i+1][2]))
			g_in = g_out
			print('{}, BP, {:.3f}, {:.5f}, {:.3f}, {:.3f}, -, -'.format(PVI_NO, CH, EL, g_in, g_out))

		# End point
		else:
			g_in = Gradient(float(PVI[i-1][1]), float(PVI[i-1][2]), CH, EL)
			g_out = g_in
			print('{}, EP, {:.3f}, {:.5f}, {:.3f}, {:.3f}, -, -'.format(PVI_NO, CH, EL, g_in, g_out))

	# กรณี Vertical Curve แบบ Parabola
	else:
		# Beginning vertical curve
		BVC_g_in = Gradient(float(PVI[i-1][1]), float(PVI[i-1][2]), CH, EL)
		BVC_g_out = Gradient(CH, EL, float(PVI[i+1][1]), float(PVI[i+1][2]))
		BVC = VerBVC(CH, EL, BVC_g_in, LVC)
		BVC_CH = BVC[0]
		BVC_EL = BVC[1]

		# Ending vertical curve
		EVC_g_in = Gradient(CH, EL, float(PVI[i+1][1]), float(PVI[i+1][2]))
		EVC_g_out = Gradient(CH, EL, float(PVI[i+1][1]), float(PVI[i+1][2]))
		EVC = VerEVC(CH, EL, EVC_g_out, LVC)
		EVC_CH = EVC[0]
		EVC_EL = EVC[1]

		R = RCurve(LVC, BVC_g_in, EVC_g_out)		

		print('-, BVC, {:.3f}, {:.5f}, {:.3f}, {:.3f}, -, -'.format(BVC_CH, BVC_EL, BVC_g_in, BVC_g_out))
		print('{}, PVI, {:.3f}, {:.5f}, -, -, {:.3f}, {:.3f}'.format(PVI_NO, CH, EL, LVC, R))
		print('-, EVC, {:.3f}, {:.5f}, {:.3f}, {:.3f}, -, -'.format(EVC_CH, EVC_EL, EVC_g_in, EVC_g_out))