from main import *
NPRO_wavelength = 1064.50*nm

y = ABCD_Matrix_of.fabry_perot_cavity(57*m, 38.79*m)
fbc = FabryPerotCavity(37*m, 57*m, 1, 1, "FBC", "yellow")
fbc.stable_waist(NPRO_wavelength/2)


devices = [BeamPathElement(16 * inch, 3 * inch, "Start", w0=217*um, z0=5*cm, wavelength=NPRO_wavelength,
                           color="red"),
           # half-wave plate (26, 3)
           Mirror(28.5 * inch, 3 * inch, "\nM1", color="red"),
           ThinLens(75 * mm, 28.5 * inch, 8 * inch, "\nL1", color="red"),
           Mirror(28.5 * inch, 9 * inch, "\nM2", color="red"),
           FaradayIsolator(1, 2.2321, 2*cm, 25.9 * inch, 9 * inch, "\nFI1", "red"),
           # half-wave plate (21, 9)
           ThinLens(250 * mm, 17 * inch, 9 * inch, "\nL2", color="red"),
           ThinLens(250 * mm, 9 * inch, 8 * inch, "\nL3", color="red"),
           Mirror(4.5 * inch, 9 * inch, "\nM3", color="red"),
           Mirror(4.5 * inch, 12 * inch, "\nM4", color="red"),
           SHG(1, 1.7379, 1.7779, 1.5*cm, 1.5*cm, 11*inch, 12*inch, "\nSHG", color="green"),
           # DSP (14.5, 12) in 1064nm
           ThinLens(250 * mm, 17.2 * inch, 12.4 * inch, "\nL4", color="green"),
           # half-wave plate (21, 12.5)
           Mirror(25.5 * inch, 12.9 * inch, "\nM5", color="green"),
           ThinLens(200 * mm, 25.5 * inch, 15 * inch, "\nL5", color="green"),
           FaradayIsolator(1, 2.3232, 2 * cm, 26 * inch, 24 * inch, "\nFI2", "green"),
           ThinLens(-125 * mm, 26 * inch, 29 * inch, "\nL6", color="green"),
           # half-wave plate (26, 39)
           Mirror(26 * inch, 41 * inch, "\nM6", color="green"),
           Mirror(21 * inch, 41 * inch, "\nM7", color="green"),
           Mirror(21 * inch, 35 * inch, "\nM8", color="green"),
           ThinLens(125 * mm, 28.5 * inch, 34 * inch, "\nL7", color="green"),
           Mirror(30.5 * inch, 34 * inch, "\nM9", color="green"),
           # open shutter (30, 29)
           # open shutter? (30, 25)
           Mirror(30 * inch, 21.5 * inch, "\nM10", color="green"),
           # half-wave plate, open? (18, 21.5)
           # DSP (14.5, 21.75)
           ThinLens(1000 * mm, 11 * inch, 21.75 * inch, "\nL8", color="green"),
           # half-wave plate, open? (6, 21.5)
           BeamPathElement(117 * cm, 0, "\nETMY")]

beam_path = LaserBeamSetup(devices)
beam_path.plot_beam_path()

beam_path = LaserBeamSetup(devices)
beam_path.plot_beam_path()

