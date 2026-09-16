# -------------------------------------------------------------
# Left electrode
# -------------------------------------------------------------

# Set up lattice
vector_a = [5.4306, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.4306, 0.0]*Angstrom
vector_c = [0.0, 0.0, 5.4306]*Angstrom
left_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
left_electrode_elements = [Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                           Silicon]

# Define coordinates
left_electrode_coordinates = [[ 0.      ,  0.      ,  0.678825],
                              [ 2.7153  ,  2.7153  ,  0.678825],
                              [ 1.35765 ,  1.35765 ,  2.036475],
                              [ 4.07295 ,  4.07295 ,  2.036475],
                              [ 2.7153  ,  0.      ,  3.394125],
                              [ 0.      ,  2.7153  ,  3.394125],
                              [ 4.07295 ,  1.35765 ,  4.751775],
                              [ 1.35765 ,  4.07295 ,  4.751775]]*Angstrom

# Set up configuration
left_electrode = BulkConfiguration(
    bravais_lattice=left_electrode_lattice,
    elements=left_electrode_elements,
    cartesian_coordinates=left_electrode_coordinates
    )

# -------------------------------------------------------------
# Right electrode
# -------------------------------------------------------------

# Set up lattice
vector_a = [5.4306, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.4306, 0.0]*Angstrom
vector_c = [0.0, 0.0, 5.4306]*Angstrom
right_electrode_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
right_electrode_elements = [Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                            Silicon]

# Define coordinates
right_electrode_coordinates = [[ 0.      ,  0.      ,  0.678825],
                               [ 2.7153  ,  2.7153  ,  0.678825],
                               [ 1.35765 ,  1.35765 ,  2.036475],
                               [ 4.07295 ,  4.07295 ,  2.036475],
                               [ 2.7153  ,  0.      ,  3.394125],
                               [ 0.      ,  2.7153  ,  3.394125],
                               [ 4.07295 ,  1.35765 ,  4.751775],
                               [ 1.35765 ,  4.07295 ,  4.751775]]*Angstrom

# Set up configuration
right_electrode = BulkConfiguration(
    bravais_lattice=right_electrode_lattice,
    elements=right_electrode_elements,
    cartesian_coordinates=right_electrode_coordinates
    )

# -------------------------------------------------------------
# Central region
# -------------------------------------------------------------

# Set up lattice
vector_a = [5.4306, 0.0, 0.0]*Angstrom
vector_b = [0.0, 5.4306, 0.0]*Angstrom
vector_c = [0.0, 0.0, 21.7224]*Angstrom
central_region_lattice = UnitCell(vector_a, vector_b, vector_c)

# Define elements
central_region_elements = [Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                           Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                           Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                           Silicon, Silicon, Silicon, Silicon, Silicon, Silicon, Silicon,
                           Silicon, Silicon, Silicon, Silicon]

# Define coordinates
central_region_coordinates = [[  0.      ,   0.      ,   0.678825],
                              [  2.7153  ,   2.7153  ,   0.678825],
                              [  1.35765 ,   1.35765 ,   2.036475],
                              [  4.07295 ,   4.07295 ,   2.036475],
                              [  2.7153  ,   0.      ,   3.394125],
                              [  0.      ,   2.7153  ,   3.394125],
                              [  4.07295 ,   1.35765 ,   4.751775],
                              [  1.35765 ,   4.07295 ,   4.751775],
                              [  0.      ,   0.      ,   6.109425],
                              [  2.7153  ,   2.7153  ,   6.109425],
                              [  1.35765 ,   1.35765 ,   7.467075],
                              [  4.07295 ,   4.07295 ,   7.467075],
                              [  2.7153  ,   0.      ,   8.824725],
                              [  0.      ,   2.7153  ,   8.824725],
                              [  4.07295 ,   1.35765 ,  10.182375],
                              [  1.35765 ,   4.07295 ,  10.182375],
                              [  0.      ,   0.      ,  11.540025],
                              [  2.7153  ,   2.7153  ,  11.540025],
                              [  1.35765 ,   1.35765 ,  12.897675],
                              [  4.07295 ,   4.07295 ,  12.897675],
                              [  2.7153  ,   0.      ,  14.255325],
                              [  0.      ,   2.7153  ,  14.255325],
                              [  4.07295 ,   1.35765 ,  15.612975],
                              [  1.35765 ,   4.07295 ,  15.612975],
                              [  0.      ,   0.      ,  16.970625],
                              [  2.7153  ,   2.7153  ,  16.970625],
                              [  1.35765 ,   1.35765 ,  18.328275],
                              [  4.07295 ,   4.07295 ,  18.328275],
                              [  2.7153  ,   0.      ,  19.685925],
                              [  0.      ,   2.7153  ,  19.685925],
                              [  4.07295 ,   1.35765 ,  21.043575],
                              [  1.35765 ,   4.07295 ,  21.043575]]*Angstrom

# Set up configuration
central_region = BulkConfiguration(
    bravais_lattice=central_region_lattice,
    elements=central_region_elements,
    cartesian_coordinates=central_region_coordinates
    )

device_configuration = DeviceConfiguration(
    central_region,
    [left_electrode, right_electrode]
    )

# -------------------------------------------------------------
# Calculator
# -------------------------------------------------------------

potentialSet = StillingerWeber_Si_1985()
calculator = TremoloXCalculator(parameters=potentialSet)

device_configuration.setCalculator(calculator)
device_configuration.update()

# -------------------------------------------------------------
# Phonon transmission spectrum
# -------------------------------------------------------------
# --- Y-2026.03 FIX ---------------------------------------------------------
# The SHIPPED vendor example is STALE: it calls
#     PhononTransmissionSpectrum(configuration=device_configuration, ...)
# which fails on Y-2026.03 with
#     TypeError: __init__() missing 1 required positional argument: 'dynamical_matrix'
# Manual signature (verified):
#   PhononTransmissionSpectrum(dynamical_matrix, configuration=None, energies=None,
#       qpoints=None, qpoints_weights=None, self_energy_calculator=None, infinitesimal=None)
# So a DynamicalMatrix must be built FIRST and passed in.
print("PTS_START")

dynamical_matrix = DynamicalMatrix(
    configuration=device_configuration,
    filename='si_tbc.hdf5',
    object_id='dm',
    calculator=calculator,
    repetitions=(3, 3, 1),        # transverse repeats; transport is along C
    # CRITICAL (manual-verified): use_equivalent_bulk DEFAULTS TO True, which treats a
    # DeviceConfiguration as a BulkConfiguration -> no electrode structure -> the device
    # transmission routine then hits 'NoneType has no attribute makeFull'.
    # Must be False for device phonon transmission.
    use_equivalent_bulk=False,
    constrain_electrodes=False,   # manual default; electrodes must be computed
)
dynamical_matrix.update()
print("DM_OK")

phonon_transmission_spectrum = PhononTransmissionSpectrum(
    dynamical_matrix=dynamical_matrix,
    energies=numpy.linspace(0, 0.5, 101)*eV,
    qpoints=MonkhorstPackGrid(5, 5),
    infinitesimal=1e-06*eV,
    self_energy_calculator=RecursionSelfEnergy(),
)
nlsave('si_tbc.hdf5', phonon_transmission_spectrum)
print("PTS_OK")

# --- report: transmission and the Landauer thermal conductance ---------------
T = numpy.array(phonon_transmission_spectrum.evaluate())
E = numpy.array(phonon_transmission_spectrum.energies().inUnitsOf(eV))
print("TRANS_MAX %.4f" % T.max())
print("TRANS_AT_LOW_E %.4f" % T[1])

# Landauer: G = (1/2pi) * int dE * E * T(E) * dn/dT   (per unit area handled later)
kB = 8.617333e-5      # eV/K
Temp = 300.0
x = E / (kB * Temp)
x = numpy.where(x < 1e-9, 1e-9, x)
ex = numpy.exp(numpy.clip(x, 0, 700))
dndT = (x**2) * ex / ((ex - 1.0)**2) / Temp        # dimensionless weight / K
integrand = T * dndT
G = numpy.trapz(integrand, E)                      # eV/K per mode-normalisation
print("LANDAUER_G_rel %.6e" % G)
print("PHASE1D_VENDOR_OK")
