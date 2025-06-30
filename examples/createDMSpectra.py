# code to demonstrate how to create a spectra object for dark matter models
# author: Andrea Albert (aalbert@slac.stanford.edu)
# date: Oct 26, 2016

from threeML import *

# DMFitFunction uses the Pythia-generated table from the standard Fermi Science Tools which is appropriate for 2 GeV < mass < 10 TeV
# DMSpectra combines the Pythia-generated table from the standard Fermi Science Tools (for 2 GeV < mass < 10 TeV) and the one used by the HAWC Collaboration (for 10 TeV < mass < 1 PeV).

spec = DMFitFunction()
spec2 = DMSpectra()

channels = spec.print_channel_mapping()

spec.mass = 1000.0  # mass needs to be in GeV
spec.channel = 4  # bb
spec.sigmav = 3.0e-26  # sigmav must be in cm^3 / s
spec.J = 1.0e20  # J-factor must be in GeV^2 / cm^5

spec2.mass = 50000.0  # mass needs to be in GeV
spec2.channel = 4  # bb
spec2.sigmav = 3.0e-26  # sigmav must be in cm^3 / s
spec2.J = 1.0e20  # J-factor must be in GeV^2 / cm^5

en = np.logspace(6.0, 11.0, 100)  # en is in keV.  so from 1 GeV to 100 TeV

DMspec = spec.evaluate(
    en, spec.mass.value, spec.channel.value, spec.sigmav.value, spec.J.value
)
DMspec2 = spec.evaluate(
    en, spec2.mass.value, spec2.channel.value, spec2.sigmav.value, spec2.J.value
)
