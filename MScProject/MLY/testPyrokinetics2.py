from pyrokinetics import Pyro

import matplotlib.pyplot as plt

inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc8/gxLinCbc8.in"
# inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/stepEcHd_MG2/input.in"
# inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/stepEcHd_MG2/input.used_inputs.in"
# outputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/stepEcHd_MG2/input.out.nc"
# inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/gxLinCbc8_Viking/input.in"
# inputFilePath = "/home/mly509/GithubCloneMly/Nuclear-Fusion/MScProject/MLY/stepEcHd_1/input.in"

pyro = Pyro(gk_file=inputFilePath, gk_code='GS2')

pyro.load_gk_output()

data = pyro.gk_output

gamma = data['growth_rate']

gamma.isel(kx=0).plot.line(x='time')

plt.show()