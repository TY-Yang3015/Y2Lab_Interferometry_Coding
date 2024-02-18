from Simulation import InterferogramSimulator

x, y = InterferogramSimulator.auto_mesh()
simulator = InterferogramSimulator(x, y, 50, 5)

simulator.add_square(550e-9,10e-9, 1.0)
simulator.add_square(560e-9, 10e-9, 1.0)
#simulator.add_gaussian(550e-9,0.1e-9, 1.0)
#simulator.add_gaussian(560e-9, 0.1e-9, 1.0)
simulator.inspect_simulation()
simulator.visualise(interactive=True, save=False)
simulator.show_spectrum(interactive=False, save=False, manual_lim=None);
