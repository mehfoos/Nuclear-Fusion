import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from netCDF4 import Dataset
import itertools

from t3d import Profiles, Grid


class plotter:

    def __init__(self):
        # MLY changes
        self.N_cumulative_last = 0 # The acumulative "N" up to this file (not including this file or future files). Set in calling function.
        self.time_cumulative_last = 0 # Similar, but for time
        self.global_first = False # flag if first sim of a set
        self.global_last = False # flag if last sim of a set
        # for backwards compat to before t3d packaging
        sys.modules['Profiles'] = Profiles
        sys.modules['Grid'] = Grid

        np.set_printoptions(linewidth=500)

        self.grid_lines = False

        # Dictionary of plotting functions for the legacy panel plot
        self.panel_plots = {
            "density": self.plot_density_profiles,
            "temperature": self.plot_temperature_profiles,
            "pressure": self.plot_pressure_profiles,
            "power_balance": self.plot_power_balance,
            "heat_flux": self.plot_heat_flux,
            "particle_balance": self.plot_particle_balance,
            "particle_flux": self.plot_particle_flux,
            #"fusion_power_density": self.plot_fusion_power_density,
            #"radiation": self.plot_radiation_profiles,
            "collisional_equilibration": self.plot_collisional_equilibration,
            "power_source": self.plot_auxiliary_power_source,
            "q": self.plot_q,
            "particle_source": self.plot_auxiliary_particle_source,
            "gamma": self.plot_gamma,
            #"state_profiles": self.plot_state_profiles,
            #"flux": self.plot_flux,
            #"source_total": self.plot_source_total,
            #"sink_terms": self.plot_sink_terms,
        }

        # Dictionary of all available plotting functions
        self.plotting_functions = self.panel_plots.copy()
        self.plotting_functions["state_profiles"] = self.plot_state_profiles
        self.plotting_functions["flux"] = self.plot_flux
        self.plotting_functions["source_total"] = self.plot_source_total
        self.plotting_functions["sink_terms"] = self.plot_sink_terms

    def list_available_plots(self):
        ''' Function to list available plots '''
        print("\nList of available plots:")
        for plot in self.plotting_functions:
            print(f"  {plot}")

    def read_data(self, fin):
        ''' Function to read the Trinity3D npy file '''

        # Load the pickle file
        data = np.load(fin, allow_pickle=True).tolist()

        self.t_ref = data['norms']['t_ref']
        try:
            self.rho_star = data['norms']['rho_star']
        except:
            self.rho_star = 1 / data['norms']['gyro_scale']
        self.t_scale = self.t_ref * self.rho_star**-2
        self.time = np.array(data['t']) * self.t_scale
        self.tags = data['species_tags']
        self.ion_tag = data['bulk_ion_tag']

        self.N = len(self.time)

        # set up color maps
        self.warm_map = pylab.cm.autumn(np.linspace(1, 0.25, self.N))
        self.cool_map = pylab.cm.Blues(np.linspace(0.25, 1, self.N))
        self.green_map = pylab.cm.YlGn(np.linspace(0.25, 1, self.N))
        # purple_map = pylab.cm.Purples(np.linspace(0.25, 1, self.N))

        # MLY adding older version stuff Placing here so that it's overriden if there are changes
        self.beta = np.asarray(data['beta_ref'])
        self.Qi_turb = np.asarray(data['qflux_turb_' + self.ion_tag])
        self.Qe_turb = np.asarray(data['qflux_turb_e'])
        self.Q_turb_MW_i = np.asarray(data['Q_turb_MW_' + self.ion_tag])
        self.Q_turb_MW_e = np.asarray(data['Q_turb_MW_e'])
        self.Sp_tot_int_MW_i = np.asarray(data['Sp_tot_int_MW_' + self.ion_tag])
        self.Sp_tot_int_MW_e = np.asarray(data['Sp_tot_int_MW_e'])
        self.Q_MW_i = np.asarray(data['Q_MW_' + self.ion_tag])
        self.Q_neo_MW_i = np.asarray(data['Q_neo_MW_' + self.ion_tag])
        self.Q_MW_e = np.asarray(data['Q_MW_e'])
        self.Q_neo_MW_e = np.asarray(data['Q_neo_MW_e'])


        # Read the species profile data
        self.n = np.asarray(data['n_e'])
        self.ni = np.asarray(data['n_i'])
        self.pi = np.asarray(data['p_' + self.ion_tag])
        self.pe = np.asarray(data['p_e'])
        self.Ti = np.asarray(data['T_' + self.ion_tag])
        self.Te = np.asarray(data['T_e'])
        self.Gamma = {}
        self.Qi = {}
        self.Qe = {}
        self.Gamma['tot'] = np.asarray(data['pflux_' + self.ion_tag])
        self.Qi['tot'] = np.asarray(data['qflux_' + self.ion_tag])
        self.Qe['tot'] = np.asarray(data['qflux_e'])
        self.flux_model_labels = data['flux_model_labels']
        for label in self.flux_model_labels:
            self.Gamma[label] = np.asarray(data['pflux_' + label + '_' + self.ion_tag])
            self.Qi[label] = np.asarray(data['qflux_' + label + '_' + self.ion_tag])
            self.Qe[label] = np.asarray(data['qflux_' + label + '_e'])
        self.aLn = np.asarray(data['aLn_' + self.ion_tag])
        self.aLpi = np.asarray(data['aLp_' + self.ion_tag])
        self.aLpe = np.asarray(data['aLp_e'])
        self.aLTi = np.asarray(data['aLT_' + self.ion_tag])
        self.aLTe = np.asarray(data['aLT_e'])

        

        # Convert pressure sources from code units to MW/m^3 by dividing by p_source_scale
        self.P_ref_MWm3 = data['norms']['P_ref_MWm3']
        self.Sn_ref_SI20 = data['norms']['Sn_ref_SI20']
        try:
            self.a_minor = data['norms']['a_minor']
        except:
            self.a_minor = data['norms']['a_ref']

        # Density sources, in code units
        self.Sn_aux_i = np.asarray(data['Sn_aux_' + self.ion_tag])
        self.Sn_aux_e = np.asarray(data['Sn_aux_e'])

        # Pressure sources, in code units, multiply by P_ref_MWm3 to get MW/m^3
        self.Sp_aux_i = np.asarray(data['Sp_aux_' + self.ion_tag])
        self.Sp_aux_e = np.asarray(data['Sp_aux_e'])
        self.Sp_alpha_i = np.asarray(data['Sp_alpha_' + self.ion_tag])
        self.Sp_alpha_e = np.asarray(data['Sp_alpha_e'])
        self.Sp_rad_i = np.asarray(data['Sp_rad_' + self.ion_tag])
        self.Sp_rad_e = np.asarray(data['Sp_rad_e'])
        if len(self.Sp_rad_i) == 0:  # old files use Sp_brem
            self.Sp_rad_i = np.asarray(data['Sp_brem_' + self.ion_tag])
            self.Sp_rad_e = np.asarray(data['Sp_brem_e'])
        self.Sp_heating_i = np.asarray(data['Sp_heating_' + self.ion_tag])
        self.Sp_heating_e = np.asarray(data['Sp_heating_e'])
        self.Sp_coll_i = np.asarray(data['Sp_coll_' + self.ion_tag])
        self.Sp_coll_e = np.asarray(data['Sp_coll_e'])
        self.Sp_tot_i = np.asarray(data['Sp_tot_' + self.ion_tag])
        self.Sp_tot_e = np.asarray(data['Sp_tot_e'])

        self.Sp_aux_int_MW_i = np.asarray(data['Sp_aux_int_MW_' + self.ion_tag])
        self.Sp_aux_int_MW_e = np.asarray(data['Sp_aux_int_MW_e'])
        self.Sp_aux_int_MW_tot = self.Sp_aux_int_MW_i + self.Sp_aux_int_MW_e
        self.Sp_alpha_int_MW_i = np.asarray(data['Sp_alpha_int_MW_' + self.ion_tag])
        self.Sp_alpha_int_MW_e = np.asarray(data['Sp_alpha_int_MW_e'])
        self.Sp_rad_int_MW_i = np.asarray(data['Sp_rad_int_MW_' + self.ion_tag])
        self.Sp_rad_int_MW_e = np.asarray(data['Sp_rad_int_MW_e'])
        self.Prad_tot = np.asarray(data['Sp_rad_int_MW_e'])
        self.Sp_heating_int_MW_i = np.asarray(data['Sp_heating_int_MW_' + self.ion_tag])
        self.Sp_heating_int_MW_e = np.asarray(data['Sp_heating_int_MW_e'])
        self.Sp_coll_int_MW_i = np.asarray(data['Sp_coll_int_MW_' + self.ion_tag])
        self.Sp_coll_int_MW_e = np.asarray(data['Sp_coll_int_MW_e'])

        self.Sn_tot_cumint_SI20_i = np.asarray(data['Sn_tot_cumint_SI20_' + self.ion_tag])
        self.Sn_tot_cumint_SI20_e = np.asarray(data['Sn_tot_cumint_SI20_e'])
        if len(self.Sn_tot_cumint_SI20_i) == 0:  # backwards compat
            self.Sn_tot_cumint_SI20_i = np.asarray(data['Sn_tot_int_SI20_' + self.ion_tag])
            self.Sn_tot_cumint_SI20_e = np.asarray(data['Sn_tot_int_SI20_e'])

        self.Sp_aux_cumint_MW_i = np.asarray(data['Sp_aux_cumint_MW_' + self.ion_tag])
        self.Sp_aux_cumint_MW_e = np.asarray(data['Sp_aux_cumint_MW_e'])
        self.Sp_alpha_cumint_MW_i = np.asarray(data['Sp_alpha_cumint_MW_' + self.ion_tag])
        self.Sp_alpha_cumint_MW_e = np.asarray(data['Sp_alpha_cumint_MW_e'])
        self.Sp_rad_cumint_MW_i = np.asarray(data['Sp_rad_cumint_MW_' + self.ion_tag])
        self.Sp_rad_cumint_MW_e = np.asarray(data['Sp_rad_cumint_MW_e'])
        self.Sp_heating_cumint_MW_i = np.asarray(data['Sp_heating_cumint_MW_' + self.ion_tag])
        self.Sp_heating_cumint_MW_e = np.asarray(data['Sp_heating_cumint_MW_e'])
        self.Sp_coll_cumint_MW_i = np.asarray(data['Sp_coll_cumint_MW_' + self.ion_tag])
        self.Sp_coll_cumint_MW_e = np.asarray(data['Sp_coll_cumint_MW_e'])
        self.Sp_tot_cumint_MW_i = np.asarray(data['Sp_tot_cumint_MW_' + self.ion_tag])
        self.Sp_tot_cumint_MW_e = np.asarray(data['Sp_tot_cumint_MW_e'])
        if len(self.Sp_tot_cumint_MW_i) == 0:  # backwards compat
            self.Sp_tot_cumint_MW_i = np.asarray(data['Sp_tot_int_MW_' + self.ion_tag])
            self.Sp_tot_cumint_MW_e = np.asarray(data['Sp_tot_int_MW_e'])

        self.Palpha_MWm3 = np.asarray(data['Palpha_MWm3'])
        self.Palpha_int_MW = data['Palpha_int_MW']

        # power balance
        self.Q_MW_i = {}
        self.Q_MW_e = {}
        self.Q_MW_i['tot'] = np.asarray(data['Q_MW_' + self.ion_tag])
        self.Q_MW_e['tot'] = np.asarray(data['Q_MW_e'])
        for label in self.flux_model_labels:
            self.Q_MW_i[label] = np.asarray(data['Q_' + label + '_MW_' + self.ion_tag])
            self.Q_MW_e[label] = np.asarray(data['Q_' + label + '_MW_e'])

        # particle balance
        self.Gam_SI20_i = {}
        self.Gam_SI20_e = {}
        self.Gam_SI20_i['tot'] = np.asarray(data['Gam_SI20_' + self.ion_tag])
        self.Gam_SI20_e['tot'] = np.asarray(data['Gam_SI20_e'])
        for label in self.flux_model_labels:
            self.Gam_SI20_i[label] = np.asarray(data['Gam_' + label + '_SI20_' + self.ion_tag])
            self.Gam_SI20_e[label] = np.asarray(data['Gam_' + label + '_SI20_e'])

        self.grid = data['grid']
        self.N_rho = self.grid['N_radial']
        self.axis = self.grid['rho_axis']
        self.mid_axis = self.grid['mid_axis']
        try:
            self.flux_label = self.grid['flux_label']
        except LookupError:
            self.flux_label = 'rminor'
        if self.flux_label == 'rminor':
            self.rad_label = r'$r/a$'
        elif self.flux_label == 'torflux':
            self.rad_label = r'$\rho_{tor}$'

        self.filename = data['trinity_infile']
        self.data = data


    def plot_density_profiles(self, axs):
        ''' Function to plot the density profiles '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            # plot profiles
            if (i == 0 and self.global_first) or (i == self.N - 1 and self.global_last):
                axs.plot(self.axis, self.n[i], '.-', color=self.green_map[i+self.N_cumulative_last],
                         label=f'$n_e$, $t =${t_curr+self.time_cumulative_last:.2f} s')
            else:
                axs.plot(self.axis, self.n[i], '.-', color=self.green_map[i+self.N_cumulative_last])
        nmax = np.max(self.n)
        axs.set_ylim(bottom=0, top=1.5 * nmax)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'$\bf a)$ density [10$^{20}$ m$^{-3}$]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_temperature_profiles(self, axs):
        ''' Function to plot the temperature profiles '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            # plot profiles
            if (i == 0 and self.global_first):
                axs.plot(self.axis, self.Ti[i], '.-', color=self.warm_map[i+self.N_cumulative_last], label=f'$T_i$, $t =${t_curr:.2f} s')
                axs.plot(self.axis, self.Te[i], '.-', color=self.cool_map[i+self.N_cumulative_last], label=f'$T_e$, $t =${t_curr:.2f} s')
            elif (i == self.N - 1 and self.global_last):
                axs.plot(self.axis, self.Ti[i], '.-', color=self.warm_map[i+self.N_cumulative_last], label=f'$T_i$, $t =${t_curr:.2f} s')
                if self.plot_electrons:
                    axs.plot(self.axis, self.Te[i], '.-', color=self.cool_map[i+self.N_cumulative_last], label=f'$T_e$, $t =${t_curr:.2f} s')
            else:
                axs.plot(self.axis, self.Ti[i], '.-', color=self.warm_map[i+self.N_cumulative_last])
                if self.plot_electrons:
                    axs.plot(self.axis, self.Te[i], '.-', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_ylim(bottom=0)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'$\bf b)$ temperature [keV]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_pressure_profiles(self, axs):
        ''' Function to plot the pressure profiles '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            # plot profiles
            axs.plot(self.axis, self.pi[i], '.-', color=self.warm_map[i+self.N_cumulative_last])
            if self.plot_electrons:
                axs.plot(self.axis, self.pe[i], '.-', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_ylim(bottom=0)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'$\bf c)$ pressure [10$^{20}$m$^{-3}$ keV]')
        axs.grid(self.grid_lines)

    def plot_power_balance(self, axs):
        ''' Function to plot the power balance '''
        if self.global_last:
            i = self.N - 1
            axs.plot(self.mid_axis, self.Q_MW_i['tot'][i], 'o-', color=self.warm_map[i+self.N_cumulative_last], fillstyle='none', label='$Q_i$ (total)')
            marker = itertools.cycle(('d', 's', 'v', 'P', '*', 'X'))
            for label in self.flux_model_labels:
                axs.plot(self.mid_axis, self.Q_MW_i[label][i], '-', marker=next(marker), color=self.warm_map[i+self.N_cumulative_last], fillstyle='none', label=f'$Q_i$ ({label})')
            axs.plot(self.mid_axis, self.Sp_tot_cumint_MW_i[i][:-1], 'o:', color=self.warm_map[i+self.N_cumulative_last], label=r'$\int S_{p,i}$')
            axs.plot(self.mid_axis, self.Q_MW_e['tot'][i], 'o-', color=self.cool_map[i+self.N_cumulative_last], fillstyle='none', label='$Q_e$ (total)')
            marker = itertools.cycle(('d', 's', 'v', 'P', '*', 'X'))
            for label in self.flux_model_labels:
                axs.plot(self.mid_axis, self.Q_MW_e[label][i], '-', marker=next(marker), color=self.cool_map[i+self.N_cumulative_last], fillstyle='none', label=f'$Q_e$ ({label})')
            axs.plot(self.mid_axis, self.Sp_tot_cumint_MW_e[i][:-1], 'o:', color=self.cool_map[i+self.N_cumulative_last], label=r'$\int S_{p,e}$')
            _, top = axs.get_ylim()
            axs.set_ylim(bottom=0, top=2 * top)
            axs.set_xlabel(self.rad_label)
            axs.set_xlim(left=0.0)
            axs.set_title(r'$\bf d)$ power balance [MW]')
            leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(0.65)
            axs.grid(self.grid_lines)
        else:
            i = self.N - 1
            axs.plot(self.mid_axis, self.Q_MW_i['tot'][i], 'o-', color=self.warm_map[i+self.N_cumulative_last], fillstyle='none')
            marker = itertools.cycle(('d', 's', 'v', 'P', '*', 'X'))
            for label in self.flux_model_labels:
                axs.plot(self.mid_axis, self.Q_MW_i[label][i], '-', marker=next(marker), color=self.warm_map[i+self.N_cumulative_last], fillstyle='none')
            axs.plot(self.mid_axis, self.Sp_tot_cumint_MW_i[i][:-1], 'o:', color=self.warm_map[i+self.N_cumulative_last])
            axs.plot(self.mid_axis, self.Q_MW_e['tot'][i], 'o-', color=self.cool_map[i+self.N_cumulative_last], fillstyle='none')
            marker = itertools.cycle(('d', 's', 'v', 'P', '*', 'X'))
            for label in self.flux_model_labels:
                axs.plot(self.mid_axis, self.Q_MW_e[label][i], '-', marker=next(marker), color=self.cool_map[i+self.N_cumulative_last], fillstyle='none')
            axs.plot(self.mid_axis, self.Sp_tot_cumint_MW_e[i][:-1], 'o:', color=self.cool_map[i+self.N_cumulative_last])
            axs.set_xlabel(self.rad_label)
            axs.set_xlim(left=0.0)
            axs.grid(self.grid_lines)

    def plot_heat_flux(self, axs):
        ''' Function to plot the heat flux '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            axs.plot(self.mid_axis, self.Qi['tot'][i], 'x-', color=self.warm_map[i+self.N_cumulative_last])
            axs.plot(self.mid_axis, self.Qe['tot'][i], 'x-', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title(r'$\bf e)$ heat flux [GB]')
        axs.grid(self.grid_lines)

    def plot_particle_balance(self, axs):
        ''' Function to plot the particle balance '''
        if self.global_last:
            i = self.N - 1
            axs.plot(self.mid_axis, self.Gam_SI20_e['tot'][i], 'o-', color=self.green_map[i+self.N_cumulative_last], fillstyle='none', label=r'$\Gamma$')
            axs.plot(self.mid_axis, self.Sn_tot_cumint_SI20_e[i][:-1], 'o:', color=self.green_map[i+self.N_cumulative_last], label=r'$\int S_n$')
            # axs.plot(self.mid_axis, self.Gam_SI20_e[i], '.-', color=self.cool_map[i+self.N_cumulative_last], label=r'$\Gamma_e$')
            # axs.plot(self.axis, self.Sp_tot_int_MW_e[i], '.:', color=self.cool_map[i+self.N_cumulative_last])
            axs.set_xlabel(self.rad_label)
            axs.set_xlim(left=0.0)
            axs.set_title(r'$\bf f)$ particle balance [10$^{20}$/s]')
            leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(0.65)
            axs.grid(self.grid_lines)
        else:            
            i = self.N - 1
            axs.plot(self.mid_axis, self.Gam_SI20_e['tot'][i], 'o-', color=self.green_map[i+self.N_cumulative_last], fillstyle='none')
            axs.plot(self.mid_axis, self.Sn_tot_cumint_SI20_e[i][:-1], 'o:', color=self.green_map[i+self.N_cumulative_last])
            # axs.plot(self.mid_axis, self.Gam_SI20_e[i], '.-', color=self.cool_map[i+self.N_cumulative_last], label=r'$\Gamma_e$')
            # axs.plot(self.axis, self.Sp_tot_int_MW_e[i], '.:', color=self.cool_map[i+self.N_cumulative_last])
            axs.set_xlabel(self.rad_label)
            axs.set_xlim(left=0.0)
            axs.set_title(r'$\bf f)$ particle balance [10$^{20}$/s]')
            axs.grid(self.grid_lines)

    def plot_particle_flux(self, axs):
        ''' Function to plot the particle flux '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            axs.plot(self.mid_axis, self.Gamma['tot'][i], 'x-', color=self.green_map[i+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title(r'$\bf g)$ particle flux [GB]')
        axs.grid(self.grid_lines)

    def plot_fusion_power_density(self, axs):
        ''' Function fusion power density '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            if (i == self.N - 1) and self.global_last:
                axs.plot(self.axis, self.Sp_alpha_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_alpha_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last], label='$S_{p,e}$')
                axs.plot(self.axis, self.Palpha_MWm3[i], '.-', color=self.green_map[i+self.N_cumulative_last], label='$S_{p,{tot}}$')
            else:
                axs.plot(self.axis, self.Sp_alpha_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_alpha_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last])
                axs.plot(self.axis, self.Palpha_MWm3[i], '.-', color=self.green_map[i+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title(r'$\bf h)$ fusion power density [MW/m$^{3}$]')
        title = (r'$P_{\alpha} = $' + f'{self.Palpha_int_MW[i]:.2f} MW\n' +
                 r'$P_{fus} = $' + f'{5.03 * self.Palpha_int_MW[i]:.2f} MW')
        leg = axs.legend(loc='best', title=title, fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_radiation_profiles(self, axs):
        ''' Function to plot the radiation profiles '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            if i == self.N - 1:
                axs.plot(self.axis, self.Sp_rad_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_rad_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last], label='$S_{p,e}$')
            else:
                axs.plot(self.axis, self.Sp_rad_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_rad_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_title(r'$\bf i)$ radiation [MW/m$^{3}$]')
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        try:
            leg = axs.legend(loc='best', title=r'$P_{rad} = $' + f'{self.Prad_tot[-1]:.2f} MW', fancybox=False, shadow=False, ncol=1, fontsize=8)
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(0.65)
        except:
            pass
        axs.grid(self.grid_lines)

    def plot_collisional_equilibration(self, axs):
        ''' Function to plot the collisional equilibration '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            if (i == self.N - 1) and self.global_last:
                axs.plot(self.axis, self.Sp_coll_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_coll_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last], label='$S_{p,e}$')
            else:
                axs.plot(self.axis, self.Sp_coll_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_coll_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_title(r'$\bf h)$ collisional equilibration [MW/m$^{3}$]')
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_auxiliary_power_source(self, axs):
        ''' Function to plot auxiliary power source, convert from Trinity units to MW/m3 '''
        if self.global_last:
            i = self.N - 1
            axs.plot(self.axis, self.Sp_aux_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last], label='$S_{p,i}$')
            axs.plot(self.axis, self.Sp_aux_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last], label='$S_{p,e}$')
            _, top = axs.get_ylim()
            axs.set_ylim(bottom=0, top=1.5 * top)
            axs.set_xlim(left=0.0)
            axs.set_title(r'$\bf i)$ auxiliary power source [MW/m$^{3}$]')
            axs.set_xlabel(self.rad_label)
            leg = axs.legend(loc='best', title='$P_{aux} =$' + f'{self.Sp_aux_int_MW_tot[i]:.2f} MW', fancybox=False, shadow=False, ncol=1, fontsize=8)
            leg.get_frame().set_edgecolor('k')
            leg.get_frame().set_linewidth(0.65)
            axs.grid(self.grid_lines)
        else:
            i = self.N - 1
            axs.plot(self.axis, self.Sp_aux_i[i] * self.P_ref_MWm3, '.-', color=self.warm_map[i+self.N_cumulative_last])
            axs.plot(self.axis, self.Sp_aux_e[i] * self.P_ref_MWm3, '.-', color=self.cool_map[i+self.N_cumulative_last])
            axs.set_xlim(left=0.0)
            axs.set_title(r'$\bf i)$ auxiliary power source [MW/m$^{3}$]')
            axs.set_xlabel(self.rad_label)
            axs.grid(self.grid_lines)
        

    def plot_q(self, axs):
        ''' Function to plot q '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            axs.plot(self.aLpi[i] - self.aLn[i], self.Qi['tot'][i], '.', color=self.warm_map[i+self.N_cumulative_last])
            axs.plot(self.aLpe[i] - self.aLn[i], self.Qe['tot'][i], '.', color=self.cool_map[i+self.N_cumulative_last])
        axs.set_title(r'$\bf j)$ $Q(L_T)$ [GB]')
        axs.set_xlabel('$a/L_T$')
        # leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        # leg.get_frame().set_edgecolor('k')
        # leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_auxiliary_particle_source(self, axs):
        ''' Function to plot the auxiliary particle source '''
        i = self.N - 1
        axs.plot(self.axis, self.Sn_aux_e[i] * self.Sn_ref_SI20, '.-', color=self.green_map[i+self.N_cumulative_last], label='$S_{n}$')
        _, top = axs.get_ylim()
        if self.global_last:
            axs.set_ylim(bottom=0, top=1.5 * top)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'$\bf k)$ auxiliary particle source [10$^{20}$/(m$^{3}$ s)]')
        axs.grid(self.grid_lines)

    def plot_gamma(self, axs):
        ''' Function to plot gamma '''
        for i, t_curr in enumerate(self.time):
            if i > 0:
                if t_curr == self.time[i - 1]:
                    # skip repeated time from Newton iteration
                    continue
            axs.plot(self.aLn[i], self.Gamma['tot'][i], '.', color=self.green_map[i+self.N_cumulative_last])
        axs.set_xlabel('$a/L_n$')
        axs.set_title(r'$\bf l)$ $\Gamma(L_n)$ [GB]')
        axs.grid(self.grid_lines)

    def plot_state_profiles(self, axs, profile='Ti', t_stop=-1):
        ''' Function to plot the density or temperature profiles '''

        if t_stop < 0:
            t_stop = self.N

        t_old = -1
        for ti in np.arange(t_stop):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            # plot profiles
            if profile == 'Ti':
                if t == 0 or t == t_stop-1:
                    axs.plot(self.axis, self.Ti[t], '.-', color=self.warm_map[t], label=f'$T_i$, $t =${self.time[t]:.2f} s')
                else:
                    axs.plot(self.axis, self.Ti[t], '.-', color=self.warm_map[t])
                axs.set_title('temperature [keV]')

            if profile == 'Te':
                if t == 0 or t == t_stop-1:
                    axs.plot(self.axis, self.Te[t], '.-', color=self.cool_map[t], label=f'$T_e$, $t =${self.time[t]:.2f} s')
                else:
                    axs.plot(self.axis, self.Te[t], '.-', color=self.cool_map[t])
                axs.set_title('temperature [keV]')

            if profile == 'ne':
                if t == 0 or t == t_stop-1:
                    axs.plot(self.axis, self.n[t], '.-', color=self.green_map[t], label=f'$n_e$, $t =${self.time[t]:.2f} s')
                else:
                    axs.plot(self.axis, self.n[t], '.-', color=self.green_map[t])
                axs.set_title(r'density [10$^{20}$ m$^{-3}$]')
            
            if profile == 'ni':
                if t == 0 or t == t_stop-1:
                    axs.plot(self.axis, self.ni[t], '.-', color=self.green_map[t], label=f'$n_e$, $t =${self.time[t]:.2f} s')
                else:
                    axs.plot(self.axis, self.ni[t], '.-', color=self.green_map[t])
                axs.set_title(r'density [10$^{20}$ m$^{-3}$]')

        _, top = axs.get_ylim()
        axs.set_ylim(bottom=0, top=1.5*top)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_flux(self, axs, profile='Ti', t_stop=-1):
        ''' Function to plot fluxes from GX'''

        if t_stop < 0:
            t_stop = self.N

        t_old = -1
        for ti in np.arange(t_stop):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti

            if profile == 'Ti':
                axs.plot(self.aLpi[t] - self.aLn[t], self.Qi['tot'][t], '.', color=self.warm_map[t])
                axs.set_title(r'$Q(L_{T_i})$ [GB]')
                axs.set_xlabel('$a/L_{{T_i}}$')

            if profile == 'Te':
                axs.plot(self.aLpe[t] - self.aLn[t], self.Qe['tot'][t], '.', color=self.cool_map[t])
                axs.set_title(r'$Q(L_{T_e})$ [GB]')
                axs.set_xlabel('$a/L_{{T_e}}$')

            if profile == 'ne':
                axs.plot(self.aLn[t], self.Gamma[t], '.', color=self.green_map[t])
                axs.set_xlabel('$a/L_n$')
                axs.set_title(r'$\Gamma(L_n)$ [GB]')
        # leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        # leg.get_frame().set_edgecolor('k')
        # leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_source_total(self, axs, profile='Ti', t=-1):
        ''' Plots total source and sink terms, for each profile '''

        if t < 0:
            t = self.N - 1
        elif t == 0:
            return
        else:
            t = t-1

        if profile == 'Ti':

            Spi_aux = self.Sp_aux_i[t] * self.P_ref_MWm3
            Spi_tot = self.Sp_tot_int_MW_i[t]  # for some reason Sp_aux is excluded?
            axs.plot(self.axis, Spi_tot + Spi_aux, 'x--', color=self.warm_map[t], label=r'$S_{p,i}^{tot}$')

            axs.plot(self.mid_axis, self.Q_MW_i['tot'], 'o:', color=self.warm_map[t], fillstyle='none', label='$-Q_i^{tot}$')
            axs.set_title('ion total source [MW]')

        if profile == 'Te':
            Spe_aux = self.Sp_aux_e[t] * self.P_ref_MWm3
            Spe_tot = self.Sp_tot_int_MW_e[t]  # for some reason Sp_aux is excluded?
            axs.plot(self.axis, Spe_tot + Spe_aux, 'x--', color=self.cool_map[t], label=r'$S_{p,e}^{tot}$')

            axs.plot(self.mid_axis, self.Q_MW_e[t], 'o:', color=self.cool_map[t], fillstyle='none', label='$-Q_e^{tot}$')
            axs.set_title('electron total source [MW]')

        if profile == 'ne':
            axs.plot(self.mid_axis, self.Sn_tot_int_SI20_e[t][:-1], 'x--', color=self.green_map[t], label=r'$S_n$')
            axs.plot(self.mid_axis, self.Gam_SI20_e[t], 'o:', color=self.green_map[t], fillstyle='none', label=r'$-\Gamma$')
            axs.set_title('ne total source [10$^{20}$/s]')

        _, top = axs.get_ylim()
        axs.set_ylim(bottom=0, top=2*top)
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_sink_terms(self, axs, profile='Ti', t=-1):
        ''' Plots source and sink term-by-term, for each profile '''

        if t < 0:
            t = self.N - 1
        elif t == 0:
            return
        else:
            t = t-1

        if profile == 'Ti':
            axs.plot(self.mid_axis, self.Q_turb_MW_i[t], 's--', color=self.warm_map[t], fillstyle='full', label=r'$Q_i^{\mathrm{turb}}$')
            if np.any(self.Q_neo_MW_i[t] > 1e-15):
                axs.plot(self.mid_axis, self.Q_neo_MW_i[t], 'd--', color=self.warm_map[t], fillstyle='full', label=r'$Q_i^{\mathrm{neo}}$')

            if np.any(self.Q_alpha_int_MW_i[t] > 1e-15):
                axs.plot(self.axis, self.Q_alpha_int_MW_i[t], 'H--', color=self.warm_map[t], fillstyle='none', label=r'$-S_{p,i}^{\mathrm{alpha}}$')

            if np.any(self.Q_rad_int_MW_i[t] > 1e-15):
                axs.plot(self.axis, self.Q_rad_int_MW_i[t], '^--', color=self.warm_map[t], fillstyle='full', label=r'$Q_{p,i}^{\mathrm{rad}}$')

            if np.mean(self.Q_coll_int_MW_i[t] > 0):
                axs.plot(self.axis, self.Q_coll_int_MW_i[t], 'v--', color=self.cool_map[t], fillstyle='none', label=r'$-S_{p,i}^{\mathrm{coll}}$')
            elif np.mean(self.Q_coll_int_MW_i[t] < 0):
                axs.plot(self.axis, -self.Q_coll_int_MW_i[t], 'v--', color=self.cool_map[t], fillstyle='full', label=r'$Q_{p,i}^{\mathrm{coll}}$')

            axs.plot(self.axis, self.Sp_aux_i[t] * self.P_ref_MWm3, 'x--', color=self.warm_map[t], label=r'$-S_{p,i}^{\mathrm{aux}}$')

            axs.set_title('ion sink terms [MW]')

        if profile == 'Te':
            axs.plot(self.mid_axis, self.Q_turb_MW_e[t], 's--', color=self.cool_map[t], fillstyle='full', label=r'$Q_e^{\mathrm{turb}}$')
            if np.any(self.Q_neo_MW_e[t] > 1e-15):
                axs.plot(self.mid_axis, self.Q_neo_MW_e[t], 'd--', color=self.cool_map[t], fillstyle='full', label=r'$Q_e^{\mathrm{neo}}$')

            if np.any(self.Q_rad_int_MW_e[t] != 0):
                axs.plot(self.axis, -self.Q_rad_int_MW_e[t], '^--', color=self.cool_map[t], fillstyle='full', label=r'$Q_{p,e}^{\mathrm{rad}}$')

            if np.any(self.Q_alpha_int_MW_e[t] != 0):
                axs.plot(self.axis, self.Q_alpha_int_MW_e[t], 'H--', color=self.cool_map[t], fillstyle='none', label=r'$-S_{p,e}^{\mathrm{alpha}}$')

            if np.any(self.Q_coll_int_MW_e[t] > 1e-15):
                axs.plot(self.axis, self.Q_coll_int_MW_e[t], 'v--', color=self.cool_map[t], fillstyle='none', label=r'$-S_{p,e}^{\mathrm{coll}}$')

            axs.plot(self.axis, self.Sp_aux_e[t] * self.P_ref_MWm3, 'x--', color=self.cool_map[t], label=r'$-S_{p,e}^{\mathrm{aux}}$')

            axs.set_title('electron sink terms [MW]')

        if profile == 'ne':
            axs.plot(self.mid_axis, self.Gam_SI20_e[t], 'o-', color=self.green_map[t], fillstyle='none', label=r'$\Gamma^{\mathrm{turb}}$')
            axs.set_title('ne sink terms [10$^{20}$/s]')

        _, top = axs.get_ylim()
        axs.set_xlabel(self.rad_label)
        # axs.set_ylim(bottom=0, top=2*top)
        axs.set_xlim(left=0.0)
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def make_plots(self, plots):
        ''' Create requested plots '''

        self.plot_electrons = True

        for plot in plots:
            _, axs = plt.subplots(1, 1)
            plotting_function = self.plotting_functions.get(plot.lower())
            plotting_function(axs)

        plt.show()

    def plot_panel(self):
        ''' Create the legacy plot panel '''

        rlabel = rf'[{self.filename}]'
        _, axs = plt.subplots(2, 7, figsize=(65, 8))

        plt.suptitle(rlabel)

        self.plot_electrons = True

        for i, (_, plot) in enumerate(self.panel_plots.items()):
            plot(axs[int(i / 7), i % 7])

        # Create the plot panel
        plt.tight_layout()
        plt.subplots_adjust(left=0.05,
                            bottom=0.1,
                            right=0.95,
                            top=0.9,
                            wspace=0.4,
                            hspace=0.4)

        plt.show()

    def plot_panel_MY(self, axs):
        ''' Create the legacy plot panel '''

        self.plot_electrons = True

        for i, (_, plot) in enumerate(self.panel_plots.items()):
            plot(axs[int(i / axs.shape[1]), i % axs.shape[1]])

        # # Create the plot panel
        # plt.tight_layout()
        # plt.subplots_adjust(left=0.05,
        #                     bottom=0.1,
        #                     right=0.95,
        #                     top=0.9,
        #                     wspace=0.4,
        #                     hspace=0.4)



def main():

    import argparse

    plots = plotter()

    description = "This program plots Trinity3D .npy or .nc output file"
    epilog = '''\
Usage: to create the plot panel
  t3d-plot [trinity.log.npy]
  t3d-plot [trinity.nc]

Usage: to create indvidual plots
  t3d-plot [trinity.log.npy] -p density temperature etc
  t3d-plot [trinity.nc] -p density temperature etc
'''
    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("filename",
                        help="Specify plot file ")
    parser.add_argument("-p", "--plots",
                        nargs="+",
                        help="Space separated list of plots",
                        default=None)
    parser.add_argument("-l", "--list",
                        help="List available plot routines",
                        action="store_true")
    parser.add_argument("-g", "--grid",
                        help="Include grid lines in plot",
                        action="store_true")

    args = parser.parse_args()

    if args.list:
        plots.list_available_plots()
        exit(0)

    if args.grid:
        plots.grid_lines = True

    if args.filename.split('.')[-1] == "npy":
        plots.read_data(args.filename)
    else:
        plots.read_netcdf(args.filename)

    if args.plots is None:
        plots.plot_panel()
    else:
        plots.make_plots(args.plots)


if __name__ == '__main__':
    '''
        This program plots Trinity3D output read from the .npy output
        file and stored in a python dictionary.
    '''

    main()
