import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

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

        # Dictionary of plotting functions
        self.plotting_functions = {
            "density": self.plot_density_profiles,
            "temperature": self.plot_temperature_profiles,
            "pressure": self.plot_pressure_profiles,
            "power_balance": self.plot_power_balance,
            "heat_flux": self.plot_heat_flux,
            "particle_balance": self.plot_particle_balance,
            "particle_flux": self.plot_particle_flux,
            "fusion_power_density": self.plot_fusion_power_density,
            "radiation": self.plot_radiation_profiles,
            "collisional_equilibration": self.plot_collisional_equilibration,
            "power_source": self.plot_auxiliary_power_source,
            "q": self.plot_q,
            "particle_source": self.plot_auxiliary_particle_source,
            "gamma": self.plot_gamma,
        }

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
        self.rho_star = 1 / data['norms']['gyro_scale']
        self.t_scale = self.t_ref * self.rho_star**-2
        self.time = np.array(data['t'])*self.t_scale
        self.tags = data['species_tags']
        self.ion_tag = data['bulk_ion_tag']

        self.N = len(self.time)

        # set up color maps
        self.warm_map = pylab.cm.autumn(np.linspace(1, 0.25, self.N))
        self.cool_map = pylab.cm.Blues(np.linspace(0.25, 1, self.N))
        self.green_map = pylab.cm.YlGn(np.linspace(0.25, 1, self.N))
        # purple_map = pylab.cm.Purples(np.linspace(0.25, 1, self.N))

        self.n = np.asarray(data['n_e'])
        self.pi = np.asarray(data['p_' + self.ion_tag])
        self.pe = np.asarray(data['p_e'])
        self.Ti = np.asarray(data['T_' + self.ion_tag])
        self.Te = np.asarray(data['T_e'])
        self.Gamma = np.asarray(data['pflux_' + self.ion_tag])
        self.Qi = np.asarray(data['qflux_' + self.ion_tag])
        self.Qe = np.asarray(data['qflux_e'])
        self.Gamma_neo = np.asarray(data['pflux_neo_' + self.ion_tag])
        self.Qi_neo = np.asarray(data['qflux_neo_' + self.ion_tag])
        self.Qe_neo = np.asarray(data['qflux_neo_e'])
        self.Gamma_turb = np.asarray(data['pflux_turb_' + self.ion_tag])
        self.Qi_turb = np.asarray(data['qflux_turb_' + self.ion_tag])
        self.Qe_turb = np.asarray(data['qflux_turb_e'])
        self.aLn = np.asarray(data['aLn_' + self.ion_tag])
        self.aLpi = np.asarray(data['aLp_' + self.ion_tag])
        self.aLpe = np.asarray(data['aLp_e'])
        self.aLTi = np.asarray(data['aLT_' + self.ion_tag])
        self.aLTe = np.asarray(data['aLT_e'])

        # Convert pressure sources from code units to MW/m^3 by dividing by p_source_scale
        self.p_source_scale = data['norms']['pressure_source_scale']*1e6
        self.P_ref_MWm3 = data['norms']['P_ref_MWm3']
        self.Sn_ref_SI20 = data['norms']['Sn_ref_SI20']
        self.a_ref = data['norms']['a_ref']

        # Pressure sources are in code units, multiply by P_ref_MWm3 to get MW/m^3
        self.Sp_aux_i = np.asarray(data['Sp_aux_' + self.ion_tag])
        self.Sp_aux_e = np.asarray(data['Sp_aux_e'])
        self.Sp_aux_int_MW_i = np.asarray(data['Sp_aux_int_MW_' + self.ion_tag])
        self.Sp_aux_int_MW_e = np.asarray(data['Sp_aux_int_MW_e'])
        self.Sp_aux_int_MW_tot = self.Sp_aux_int_MW_i + self.Sp_aux_int_MW_e

        self.Sn_aux_i = np.asarray(data['Sn_aux_' + self.ion_tag])
        self.Sn_aux_e = np.asarray(data['Sn_aux_e'])

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

        self.Palpha_MWm3 = np.asarray(data['Palpha_MWm3'])
        self.Palpha_int_MW = data['Palpha_int_MW']

        # power balance
        self.Q_MW_i = np.asarray(data['Q_MW_' + self.ion_tag])
        self.Q_neo_MW_i = np.asarray(data['Q_neo_MW_' + self.ion_tag])
        self.Q_turb_MW_i = np.asarray(data['Q_turb_MW_' + self.ion_tag])
        self.Sp_tot_int_MW_i = np.asarray(data['Sp_tot_int_MW_' + self.ion_tag])
        self.Q_MW_e = np.asarray(data['Q_MW_e'])
        self.Q_neo_MW_e = np.asarray(data['Q_neo_MW_e'])
        self.Q_turb_MW_e = np.asarray(data['Q_turb_MW_e'])
        self.Sp_tot_int_MW_e = np.asarray(data['Sp_tot_int_MW_e'])

        # particle balance
        self.Gam_SI20_i = np.asarray(data['Gam_SI20_' + self.ion_tag])
        self.Gam_SI20_e = np.asarray(data['Gam_SI20_e'])
        self.Sn_tot_int_SI20_i = np.asarray(data['Sn_tot_int_SI20_' + self.ion_tag])
        self.Sn_tot_int_SI20_e = np.asarray(data['Sn_tot_int_SI20_e'])

        self.grid = data['grid']
        self.N_rho = self.grid['N_radial']
        self.rho_edge = self.grid['rho_edge']
        self.axis = self.grid['rho_axis']
        self.mid_axis = self.grid['mid_axis']
        self.drho = self.axis[1] - self.axis[0]
        try:
            self.flux_label = self.grid['flux_label']
        except LookupError:
            self.flux_label = 'rminor'
        if self.flux_label == 'rminor':
            self.rad_label = r'$r/a$'
        elif self.flux_label == 'torflux':
            self.rad_label = r'$\rho_{tor}$'

        geo = data['geo_obj']
        self.area = geo.area.profile
        self.grho = geo.grho.profile
        self.area_grid = geo.area.toGridProfile().profile
        self.grho_grid = geo.grho.toGridProfile().profile
        self.G_fac = geo.G_fac.profile
        self.B_ref = geo.B_ref.profile

        self.grid = data['grid_obj']

        self.Prad_tot = geo.volume_integrate(Profiles.GridProfile(self.Sp_rad_e[-1], self.grid))*self.P_ref_MWm3

        self.dtau = data['time']['dtau']

        # jacobian = np.round(data['solver_jacobian'][0]*dtau, 14)
        # rhs = np.round(data['solver_rhs'][0]*dtau, 14)
        # print("JACOBIAN")
        # print(jacobian)
        # print("RHS")
        # print(rhs)

        # run settings
        self.alpha = data['time']['alpha']
        self.dtau = data['time']['dtau']
        self.N_steps = data['time']['N_steps']

        self.filename = data['trinity_infile']

    def plot_density_profiles(self, axs):
        ''' Function to plot the density profiles '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            # plot profiles
            if (t == 0 and self.global_first) or (t == self.N - 1 and self.global_last):
                axs.plot(self.axis, self.n[t], '.-', color=self.green_map[t+self.N_cumulative_last],
                         label=f'$n_e$, $t =${self.time[t]+self.time_cumulative_last:.2f} s')
            else:
                axs.plot(self.axis, self.n[t], '.-', color=self.green_map[t+self.N_cumulative_last])
        nmax = np.max(self.n)
        axs.set_ylim(bottom=0, top=1.5*nmax)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'density [10$^{20}$ m$^{-3}$]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_temperature_profiles(self, axs):
        ''' Function to plot the temperature profiles '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            # plot profiles
            if t == 0:
                axs.plot(self.axis, self.Ti[t], '.-', color=self.warm_map[t+self.N_cumulative_last], label=f'$T_i$, $t =${self.time[t]:.2f} s')
                axs.plot(self.axis, self.Te[t], '.-', color=self.cool_map[t+self.N_cumulative_last], label=f'$T_e$, $t =${self.time[t]:.2f} s')
            elif t == self.N-1:
                axs.plot(self.axis, self.Ti[t], '.-', color=self.warm_map[t+self.N_cumulative_last], label=f'$T_i$, $t =${self.time[t]:.2f} s')
                if self.plot_electrons:
                    axs.plot(self.axis, self.Te[t], '.-', color=self.cool_map[t+self.N_cumulative_last], label=f'$T_e$, $t =${self.time[t]:.2f} s')
            else:
                axs.plot(self.axis, self.Ti[t], '.-', color=self.warm_map[t+self.N_cumulative_last])
                if self.plot_electrons:
                    axs.plot(self.axis, self.Te[t], '.-', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_ylim(bottom=0)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title('temperature [keV]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_pressure_profiles(self, axs):
        ''' Function to plot the pressure profiles '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            axs.plot(self.axis, self.pi[t], '.-', color=self.warm_map[t+self.N_cumulative_last])
            if self.plot_electrons:
                axs.plot(self.axis, self.pe[t], '.-', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_ylim(bottom=0)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title(r'pressure [10$^{20}$m$^{-3}$ keV]')
        axs.grid(self.grid_lines)

    def plot_power_balance(self, axs):
        ''' Function to plot the power balance '''
        t = self.N - 1
        axs.plot(self.mid_axis, self.Q_MW_i[t], 'o-', color=self.warm_map[t+self.N_cumulative_last], fillstyle='none', label='$Q_i$')
        if np.any(self.Q_neo_MW_i[t] > 1e-15):
            axs.plot(self.mid_axis, self.Q_neo_MW_i[t], 'd--', color=self.warm_map[t+self.N_cumulative_last], fillstyle='none', label=r'$Q_i^{\mathrm{neo}}$')
            axs.plot(self.mid_axis, self.Q_turb_MW_i[t], 's--', color=self.warm_map[t+self.N_cumulative_last], fillstyle='none', label=r'$Q_i^{\mathrm{turb}}$')
        axs.plot(self.mid_axis, self.Sp_tot_int_MW_i[t][:-1], 'o:', color=self.warm_map[t+self.N_cumulative_last], label=r'$\int S_{p,i}$')
        axs.plot(self.mid_axis, self.Q_MW_e[t], 'o-', color=self.cool_map[t+self.N_cumulative_last], fillstyle='none', label='$Q_e$')
        if np.any(self.Q_neo_MW_e[t] > 1e-15):
            axs.plot(self.mid_axis, self.Q_neo_MW_e[t], 'd--', color=self.cool_map[t+self.N_cumulative_last], fillstyle='none', label=r'$Q_e^{\mathrm{neo}}$')
            axs.plot(self.mid_axis, self.Q_turb_MW_e[t], 's--', color=self.cool_map[t+self.N_cumulative_last], fillstyle='none', label=r'$Q_e^{\mathrm{turb}}$')
        axs.plot(self.mid_axis, self.Sp_tot_int_MW_e[t][:-1], 'o:', color=self.cool_map[t+self.N_cumulative_last], label=r'$\int S_{p,e}$')
        _, top = axs.get_ylim()
        axs.set_ylim(bottom=0, top=2*top)
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title('power balance [MW]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_heat_flux(self, axs):
        ''' Function to plot the heat flux '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            axs.plot(self.mid_axis, self.Qi[t], 'x-', color=self.warm_map[t+self.N_cumulative_last])
            axs.plot(self.mid_axis, self.Qe[t], 'x-', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title('heat flux [GB]')
        axs.grid(self.grid_lines)

    def plot_particle_balance(self, axs):
        ''' Function to plot the particle balance '''
        t = self.N - 1
        axs.plot(self.mid_axis, self.Gam_SI20_e[t], 'o-', color=self.green_map[t+self.N_cumulative_last], fillstyle='none', label=r'$\Gamma$')
        axs.plot(self.mid_axis, self.Sn_tot_int_SI20_e[t][:-1], 'o:', color=self.green_map[t+self.N_cumulative_last], label=r'$\int S_n$')
        # axs.plot(self.mid_axis, self.Gam_SI20_e[t], '.-', color=self.cool_map[t+self.N_cumulative_last], label=r'$\Gamma_e$')
        # axs.plot(self.axis, self.Sp_tot_int_MW_e[t], '.:', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title('particle balance [10$^{20}$/s]')
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_particle_flux(self, axs):
        ''' Function to plot the particle flux '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            axs.plot(self.mid_axis, self.Gamma[t], 'x-', color=self.green_map[t+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title('particle flux [GB]')
        axs.grid(self.grid_lines)

    def plot_fusion_power_density(self, axs):
        ''' Function fusion power density '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            if t == self.N-1:
                axs.plot(self.axis, self.Sp_alpha_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_alpha_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last], label='$S_{p,e}$')
                axs.plot(self.axis, self.Palpha_MWm3[t], '.-', color=self.green_map[t+self.N_cumulative_last], label='$S_{p,{tot}}$')
            else:
                axs.plot(self.axis, self.Sp_alpha_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_alpha_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last])
                axs.plot(self.axis, self.Palpha_MWm3[t], '.-', color=self.green_map[t+self.N_cumulative_last])
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        axs.set_title('fusion power density \n [MW/m$^{3}$]')
        title = (r'$P_{\alpha} = $' + f'{self.Palpha_int_MW[t]:.2f} MW\n' +
                 r'$P_{fus} = $' + f'{5.03*self.Palpha_int_MW[t]:.2f} MW')
        leg = axs.legend(loc='best', title=title, fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_radiation_profiles(self, axs):
        ''' Function to plot the radiation profiles '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            if t == self.N-1:
                axs.plot(self.axis, self.Sp_rad_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_rad_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last], label='$S_{p,e}$')
            else:
                axs.plot(self.axis, self.Sp_rad_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_rad_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_title('radiation [MW/m$^{3}$]')
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        leg = axs.legend(loc='best', title=r'$P_{rad} = $' + f'{self.Prad_tot:.2f} MW', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_collisional_equilibration(self, axs):
        ''' Function to plot the collisional equilibration '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            if t == self.N-1:
                axs.plot(self.axis, self.Sp_coll_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last], label='$S_{p,i}$')
                axs.plot(self.axis, self.Sp_coll_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last], label='$S_{p,e}$')
            else:
                axs.plot(self.axis, self.Sp_coll_i[t]*self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last])
                axs.plot(self.axis, self.Sp_coll_e[t]*self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_title('collisional equilibration \n [MW/m$^{3}$]')
        axs.set_xlabel(self.rad_label)
        axs.set_xlim(left=0.0)
        leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_auxiliary_power_source(self, axs):
        ''' Function to plot auxiliary power source, convert from Trinity units to MW/m3 '''
        t = self.N - 1
        axs.plot(self.axis, self.Sp_aux_i[t] * self.P_ref_MWm3, '.-', color=self.warm_map[t+self.N_cumulative_last], label='$S_{p,i}$')
        axs.plot(self.axis, self.Sp_aux_e[t] * self.P_ref_MWm3, '.-', color=self.cool_map[t+self.N_cumulative_last], label='$S_{p,e}$')
        _, top = axs.get_ylim()
        axs.set_ylim(bottom=0, top=1.5*top)
        axs.set_xlim(left=0.0)
        axs.set_title('auxiliary power source \n [MW/m$^{3}$]')
        axs.set_xlabel(self.rad_label)
        leg = axs.legend(loc='best', title='$P_{aux} =$' + f'{self.Sp_aux_int_MW_tot[t]:.2f} MW', fancybox=False, shadow=False, ncol=1, fontsize=8)
        leg.get_frame().set_edgecolor('k')
        leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_q(self, axs):
        ''' Function to plot q '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            axs.plot(self.aLpi[t] - self.aLn[t], self.Qi[t], '.', color=self.warm_map[t+self.N_cumulative_last])
            axs.plot(self.aLpe[t] - self.aLn[t], self.Qe[t], '.', color=self.cool_map[t+self.N_cumulative_last])
        axs.set_title(r'$Q(L_{T})$ [GB]')
        axs.set_xlabel('$a/L_{{T}}$')
        # leg = axs.legend(loc='best', fancybox=False, shadow=False, ncol=1, fontsize=8)
        # leg.get_frame().set_edgecolor('k')
        # leg.get_frame().set_linewidth(0.65)
        axs.grid(self.grid_lines)

    def plot_auxiliary_particle_source(self, axs):
        ''' Function to plot the auxiliary power source '''
        t = self.N - 1
        axs.plot(self.axis, self.Sn_aux_e[t] * self.Sn_ref_SI20, '.-', color=self.green_map[t+self.N_cumulative_last], label='$S_{n}$')
        _, top = axs.get_ylim()
        axs.set_ylim(bottom=0, top=1.5*top)
        axs.set_xlim(left=0.0)
        axs.set_xlabel(self.rad_label)
        axs.set_title('auxiliary particle source \n [10$^{20}$/(m$^{3}$ s)]')
        axs.grid(self.grid_lines)

    def plot_gamma(self, axs):
        ''' Function to plot gamma '''
        t_old = -1
        for ti in np.arange(self.N):
            t_curr = self.time[ti]
            if (t_curr == t_old):
                # skip majority of prints
                continue
            else:
                t_old = t_curr
                t = ti
            axs.plot(self.aLn[t], self.Gamma[t], '.', color=self.green_map[t+self.N_cumulative_last])
        axs.set_xlabel('$a/L_n$')
        axs.set_title(r'$\Gamma(L_n)$ [GB]')
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

        for i, (_, plot) in enumerate(self.plotting_functions.items()):
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


def main():

    import argparse

    plots = plotter()

    description = "This program plots Trinity3D .npy output file"
    epilog = '''\
Usage: to create the plot panel
  t3d-plot [trinity-log.npy]

Usage: to create indvidual plots
  t3d-plot [trinity-log.npy] -p density temperature etc
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

    plots.read_data(args.filename)
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
