import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, asdict

# ============================================================
# constants
# ============================================================
c = 2.99792458e8
eps0 = 8.8541878128e-12
h = 6.62607015e-34
hbar = 1.054571817e-34
eV = 1.602176634e-19
kB = 1.380649e-23
m_e = 9.1093837015e-31
e_charge = 1.602176634e-19
pi = np.pi

# ============================================================
# dataclasses
# ============================================================
@dataclass
class LaserConfig:
    wavelength_m: float = 795e-9
    waist_m: float = 22e-6                 # focused 1/e^2 intensity radius
    tau_fwhm_s: float = 30e-15             # intensity FWHM
    avg_power_W: float = 8.0
    rep_rate_Hz: float = 1.0e3


@dataclass
class GasConfig:
    species: str = "Ar"
    pressure_mbar: float = 50.0
    temperature_K: float = 300.0
    medium_length_m: float = 1.2e-3
    jet_center_m: float = 0.0
    profile: str = "gaussian"              # "gaussian" or "tophat"
    sigma_z_m: float = 0.45e-3             # used if profile="gaussian"


@dataclass
class HarmonicConfig:
    q_min: int = 11
    q_max: int = 51
    selected_q_plot: int = 21
    Ip_eV: float = 15.7596                 # argon ionization potential


@dataclass
class GridConfig:
    Nr: int = 140
    Nz: int = 260
    Nt: int = 220
    r_max_m: float = 140e-6
    t_max_s: float = 90e-15


@dataclass
class IonizationConfig:
    model: str = "surrogate"
    I_sat_Wm2: float = 1.2e18
    gamma: float = 6.0
    eta_max: float = 0.98


@dataclass
class DipoleConfig:
    short_weight: float = 1.0
    long_weight: float = 0.35

    # phase coefficient alpha in phi = alpha * I
    alpha_short_m2_per_W_at_q21: float = 0.8e-18
    alpha_long_m2_per_W_at_q21: float = 3.0e-18

    # amplitude response
    p_base: float = 4.0
    p_slope_per_harmonic: float = 0.015
    ion_suppression_beta: float = 10.0

    # spectral rolloff softness near cutoff
    rolloff_width_harmonics: float = 5.0

    # long trajectory usually somewhat weaker after macroscopic averaging
    long_amplitude_factor: float = 0.6


@dataclass
class AbsorptionConfig:
    model: str = "powerlaw"
    L_abs_ref_m_at_q21: float = 0.6e-3
    power_law_exponent: float = 0.5


@dataclass
class CalibrationConfig:
    """
    Maps dimensionless coherent surrogate yield into physical XUV pulse energy.
    This is the one remaining empirical normalization.

    U_q = scale * Y_q * E_pulse
    """
    scale: float = 2.0e-16


# ============================================================
# utility / physics helpers
# ============================================================
def pulse_energy_from_avg_power(avg_power_W, rep_rate_Hz):
    return avg_power_W / rep_rate_Hz


def rayleigh_length(wavelength_m, waist_m):
    return pi * waist_m**2 / wavelength_m


def gaussian_peak_power_from_pulse_energy(E_pulse_J, tau_fwhm_s):
    factor = np.sqrt(pi / (4.0 * np.log(2.0)))
    return E_pulse_J / (tau_fwhm_s * factor)


def gaussian_peak_intensity_from_peak_power(P_peak_W, waist_m):
    return 2.0 * P_peak_W / (pi * waist_m**2)


def electric_field_from_intensity(I_Wm2):
    return np.sqrt(2.0 * I_Wm2 / (c * eps0))


def photon_energy_J_from_q(q, wavelength_m):
    omega0 = 2.0 * pi * c / wavelength_m
    return hbar * q * omega0


def photon_energy_eV_from_q(q, wavelength_m):
    return photon_energy_J_from_q(q, wavelength_m) / eV


def argon_density_m3(pressure_mbar, temperature_K):
    pressure_Pa = pressure_mbar * 100.0
    return pressure_Pa / (kB * temperature_K)


def ponderomotive_energy_eV(I_Wm2, wavelength_m):
    # Up[eV] = 9.337e-14 * I[W/cm^2] * lambda[um]^2
    I_Wcm2 = I_Wm2 / 1e4
    lam_um = wavelength_m * 1e6
    return 9.337e-14 * I_Wcm2 * lam_um**2


def cutoff_energy_eV(I_Wm2, wavelength_m, Ip_eV):
    Up = ponderomotive_energy_eV(I_Wm2, wavelength_m)
    return Ip_eV + 3.17 * Up


def cutoff_harmonic_order(I_Wm2, wavelength_m, Ip_eV):
    Ecut = cutoff_energy_eV(I_Wm2, wavelength_m, Ip_eV)
    Efund = 1240.0 / (wavelength_m * 1e9)
    return Ecut / Efund


def gas_density_profile_z(z, gas: GasConfig, N0_peak):
    if gas.profile == "tophat":
        return np.where(np.abs(z - gas.jet_center_m) <= gas.medium_length_m / 2.0, N0_peak, 0.0)
    elif gas.profile == "gaussian":
        return N0_peak * np.exp(-0.5 * ((z - gas.jet_center_m) / gas.sigma_z_m)**2)
    else:
        raise ValueError(f"Unknown gas profile: {gas.profile}")


def gaussian_intensity_rt(R, T, zpos, I0_Wm2, waist_m, zR_m, tau_fwhm_s):
    wz = waist_m * np.sqrt(1.0 + (zpos / zR_m)**2)
    spatial = (waist_m / wz)**2 * np.exp(-2.0 * R**2 / wz**2)
    temporal = np.exp(-4.0 * np.log(2.0) * T**2 / tau_fwhm_s**2)
    return I0_Wm2 * spatial * temporal


def gouy_phase(zpos, zR_m):
    return np.arctan(zpos / zR_m)


def radius_of_curvature(zpos, zR_m):
    if np.isclose(zpos, 0.0):
        return np.inf
    return zpos * (1.0 + (zR_m / zpos)**2)


def ionization_fraction_surrogate(I_Wm2, ion: IonizationConfig):
    eta = (I_Wm2 / ion.I_sat_Wm2) ** ion.gamma
    return np.clip(eta, 0.0, ion.eta_max)


def absorption_length_q(q, absorption: AbsorptionConfig):
    if absorption.model == "powerlaw":
        return absorption.L_abs_ref_m_at_q21 * (21.0 / q) ** absorption.power_law_exponent
    raise ValueError(f"Unknown absorption model: {absorption.model}")


def q_scaled_alpha(q, alpha_at_q21):
    return alpha_at_q21 * (q / 21.0)


def dipole_phase(I_Wm2, alpha_m2_per_W):
    return alpha_m2_per_W * I_Wm2


def spectral_rolloff(q, q_cut, width):
    if q <= 0.85 * q_cut:
        return 1.0
    return np.exp(-(q - 0.85 * q_cut) / width)


def dipole_amplitude_surrogate(I_Wm2, eta, q, q_cut, dip: DipoleConfig):
    I_norm = I_Wm2 / np.max(I_Wm2)
    p_q = dip.p_base + dip.p_slope_per_harmonic * max(q - 15, 0)
    return (I_norm ** p_q) * np.exp(-dip.ion_suppression_beta * eta) * spectral_rolloff(
        q, q_cut, dip.rolloff_width_harmonics
    )


def build_trajectory_source(
    q,
    I_rt,
    eta_rt,
    Nn_rt,
    R,
    zpos,
    wavelength_m,
    zR_m,
    q_cut,
    dip: DipoleConfig,
):
    """
    Build the local nonlinear polarization source using short and long trajectories.
    """
    amp_base = dipole_amplitude_surrogate(I_rt, eta_rt, q, q_cut, dip)

    alpha_s = q_scaled_alpha(q, dip.alpha_short_m2_per_W_at_q21)
    alpha_l = q_scaled_alpha(q, dip.alpha_long_m2_per_W_at_q21)

    phi_g = q * gouy_phase(zpos, zR_m)

    Rz = radius_of_curvature(zpos, zR_m)
    if np.isfinite(Rz):
        phi_curv = q * (2.0 * pi / wavelength_m) * (R**2 / (2.0 * Rz))
    else:
        phi_curv = np.zeros_like(R)

    phi_s = phi_g + phi_curv + dipole_phase(I_rt, alpha_s)
    phi_l = phi_g + phi_curv + dipole_phase(I_rt, alpha_l)

    P_short = dip.short_weight * amp_base * np.exp(1j * phi_s)
    P_long = dip.long_weight * dip.long_amplitude_factor * amp_base * np.exp(1j * phi_l)

    return Nn_rt * (P_short + P_long)


def coherent_exit_field_for_q(
    q,
    laser: LaserConfig,
    gas: GasConfig,
    ion: IonizationConfig,
    dip: DipoleConfig,
    absorption: AbsorptionConfig,
    grids,
    derived,
):
    r, z, t, R, T = grids
    dz = z[1] - z[0]

    I0 = derived["I0_Wm2"]
    zR = derived["zR_m"]
    q_cut = derived["q_cut"]
    N0_peak = derived["N0_peak_m3"]

    Nz_profile = gas_density_profile_z(z, gas, N0_peak)
    z_exit = z[-1]

    E_q_rt = np.zeros((len(r), len(t)), dtype=np.complex128)
    L_abs_q = absorption_length_q(q, absorption)

    for iz, zpos in enumerate(z):
        N0_local = Nz_profile[iz]
        if N0_local <= 0.0:
            continue

        I_rt = gaussian_intensity_rt(R, T, zpos, I0, laser.waist_m, zR, laser.tau_fwhm_s)
        eta_rt = ionization_fraction_surrogate(I_rt, ion)
        Nn_rt = N0_local * (1.0 - eta_rt)

        Pq_rt = build_trajectory_source(
            q=q,
            I_rt=I_rt,
            eta_rt=eta_rt,
            Nn_rt=Nn_rt,
            R=R,
            zpos=zpos,
            wavelength_m=laser.wavelength_m,
            zR_m=zR,
            q_cut=q_cut,
            dip=dip,
        )

        absorption_factor = np.exp(-(z_exit - zpos) / (2.0 * L_abs_q))
        E_q_rt += Pq_rt * absorption_factor * dz

    return E_q_rt


def surrogate_yield_from_exit_field(E_q_rt, r, t):
    """
    Dimensionless coherent surrogate yield.
    """
    intensity_rt = np.abs(E_q_rt) ** 2
    radial_time_integrand = 2.0 * pi * r[:, None] * intensity_rt
    return np.trapz(np.trapz(radial_time_integrand, r, axis=0), t)


def photons_per_pulse_from_surrogate_yield(Y_q, E_pulse_J, q, wavelength_m, cal: CalibrationConfig):
    U_q_J = cal.scale * Y_q * E_pulse_J
    Eph_J = photon_energy_J_from_q(q, wavelength_m)
    Nphot = U_q_J / Eph_J if Eph_J > 0 else 0.0
    return U_q_J, Nphot


# ============================================================
# main simulator
# ============================================================
class HHGArgonModel:
    def __init__(
        self,
        laser: LaserConfig,
        gas: GasConfig,
        harmonic: HarmonicConfig,
        grid: GridConfig,
        ion: IonizationConfig,
        dip: DipoleConfig,
        absorption: AbsorptionConfig,
        calibration: CalibrationConfig,
    ):
        self.laser = laser
        self.gas = gas
        self.harmonic = harmonic
        self.grid = grid
        self.ion = ion
        self.dip = dip
        self.absorption = absorption
        self.calibration = calibration

        self._build_grids()
        self._build_derived()

    def _build_grids(self):
        self.r = np.linspace(0.0, self.grid.r_max_m, self.grid.Nr)
        self.z = np.linspace(
            self.gas.jet_center_m - self.gas.medium_length_m,
            self.gas.jet_center_m + self.gas.medium_length_m,
            self.grid.Nz,
        )
        self.t = np.linspace(-self.grid.t_max_s, self.grid.t_max_s, self.grid.Nt)
        self.R, self.T = np.meshgrid(self.r, self.t, indexing="ij")

    def _build_derived(self):
        E_pulse = pulse_energy_from_avg_power(self.laser.avg_power_W, self.laser.rep_rate_Hz)
        P_peak = gaussian_peak_power_from_pulse_energy(E_pulse, self.laser.tau_fwhm_s)
        I0 = gaussian_peak_intensity_from_peak_power(P_peak, self.laser.waist_m)
        zR = rayleigh_length(self.laser.wavelength_m, self.laser.waist_m)
        N0_peak = argon_density_m3(self.gas.pressure_mbar, self.gas.temperature_K)
        q_cut = cutoff_harmonic_order(I0, self.laser.wavelength_m, self.harmonic.Ip_eV)
        Ecut = cutoff_energy_eV(I0, self.laser.wavelength_m, self.harmonic.Ip_eV)

        self.derived = {
            "E_pulse_J": E_pulse,
            "P_peak_W": P_peak,
            "I0_Wm2": I0,
            "I0_Wcm2": I0 / 1e4,
            "E0_Vpm": electric_field_from_intensity(I0),
            "zR_m": zR,
            "N0_peak_m3": N0_peak,
            "q_cut": q_cut,
            "Ecut_eV": Ecut,
        }

    def simulate(self):
        q_values = list(range(self.harmonic.q_min, self.harmonic.q_max + 1, 2))

        photons = []
        energies_J = []
        energies_eV = []
        surrogate_yields = []

        selected_exit_field = None
        selected_radial_profile = None

        for q in q_values:
            E_q_rt = coherent_exit_field_for_q(
                q=q,
                laser=self.laser,
                gas=self.gas,
                ion=self.ion,
                dip=self.dip,
                absorption=self.absorption,
                grids=(self.r, self.z, self.t, self.R, self.T),
                derived=self.derived,
            )

            Y_q = surrogate_yield_from_exit_field(E_q_rt, self.r, self.t)
            U_q, N_q = photons_per_pulse_from_surrogate_yield(
                Y_q=Y_q,
                E_pulse_J=self.derived["E_pulse_J"],
                q=q,
                wavelength_m=self.laser.wavelength_m,
                cal=self.calibration,
            )

            photons.append(N_q)
            energies_J.append(U_q)
            energies_eV.append(photon_energy_eV_from_q(q, self.laser.wavelength_m))
            surrogate_yields.append(Y_q)

            if q == self.harmonic.selected_q_plot:
                selected_exit_field = E_q_rt
                selected_radial_profile = np.trapz(np.abs(E_q_rt) ** 2, self.t, axis=1)

        z_focus = self.gas.jet_center_m
        I_focus = gaussian_intensity_rt(
            self.R,
            self.T,
            z_focus,
            self.derived["I0_Wm2"],
            self.laser.waist_m,
            self.derived["zR_m"],
            self.laser.tau_fwhm_s,
        )
        eta_focus = ionization_fraction_surrogate(I_focus, self.ion)

        self.results = {
            "q": np.array(q_values),
            "photons_per_pulse": np.array(photons),
            "pulse_energy_J": np.array(energies_J),
            "photon_energy_eV": np.array(energies_eV),
            "surrogate_yield": np.array(surrogate_yields),
            "selected_q": self.harmonic.selected_q_plot,
            "selected_exit_field": selected_exit_field,
            "selected_radial_profile": selected_radial_profile,
            "I_focus_Wm2": I_focus,
            "eta_focus": eta_focus,
        }
        return self.results

    def print_summary(self):
        d = self.derived
        print("=== Laser ===")
        print(f"wavelength          : {self.laser.wavelength_m*1e9:.1f} nm")
        print(f"waist               : {self.laser.waist_m*1e6:.2f} um")
        print(f"pulse FWHM          : {self.laser.tau_fwhm_s*1e15:.2f} fs")
        print(f"average power       : {self.laser.avg_power_W:.3f} W")
        print(f"rep rate            : {self.laser.rep_rate_Hz:.3e} Hz")
        print(f"pulse energy        : {d['E_pulse_J']:.3e} J")
        print(f"peak power          : {d['P_peak_W']:.3e} W")
        print(f"peak intensity      : {d['I0_Wm2']:.3e} W/m^2")
        print(f"peak intensity      : {d['I0_Wcm2']:.3e} W/cm^2")
        print(f"Rayleigh length     : {d['zR_m']*1e3:.3f} mm")
        print()
        print("=== Gas / HHG ===")
        print(f"gas species         : {self.gas.species}")
        print(f"pressure            : {self.gas.pressure_mbar:.1f} mbar")
        print(f"temperature         : {self.gas.temperature_K:.1f} K")
        print(f"peak density        : {d['N0_peak_m3']:.3e} m^-3")
        print(f"estimated cutoff    : q ~ {d['q_cut']:.1f}")
        print(f"cutoff energy       : {d['Ecut_eV']:.1f} eV")
        print()
        if hasattr(self, "results"):
            print("=== Harmonics ===")
            for qq, ee, nn in zip(
                self.results["q"],
                self.results["photon_energy_eV"],
                self.results["photons_per_pulse"],
            ):
                print(f"q={qq:2d}   E={ee:6.2f} eV   photons/pulse={nn:.3e}")

    def plot(self):
        if not hasattr(self, "results"):
            raise RuntimeError("Run simulate() before plot().")

        q = self.results["q"]
        Nph = self.results["photons_per_pulse"]
        Uq = self.results["pulse_energy_J"]
        EeV = self.results["photon_energy_eV"]

        I_focus = self.results["I_focus_Wm2"]
        eta_focus = self.results["eta_focus"]
        E_q_rt = self.results["selected_exit_field"]
        Iq_r = self.results["selected_radial_profile"]
        q_sel = self.results["selected_q"]

        fig, axes = plt.subplots(2, 2, figsize=(12, 9))

        axes[0, 0].semilogy(q, np.maximum(Nph, 1e-30), "o-", lw=2)
        axes[0, 0].set_xlabel("Harmonic order")
        axes[0, 0].set_ylabel("Photons / pulse")
        axes[0, 0].set_title("Estimated HHG photon yield")
        axes[0, 0].grid(alpha=0.3)

        axes[0, 1].semilogy(EeV, np.maximum(Uq, 1e-30), "s-", lw=2)
        axes[0, 1].set_xlabel("Photon energy [eV]")
        axes[0, 1].set_ylabel("Pulse energy per harmonic [J]")
        axes[0, 1].set_title("Estimated harmonic pulse energy")
        axes[0, 1].grid(alpha=0.3)

        im0 = axes[1, 0].imshow(
            (I_focus / 1e4).T,
            aspect="auto",
            origin="lower",
            extent=[self.r[0]*1e6, self.r[-1]*1e6, self.t[0]*1e15, self.t[-1]*1e15],
            cmap="inferno",
        )
        axes[1, 0].set_xlabel("r [um]")
        axes[1, 0].set_ylabel("t [fs]")
        axes[1, 0].set_title("Driving intensity at jet center [W/cm^2]")
        plt.colorbar(im0, ax=axes[1, 0])

        im1 = axes[1, 1].imshow(
            eta_focus.T,
            aspect="auto",
            origin="lower",
            extent=[self.r[0]*1e6, self.r[-1]*1e6, self.t[0]*1e15, self.t[-1]*1e15],
            cmap="viridis",
        )
        axes[1, 1].set_xlabel("r [um]")
        axes[1, 1].set_ylabel("t [fs]")
        axes[1, 1].set_title("Effective ionization fraction")
        plt.colorbar(im1, ax=axes[1, 1])

        plt.tight_layout()

        if E_q_rt is not None:
            plt.figure(figsize=(11, 4.5))

            plt.subplot(1, 2, 1)
            plt.imshow(
                np.abs(E_q_rt).T**2,
                aspect="auto",
                origin="lower",
                extent=[self.r[0]*1e6, self.r[-1]*1e6, self.t[0]*1e15, self.t[-1]*1e15],
                cmap="magma",
            )
            plt.xlabel("r [um]")
            plt.ylabel("t [fs]")
            plt.title(f"Exit-plane |E_q(r,t)|^2, q={q_sel}")
            plt.colorbar(label="arb. units")

            plt.subplot(1, 2, 2)
            if np.max(Iq_r) > 0:
                plt.plot(self.r * 1e6, Iq_r / np.max(Iq_r), lw=2)
            else:
                plt.plot(self.r * 1e6, Iq_r, lw=2)
            plt.xlabel("r [um]")
            plt.ylabel("Normalized near-field")
            plt.title(f"Near-field radial profile, q={q_sel}")
            plt.grid(alpha=0.3)

            plt.tight_layout()

        plt.show()


# ============================================================
# run
# ============================================================
if __name__ == "__main__":
    laser = LaserConfig(
        wavelength_m=795e-9,
        waist_m=20e-6,
        tau_fwhm_s=30e-15,
        avg_power_W=8.0,
        rep_rate_Hz=1e3,
    )

    gas = GasConfig(
        species="Ar",
        pressure_mbar=50.0,
        temperature_K=300.0,
        medium_length_m=1.2e-3,
        jet_center_m=0.0,
        profile="gaussian",
        sigma_z_m=0.35e-3,
    )

    harmonic = HarmonicConfig(
        q_min=11,
        q_max=51,
        selected_q_plot=21,
        Ip_eV=15.7596,
    )

    grid = GridConfig(
        Nr=140,
        Nz=260,
        Nt=220,
        r_max_m=120e-6,
        t_max_s=90e-15,
    )

    ion = IonizationConfig(
        model="surrogate",
        I_sat_Wm2=1.2e18,
        gamma=6.0,
        eta_max=0.98,
    )

    dip = DipoleConfig(
        short_weight=1.0,
        long_weight=0.35,
        alpha_short_m2_per_W_at_q21=0.8e-18,
        alpha_long_m2_per_W_at_q21=3.0e-18,
        p_base=4.0,
        p_slope_per_harmonic=0.015,
        ion_suppression_beta=10.0,
        rolloff_width_harmonics=5.0,
        long_amplitude_factor=0.6,
    )

    absorption = AbsorptionConfig(
        model="powerlaw",
        L_abs_ref_m_at_q21=0.6e-3,
        power_law_exponent=0.5,
    )

    calibration = CalibrationConfig(
        scale=2.0e-16
    )

    model = HHGArgonModel(
        laser=laser,
        gas=gas,
        harmonic=harmonic,
        grid=grid,
        ion=ion,
        dip=dip,
        absorption=absorption,
        calibration=calibration,
    )

    model.simulate()
    model.print_summary()
    model.plot()