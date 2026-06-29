from __future__ import print_function
import math
import numpy as np


class LAMMPSPairPotential(object):
    """Base class for LAMMPS pair potential models.

    Provides common functionality for mapping particle types and checking unit consistency.
    """

    def __init__(self):
        self.pmap = dict()
        self.units = "lj"

    def map_coeff(self, name, ltype):
        self.pmap[ltype] = name

    def check_units(self, units):
        if units != self.units:
            raise Exception("Conflicting units: %s vs. %s" % (self.units, units))


class Harmonic(LAMMPSPairPotential):
    """Implements the Harmonic pair potential model for particle interactions.

    This class defines the Harmonic potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(Harmonic, self).__init__()
        self.units = "real"
        # set coeffs: K, r0
        self.coeff = {
            "A": {"A": (0.2, 9.0), "B": (math.sqrt(0.2 * 0.4), 9.0)},
            "B": {"A": (math.sqrt(0.2 * 0.4), 9.0), "B": (0.4, 9.0)},
        }

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1] - r
        if r <= coeff[1]:
            return 2.0 * delta * coeff[0] / r
        else:
            return 0.0

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r = math.sqrt(rsq)
        delta = coeff[1] - r
        if r <= coeff[1]:
            return delta * delta * coeff[0]
        else:
            return 0.0


class LJCutMelt(LAMMPSPairPotential):
    """Implements the Lennard-Jones Cutoff (Melt) pair potential model for particle interactions.

    This class defines the LJCutMelt potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(LJCutMelt, self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = "lj"
        self.coeff = {"lj": {"lj": (48.0, 24.0, 4.0, 4.0)}}

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1 * r6inv - lj2)) * r2inv

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return r6inv * (lj3 * r6inv - lj4)


class LJCutSPCE(LAMMPSPairPotential):
    """Implements the Lennard-Jones Cutoff (SPCE) pair potential model for particle interactions.

    This class defines the LJCutSPCE potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(LJCutSPCE, self).__init__()
        self.units = "real"
        # SPCE oxygen LJ parameters in real units
        eps = 0.15535
        sig = 3.166
        self.coeff = {
            "OW": {
                "OW": (
                    48.0 * eps * sig**12,
                    24.0 * eps * sig**6,
                    4.0 * eps * sig**12,
                    4.0 * eps * sig**6,
                ),
                "HW": (0.0, 0.0, 0.0, 0.0),
            },
            "HW": {"OW": (0.0, 0.0, 0.0, 0.0), "HW": (0.0, 0.0, 0.0, 0.0)},
        }

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        lj1 = coeff[0]
        lj2 = coeff[1]
        return (r6inv * (lj1 * r6inv - lj2)) * r2inv

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r6inv = r2inv * r2inv * r2inv
        lj3 = coeff[2]
        lj4 = coeff[3]
        return r6inv * (lj3 * r6inv - lj4)


class Patchy(LAMMPSPairPotential):
    """Implements the Patchy pair potential model for particle interactions.

    [1] I. Palaia and A. Šarić, “Controlling cluster size in 2D phase-separating binary mixtures with specific interactions,” J. Chem. Phys., vol. 156, no. 19, p. 194902, May 2022, doi: 10.1063/5.0087769.

    This class defines the Patchy potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(Patchy, self).__init__()
        self.units = "real"
        siglj = 12.631578947368421
        Rclj = 14.178467978644711
        epslj = 2.0
        epsp = 0.4
        sigp = 0.1 * Rclj
        Rcp = 0.3 * Rclj
        # Patchy LJ parameters in real units
        # set coeffs: K, r0
        self.coeff = {
            "Ac": {
                "Ac": (siglj, Rclj, 0, epslj, epsp),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bc": (siglj, Rclj, 2, epslj, epsp),
                "Bp": (0.0, 0.0, 1, 0, 0),
            },
            "Ap": {
                "Ac": (0.0, 0.0, 1, 0, 0),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bc": (0.0, 0.0, 1, 0, 0),
                "Bp": (sigp, Rcp, 3, epslj, epslj),
            },
            "Bc": {
                "Ac": (siglj, Rclj, 2, epslj, epsp),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bp": (0.0, 0.0, 1, 0, 0),
                "Bc": (siglj, Rclj, 4, epslj, epsp),
            },
            "Bp": {
                "Ac": (0.0, 0.0, 1, 0, 0),
                "Ap": (sigp, Rcp, 3, epslj, epslj),
                "Bc": (0.0, 0.0, 1, 0, 0),
                "Bp": (0.0, 0.0, 1, 0, 0),
            },
        }

    def compute_force(self, rsq, itype, jtype):
        """Calculates the force between two particles using the Patchy model.

        Computes the force based on the squared distance and particle types, applying the Patchy potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed force.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        sigma = coeff[0]
        Rc = coeff[1]
        typ = coeff[2]
        eps = coeff[3]
        epsp = coeff[4]
        force = 0.0
        match typ:
            case 0:
                force = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                return force
            case 1:
                return 0.0
            case 2:
                force = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                return force
            case 3:
                r = math.sqrt(rsq)
                force = 0.0
                if r <= Rc and r > sigma:
                    force = (
                        -2
                        * (
                            eps
                            * np.cos(math.pi * (r - sigma) / (Rc - sigma) / 2)
                            * np.sin(math.pi * (r - sigma) / (Rc - sigma) / 2)
                            * (math.pi)
                            / (2 * (Rc - sigma))
                        )
                        / r
                    )
                return force
            case 4:
                force = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                return force
            case _:
                return 0.0

    def compute_energy(self, rsq, itype, jtype):
        """Calculates the potential energy between two particles using the Patchy model.

        Computes the energy based on the squared distance and particle types, applying the Patchy potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed potential energy.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        sigma = coeff[0]
        Rc = coeff[1]
        typ = coeff[2]
        eps = coeff[3]
        epsp = coeff[4]
        match typ:
            case 0:
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                ener = 0.0
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2
                    ) + ener
                return ener
            case 1:
                return 0.0
            case 2:
                ener = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2
                    ) + ener
                return ener
            case 3:
                ener = 0.0
                r = math.sqrt(rsq)
                if r < Rc:
                    if r <= sigma:
                        ener = -eps
                    else:
                        ener = (
                            -eps * np.cos(math.pi * (r - sigma) / (Rc - sigma) / 2) ** 2
                        )
                return ener
            case 4:
                ener = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -(epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2)
                        + ener
                    )
                return ener
            case _:
                return 0.0


class PatchyYuk(LAMMPSPairPotential):
    """Implements the Patchy pair potential model for particle interactions with a Yukawa-like potential added to B-B interactions.

    [1] I. Palaia and A. Šarić, “Controlling cluster size in 2D phase-separating binary mixtures with specific interactions,” J. Chem. Phys., vol. 156, no. 19, p. 194902, May 2022, doi: 10.1063/5.0087769.

    This class defines the Patchy potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(PatchyYuk, self).__init__()
        self.units = "real"
        siglj = 12.631578947368421
        Rclj = 14.178467978644711
        epslj = 2.0
        epsp = 0.4
        sigp = 0.1 * Rclj
        Rcp = 0.3 * Rclj
        kappa = 0.05 / siglj
        # Patchy LJ parameters in real units
        # set coeffs: K, r0
        self.coeff = {
            "Ac": {
                "Ac": (siglj, Rclj, 0, epslj, epsp),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bc": (siglj, Rclj, 2, epslj, epsp),
                "Bp": (0.0, 0.0, 1, 0, 0),
            },
            "Ap": {
                "Ac": (0.0, 0.0, 1, 0, 0),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bc": (0.0, 0.0, 1, 0, 0),
                "Bp": (sigp, Rcp, 3, epslj, epslj),
            },
            "Bc": {
                "Ac": (siglj, Rclj, 2, epslj, epsp),
                "Ap": (0.0, 0.0, 1, 0, 0),
                "Bp": (0.0, 0.0, 1, 0, 0),
                "Bc": (siglj, Rclj, 4, epslj, epsp, kappa),
            },
            "Bp": {
                "Ac": (0.0, 0.0, 1, 0, 0),
                "Ap": (sigp, Rcp, 3, epslj, epslj),
                "Bc": (0.0, 0.0, 1, 0, 0),
                "Bp": (0.0, 0.0, 1, 0, 0),
            },
        }

    def compute_force(self, rsq, itype, jtype):
        """Calculates the force between two particles using the Patchy model.

        Computes the force based on the squared distance and particle types, applying the Patchy potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed force.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        sigma = coeff[0]
        Rc = coeff[1]
        typ = coeff[2]
        eps = coeff[3]
        epsp = coeff[4]
        force = 0.0
        match typ:
            case 0:
                force = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                return force
            case 1:
                return 0.0
            case 2:
                force = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                return force
            case 3:
                r = math.sqrt(rsq)
                force = 0.0
                if r <= Rc and r > sigma:
                    force = (
                        -2
                        * (
                            eps
                            * np.cos(math.pi * (r - sigma) / (Rc - sigma) / 2)
                            * np.sin(math.pi * (r - sigma) / (Rc - sigma) / 2)
                            * (math.pi)
                            / (2 * (Rc - sigma))
                        )
                        / r
                    )
                return force
            case 4:
                force = 0.0
                kappa = coeff[5]
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                dyuk = 3*eps * math.exp(-kappa * r) * (kappa + 1 / r) / r**2
                if r <= 2 * Rc:
                    if r <= Rc:
                        force = 12 * eps * (r6inv * r6inv - r6inv) * r2inv
                    else:
                        force = -2 * (
                            epsp
                            * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * np.sin(math.pi * (r - Rc) / (2 * Rc - Rc) / 2)
                            * (math.pi)
                            / (2 * (2 * Rc - Rc))
                            / r
                        )
                if r < 6 * Rc:
                    force += dyuk
                return force
            case _:
                return 0.0

    def compute_energy(self, rsq, itype, jtype):
        """Calculates the potential energy between two particles using the Patchy model.

        Computes the energy based on the squared distance and particle types, applying the Patchy potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed potential energy.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        sigma = coeff[0]
        Rc = coeff[1]
        typ = coeff[2]
        eps = coeff[3]
        epsp = coeff[4]
        match typ:
            case 0:
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                ener = 0.0
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2
                    ) + ener
                return ener
            case 1:
                return 0.0
            case 2:
                ener = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2
                    ) + ener
                return ener
            case 3:
                ener = 0.0
                r = math.sqrt(rsq)
                if r < Rc:
                    if r <= sigma:
                        ener = -eps
                    else:
                        ener = (
                            -eps * np.cos(math.pi * (r - sigma) / (Rc - sigma) / 2) ** 2
                        )
                return ener
            case 4:
                kappa = coeff[5]
                ener = 0.0
                r = math.sqrt(rsq)
                r2inv = 1.0 / rsq
                r6inv = Rc**6 * r2inv * r2inv * r2inv
                Rcy = 6*Rc
                yukrc = eps * math.exp(-kappa * Rcy) / Rcy
                yuk = 3*(eps * math.exp(-kappa * r) / r - yukrc)
                if r <= Rc:
                    ener = -epsp + eps * (r6inv * r6inv - 2 * r6inv + 1)
                elif r <= 2 * Rc:
                    ener = (
                        -(epsp * np.cos(math.pi * (r - Rc) / (2 * Rc - Rc) / 2) ** 2)
                        + ener
                    )
                if r < Rcy:
                    ener += yuk
                return ener
            case _:
                return 0.0


class SALR(LAMMPSPairPotential):
    def __init__(self):
        """Initializes the SALR (Short-range Attractive Long-range Repulsive) potential parameters.

        [1]A. Imperio and L. Reatto, “Microphase morphology in two-dimensional fluids under lateral confinement,” Phys. Rev. E, vol. 76, p. 40402, 2007, doi: 10.1103/PhysRevE.76.040402.


            Sets up the coefficients and units for the SALR pair potential model.
        """
        super(SALR, self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = "real"
        eps = 0.1
        sig = 4.0
        sig2 = sig * sig
        sig10 = sig2 * sig2 * sig2 * sig2 * sig2
        k1 = 1.0
        a1 = 0.1 / sig
        k2 = 2.0
        a2 = 0.25 / sig
        rc = 60.0
        r2inv = 1 / rc**2
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        urc = sig10 * r10inv + k1 * np.exp(-a1 * rc) - k2 * np.exp(-a2 * rc)
        self.coeff = {
            "1": {
                "1": (eps, sig10, k1, k2, a1, a2, urc),
                "2": (0.5 * eps, sig10, k1, k2, a1, a2, urc),
            },
            "2": {
                "1": (0.5 * eps, sig10, k1, k2, a1, a2, urc),
                "2": (eps, sig10, k1, k2, a1, a2, urc),
            },
        }

    def compute_force(self, rsq, itype, jtype):
        """Calculates the force between two particles using the SALR model.

        Computes the force based on the squared distance and particle types, applying the SALR potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed force.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r = np.sqrt(rsq)
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        duslr = (
            10 * sig10 * r10inv * r2inv
            + (k1 * a1 * np.exp(-a1 * r) - k2 * a2 * np.exp(-a2 * r)) / r
        )
        return eps * duslr
        # r6inv  = r2inv*r2inv*r2inv

        # return (r6inv * (lj1*r6inv - lj2))*r2inv

    def compute_energy(self, rsq, itype, jtype):
        """Calculates the potential energy between two particles using the SALR model.

        Computes the energy based on the distance squared and particle types, applying the SALR potential parameters.

        Args:
            rsq: Squared distance between two particles.
            itype: Type of the first particle.
            jtype: Type of the second particle.

        Returns:
            float: The computed potential energy.
        """
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r = np.sqrt(rsq)
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        urc = coeff[6]
        uslr = sig10 * r10inv + k1 * np.exp(-a1 * r) - k2 * np.exp(-a2 * r) - urc
        return eps * uslr


class OPP(LAMMPSPairPotential):
    """Implements the OPP (Oscillatory Pair Potential) model for particle interactions.

    [1]M. Mihalkovič and C. L. Henley, “Empirical oscillating potentials for alloys fromab initiofits and the prediction of quasicrystal-related structures in the Al-Cu-Sc system,” Phys. Rev. B, vol. 85, no. 9, p. 92102, Mar. 2012, doi: 10.1103/physrevb.85.092102.


    This class defines the OPP potential, including initialization of parameters and methods to compute force and energy.
    """

    def __init__(self):
        super(SALR, self).__init__()
        # set coeffs: 48*eps*sig**12, 24*eps*sig**6,
        #              4*eps*sig**12,  4*eps*sig**6
        self.units = "real"
        eps = 0.1
        sig = 4.0
        sig2 = sig * sig
        sig10 = sig2 * sig2 * sig2 * sig2 * sig2
        k1 = 1.0
        a1 = 0.1 / sig
        k2 = 2.0
        a2 = 0.25 / sig
        rc = 60.0
        # Parameters for I4 phase (Glotzer, Nature Materials)
        k = 8.5 / sig
        phi = 0.47
        rmx = 2.769 * sig
        sigma = sig
        r2inv = 1 / rc**2
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        urc = sig10 * r10inv + k1 * np.exp(-a1 * rc) - k2 * np.exp(-a2 * rc)
        self.coeff = {"1": {"1": (eps, sig10, k1, k2, a1, a2, urc, k, phi, rmx, sigma)}}

    def compute_force(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r = np.sqrt(rsq)
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        k = coeff[7]
        phi = coeff[8]
        rmx = coeff[9]
        sigma = coeff[10]
        sig3 = sigma**3
        if r < rmx:
            duslr = (
                10 * sig10 * r10inv * r2inv
                + (
                    k1 * a1 * np.exp(-a1 * r)
                    - k2 * a2 * np.exp(-a2 * r)
                    + 3 * sig3 * np.cos(k * (r - 1.25 * sigma) - phi) * r2inv**2
                    + sig3 * k * np.sin(k * (r - 1.25 * sigma) - phi) / r**3
                )
                / r
            )
        else:
            duslr = (
                10 * sig10 * r10inv * r2inv
                + (k1 * a1 * np.exp(-a1 * r) - k2 * a2 * np.exp(-a2 * r)) / r
            )
        return eps * duslr

    def compute_energy(self, rsq, itype, jtype):
        coeff = self.coeff[self.pmap[itype]][self.pmap[jtype]]
        r2inv = 1.0 / rsq
        r = np.sqrt(rsq)
        r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
        eps = coeff[0]
        sig10 = coeff[1]
        k1 = coeff[2]
        k2 = coeff[3]
        a1 = coeff[4]
        a2 = coeff[5]
        urc = coeff[6]
        k = coeff[7]
        phi = coeff[8]
        rmx = coeff[9]
        sigma = coeff[10]
        sig3 = sigma**3
        if r > rmx:
            uslr = sig10 * r10inv + k1 * np.exp(-a1 * r) - k2 * np.exp(-a2 * r) - urc
        else:
            uslr = (
                sig10 * r10inv
                + k1 * np.exp(-a1 * r)
                - k2 * np.exp(-a2 * r)
                - urc
                + sig3 * np.cos(k * (r - 1.25 * sigma) - phi) / r**3
                - sig3 * np.cos(k * (rmx - 1.25 * sigma) - phi) / rmx**3
            )
        return eps * uslr
