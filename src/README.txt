https://grinch.phys.uconn.edu:2880/Gluex/beamline/phonopolar/4-2017/README.txt

   Monte Carlo dataset for Calculation of Photon Beam Properties for GlueX
                       Richard.T.Jones at uconn.edu
                        prepared on June 18, 2017

*** Contents of this document ***
      I.  Purpose
     II.  Monte Carlo model
    III.  Beamline optics model
     IV.  Directories structure
      V.  ROOT trees dictionary
     VI.  Example analysis
    VII.  Advanced techniques


   I.  Purpose

Systematic uncertainty in the flux and polarization spectra for a coherent
bremsstrahlung photon beam arises from two distinct sources: the physics
model for the coherent bremsstrahlung process, and the optics model for the
crystal and beamline. Of the two, the crystal/optics model is the more
challenging one to constrain because the beamline is still new and not very
well understood. The limited data we have suggest that the optics are not
constant over the course of a run. On the other hand, the basic physics of
coherent bremsstrahlung is quite well understood, with differential cross
sections based on leading-order QED expected to be accurate at the percent
level outside the channeling regime with Egamma > 0.1 Eelectron.

Folding the differential cross section with the crystal/optics/collimator
model to predict spectra involves integration over many degrees of freedom.
In the past, experimenters have relied on formulas in which these integrals
have been carried out analytically under certain approximations. The purpose
of this dataset is to enable studies of the dependence of photon beam spectra
on various crystal/beamline optics properties that are free of such
uncontrolled approximations by doing the integrals using Monte Carlo
integration. A high-statistics sample of CB production events is provided,
broken down into separate datasets for coherent (zero-phonon radiation),
"incoherent nuclear" (radiation involving one or more phonons in the
crystal), and "incoherent electron" (radiation with the recoil momentum
carried away by a single recoil electron) processes. Summing over the
events from all three datasets, weighted by the product of the differential
cross section and the Monte Carlo importance-sampling weight for each event,
generates predicted spectra for an arbitrary set of cuts that
represent the effects of beam optics and collimation.


  II.  Monte Carlo model

The events thrown by the Monte Carlo generator contain the full kinematics
of the beam electron and the radiated photon. The energy/momentum/mass of
the recoil system can be inferred by momentum conservation. In the 4-2017
datasets, the incident electron momentum is along the z axis with zero
momentum spread and divergence angle, with energy 12.0 GeV. The crystal
is assumed to be ideal diamond at 300 K, aligned in the PERP orientation.
As for the larger crystal rocking angle, a set of 15 different choices
for this angle: 0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80,
0.90, 0.10, 0.13, 0.16, 0.19, and 0.25 are included in this dataset. For
each of these crystal orientations, 3 samples were generated for the
3 principal physics processes: coherent radiation, incoherent nuclear
radiation, and incoherent electron radiation. In each case, a weighted
phase space generator was used to throw the event kinematics, and the
sum of all QED diagrams at tree level was used to compute the differential
cross section for the given process. Details for each process type are
given below. The differential cross section units are different for the
three cases, so attention should be paid to how each is normalized. In
each case, the total flux is obtained by histogramming the respective
event sample subject to collimation cuts, weighted by the Monte Carlo
importance weight factor multiplied by the product differential cross
section times target thickness times beam current, divided by the total
number of Monte Carlo events generated.


  II.A) coherent radiation

This process is more kinematically restrictive than the other two. The
energy of the radiated photon is a multi-valued function of the emission
angle with respect to the incident electron direction (z-axis). Each
branch of the function corresponds to one of the allowed reciprocal lattice
vectors of the diamond crystal. CB events are generated at all angles with
respect to the beam axis, so some finite cut on the photon direction is
needed to make a sensible coherent beam spectrum. The sum over reciprocal
lattice vectors that was included in the generator extends well above 100
in terms of the transverse h and k indices, and above 20 in the l index.
The difference is because the l index is cut off much more sharply by the
cross section than are the transverse indices h,k. The polarization
carried by each individual photon is stored together with the 4-momenta
in the data files. The one part of this calculation that is not uniquely
specified by QED is the atomic form factor of the carbon atoms. For this,
I used the dipole form factor in the parameterized form given by Eq. 98
of Ref. U. Timm, Fortschritte der Physik 17, pp. 765-808, 1969. To improve
upon this approximate form factor, simply multiply the recorded
differential cross section for each event in the data files by the ratio
improved / dipole form factors squared evaluated at the q2 for that event.


  II.B) incoherent nuclear radiation

This process is really an identical calculation to the coherent one
covered in the previous section, except that this one covers the case
where the coherent radiation is accompanied by the emission/absorption
of one or more phonons on the crystal lattice of the radiator. The
reason they are separated here is because of singularities that appear
in the integrand for the case of zero phonons, which bcome smooth functions
of the momentum transfer for final states with a different number of
phonons than the initial state. The incoherent atomic (nuclear) differential
cross section is computed as a sum over reciprocal lattice vectors the same
as the coherent part. Screening of the nuclear charge at low q2 by the
atomic electrons is accounted for by an atomic form factor. The same
approximate dipole atomic form factor was used for this part as was used
for the coherent part. The method for replacing the dipole formula with
a more accurate form factor that was described in the previous section
also applies here.

The incoherent cross section is somewhat dependent on the temperature of
the crystal, and on the phonon dispersion curves that are assumed for
the crystal. The model for phonon dispersion that is used here is a simple
symmetric 3D Debye model for the optical phonon branch, with a speed of
sound vP = 15 km/s for longitudinal (pressure P-wave) modes and vS = 
12 km/S for transverse (shear S-wave) modes. A simple Einstein model is
used for the optical phonon branch, with a zero-dispersion phonon frequency
adjusted to match the Debye temperature of diamond of 2200K. Any of these
parameters can be readily varied without regenerating the dataset. How
this is done is presented below under Advanced Techniques.


  II.C) incoherent electron radiation

Final states where the momentum transfer to the crystal is taken up by
a single recoiling electron are described by an elementary radiative
Moller scattering QED process. Because there are 2 electrons involved,
I decided to compute the cross section and polarization for this process
in a completely separate calculation, in which the sum is taken over all
8 tree-level diagrams. This process contributes about 15% of the total
photon yield at the energies of interest to GlueX, and has a completely
different polarization signature from the other two processes. There is
no model dependence in this calculation, beyond the fact that higher-order
QED diagrams, the ones with loops, can be neglected.


 III.  Beamline optics model

No smearing of the electron beam, either from beamline optics or from
multiple-scattering in the radiator, is included in the generator. The
event generator assumes the incident electron enters the reaction with
a momentum of 12.000 GeV directed precisely along the z axis. Effects 
coming from beam momentum spread, multiple scattering in the radiator,
and effects from non-ideal crystal structure in the radiator are all
taken into account when the data are analyzed and spectra are generated.
Many of the properties of the photon beam spectra, both polarization and
rate, are smooth functions of the kinematics whose variation on the scale
of the momentum spread and emittance of the CEBAF 12 GeV electron beam
are entirely negligible, with two outstanding exceptions. The first
exception is the effects they have in smearing out the collimation function,
transforming a rectangular cutoff on the low-energy side of the coherent
peaks into a smoothed out tail. The second exception is in the effects 
they have in smearing the coherent edges, transforming them from an abrupt
discontinuity (as seen in the simulation) to the smooth step seen in
experimental plots of the intensity or polarization vs. energy.

A straight-forward histogramming of all of the events in the generated
samples shows no truncation on the low-energy side of the coherent peaks;
they extend from the respective coherent edge all the way to zero energy
because no collimation is assumed by the generator. This gives full
freedom in the analysis to assume an arbitrary collimation function.
The simplest approach, and the one that is used in the example ROOT
TSelector described below, is to assume a simple Gaussian model for the
phase space of the electron beam at its virtual focus, and generate an
independent ray to represent the electron for each generated event. 
This provides a simple geometric transformation (translation + rotation)
to be applied to the kinematics of that event, leading to a final photon
ray for that event which can be projected downstream to the collimator
face and either accepted or rejected, based on the dimensions of the
collimator aperture. From this simple starting point, non-Gaussian
tails on the beam phase-space distribution, etc., can be introduced
for study. The exact same algoritm should be applied to all three
generated datasets (coherent nuclear, incoherent nuclear, and electron
recoil) to obtain consistent results.

Most of the work that will be needed to understand the properties of the
GlueX photon beam will be focused on getting the effects of collimation
to agree between simulation and experiment. The remaining task will be
to get the shapes of the coherent edges to agree. In GlueX, this comes
mostly from the mosaic spread of the radiator crystal, with a small
additional contributionf from the spread in the endpoint energy of the
electron beam and its temporal variation over the duration of the
measurement period. In the simulation, this endpoint energy smearing is
easily applied at the analysis step by multiplying all of the energies 
of the particles in a particular event by a Gaussian factor with mean of
1 and a sigma of a few MeV. The same technique can also be used to 
account for crystal mosaic spread as long is it is small, amounting to
less than 50 MeV of smearing in the photon energy. If the need arises to
explore larger deviations of the crystal mosaic spread corresponding to
coherent edge smearing of hundreds of MeV, new simulations will need to
be performed in which the crystal rocking curves that have been measured
at a synchrotron X-ray source are explicitly taken into account.


  IV.  Directories structure

The data from these three datasets are shared from a data server that
supports a variety of remote access protocols, including ordinary http,
gridftp, and xrootd. All access to these data requires that you have a
valid Gluex grid certificate, and that it is loaded into your working
environment or your browser. The easiest way to explore the datasets is
by directing a web browser to the following address.

https://grinch.phys.uconn.edu:2280/Gluex/beamline/phonopolar/4-2017

If you get an access error, verify that your grid certificate is still
valid, and that it is loaded into your browser's certificate store as
a client certificate. The listing under the URL shown above should look
like the following.

     eebrems.d		Sun Jun 04 16:52:36 EDT 2017
     perp.000.d		Wed Jun 07 07:39:46 EDT 2017
     perp.010.d		Sun Jun 04 09:09:58 EDT 2017
     perp.020.d		Tue Jun 06 07:04:35 EDT 2017
     perp.030.d		Sun Jun 04 10:36:13 EDT 2017
     perp.040.d		Sun Jun 04 17:26:43 EDT 2017
     perp.050.d		Wed Jun 07 07:39:55 EDT 2017
     perp.060.d		Sun Jun 04 08:39:12 EDT 2017
     perp.070.d		Sun Jun 04 14:35:41 EDT 2017
     perp.080.d		Sun Jun 04 11:06:32 EDT 2017
     perp.090.d		Sun Jun 04 16:45:23 EDT 2017
     perp.100.d		Wed Jun 07 07:40:02 EDT 2017
     perp.130.d		Wed Jun 07 07:38:59 EDT 2017
     perp.160.d		Wed Jun 07 07:38:59 EDT 2017
     perp.190.d		Wed Jun 07 07:38:59 EDT 2017
     perp.250.d		Wed Jun 07 07:38:59 EDT 2017

The names of each of these directories is encoded as follows. There
is just one recoil electron process dataset, under the eebrems.d folder.
The perp.YYY.d directories each contain a full simulation of the incoherent
nuclear process, where the crystal orientation is as follows.

     theta_X = 3.331 mrad (sets coherent edge at 9.000 GeV, PERP alignment)
     theta_Y = 0.YYY mrad (YYY is a component of the folder name)

For initial studies, you may just want to focus on one typical value
for YYY, as most properties are only weakly sensitive to theta_Y in
the PERP cyrstal orientation.

Inside each of these directories listed above is an embedded directory
called phonoZero.d, which contains the data for the coherent nuclear 
process. Here, the "Zero" in the foldername is suggestive of the fact
that the data inside that folder were generated for the zero-phonon 
process.

All of the data belonging to a single dataset reside together inside a
single folder that ends with the .d suffix. The data are contained in
root trees that are chained across many files, typically 400 per dataset.
The filenames all have similar names with a serial number embedded in 
the filename. The order of the files is insignificant. However, you do
need to know one additional piece of information that is not contained
inside these files: how many generated MC events were thrown to generate
each dataset. This is not identical to the number of rows in the TChain
because events with a differential cross section * weight product below
some cutoff have been eliminated from the output when the trees were
written. The weights have been compensated to account for this loss, so
no correction from this truncation needs to be applied during analysis.
The total generated statistics for each dataset is listed below.

    * eebrems.d   :    4 x 10^9 generated
    * perp.YYY.d  :    1 x 10^11 generated
    * phonoZero.d :    1 x 10^11 generated

These numbers will be used in the example analysis below to normalize
the computed rates for a given crystal/optics model.

Other remote access protocols use a similar URL structure to access
the same files under the same hierarchical directory structure. I include
examples of each below.

root://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/perp.030.d/phonotree_23a.root
root://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/perp.250.d/phonoZero.d/phonotree_18b.root
root://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/eebrems.d/eebrems_67c.root

gsiftp://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/perp.030.d/phonotree_23a.root
gsiftp://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/perp.250.d/phonoZero.d/phonotree_18b.root
gsiftp://nod29.phys.uconn.edu/Gluex/beamline/phonopolar/4-2017/eebrems.d/eebrems_67c.root

      V.  ROOT trees dictionary

All of the root trees in the datasets described above have an identical
row structure. Some of the columns are omitted from the recoil electron
dataset trees because they do not apply to that case.

   double k[3]         : 3-momentum of radiated photon (GeV/c)
   double qphoton[3]   : 3-momentum of exchange photon (GeV/c)
   double Qphonons[3]  : total 3-momentum of phonons, if any (GeV/c)
   double diffXS       : differential cross section (microbarns)
   double polar        : linear polarization (dimensionless)
   double wgt          : Monte Carlo weight factor (dimensionless)
   double costhetastar : internal generator variable, please ignore
   double logwgts[10]  : internal generator variables, please ignore
   int hkl[3]          : Miller indices of scattering vector, if any

Normally one expects that a differential cross section would have units
in the denominator, but here the denominator has been converted to counts
for easy histogramming of results, so the units are already in microbarns.
For the coherent and incoherent nuclear datasets, the normalization is in
microbarns per unit cell of the target. For the eebrems simulation, the
normalization is in microbarns per electron in the target. Keep that
difference in mind when summing these two different kinds of reactions.


     VI.  Example analysis

As an illustration of how these data can be used to generate a predicted
spectrum for rates and polarization is available online at the following URL.

http://zeus.phys.uconn.edu/~halld/phonopolar.C.html
http://zeus.phys.uconn.edu/~halld/phonopolar.h.html

This sample TSelector can be called on a TChain object created by chaining
together all 400 files from a given dataset, or it can be invoked from a
PROOF session for remote execution on a parallel PROOF cluster. In either
case, the results appear upon return from the TCain::Process method in the
phonopolar.root output file. The meaning of the output histograms should be
clear just from reading the code and seeing the comments therein.


    VII.  Advanced techniques

Coming soon.
