

                                READ ME


All shock_rad directories with a suffix "latest" contain the following
things:

1. Semi-analytical expression of the photon escape timescale,
calculated using the "esctime.h" self-made library program. The
program uses the publicly available but slightly modified version of
the Monte-Carlo numerical integration technique program "vegas_arh.h".

2. The self-made library programs, syn_related.h, ssc_related.h, and
pair_prod.h that have been modified to incorporate the trapezoidal
rule for carrying out numerical integration to calculate quantities
like synchrotron/SSC photon density rate, etc.

3. Lot of local variables have been renamed in these library programs.

4. All pow(x, 2.0) have been changed to SQR(x).

5. All multiplicative constants have been replaced by their actual
values wherever possible.

6. Instantaneous photon density is added to the leftover total photon
density from previous time step. The new total photon density is used
to calculate the escape photon density rates. The escape rates are
then subtracted from the total photon density and this updated photon
density is what's used in the ssc density, ssc losses, and pair
production calculations.
 
7. The distance travelled by the photons going in the forward
direction has been multiplied with the cos_theta term to make sure
that it falls along the line of sight of the observer.

8. The case of viewing the jet at oblique angles where the back side
of the jet would be visible to the observer first has also been
incorporated at the end of the shock_rad.c program.

9. Some parameters which are calculated only once have been put
outside of a loop. All files with "constparams" include that change.

10. The post processing file has been corrected. The dividing factor
is now the binsize and not the number of files.

11. All the constant parameters have been put before the while loop.

12. The value of Y has been changed to ro_i/ro_o.

13. The value of ro_i & ro_o has been changed to their comoving
values.

14. The parameter zeta_e has been replaced with epsilelec in the Qinj
calculation.

15. Gamma_fs(rs) has been replaced with (Gamma_fs(rs) - 1) when
calculating gamma_b_fs(rs) to take into account the contribution of
the protons' kinetic energy instead of their total internal energy
being transferred to relativistic electrons.

16. The condition that is gamma_min > gamma_max, then the program
should stop has been included.

17. The total internal energy of a shocked region, in the lab frame,
E_fs(rs) has been replaced with (E_fs(rs) * (1 - (1/Gamma_sh))) to
obtain the kinetic energy of the shocked fluid in the comoving frame
of the shocked shell, and use that in B & Qinj calculations.

18. Changes made in the ppc_binsize.c program. The extra while and for
loop was removed that was being used to read a file.

19. In synchrotron related .h program, the initialization of yofx[j]
to zero has been removed, to save some time.

20. The nep[j] gets zero values for which the u_i[j] has values
smaller than 1.0e-4 and if they don't lie between gamma_min and
gamma_max.  - This condition has been removed and original one is
being used.

21. The calculations in ssc gammadot and pair production sigmagg
functions are skipped with the "continue" statments for which nph[j]
values are smaller than 1.0e-15, and nep[j] values are equal to zero
for the entire subroutine.

22. The calculation of emission region distance from the central
region has been put in the code. The distance is from the position of
the CD and would be the same for both FS and RS region in the lab
frame.

23. Added extra sets of parantheses in SQR, pow, exp, and sqrt
functions (where applicable), in all .h (except for myrtnew.h,
sortfile.h, & tridag.h) and .c (except for time_int_pp.c) programs.

24. Added Eint, Gamma_m, and eta_IS calculation in shockrad_synssc.c

25. Fixed accf == FALSE to accr == FALSE in the 2nd part of the post
processing loop.
