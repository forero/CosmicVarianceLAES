Major points

1.  As the authors are aware, each observational survey for LAEs has different selection effects, and some use very different thresholds in equivalent width.  Hence you need to use the number density and autocorrelation function (ACF) of LAEs from the *same* survey for an investigation like this.  You could try to get Yamada's LAE positions to calculate the ACF, since you appear to prefer that survey for number density estimation.  Or you could pursue joint analysis of number density and ACF for smaller surveys e.g., Ouchi's.

- [] Check where are the number densities for the other surveys with
respect to Yamadas's 

- [X] Mail Yamada to see whether we can have access to their results for
the correlation function. 

2.  The major result claimed that all models have a narrow mass range of halos hosting LAEs is guaranteed by the choice of a minimum value of f_occ of 0.1.  This should be obvious to the authors, as there is a degeneracy between f_occ and the mass range of halos that you have hidden by forcing f_occ to values larger than those previously claimed in the literature.

- [] Clarify the point. Not all models have a narrow mass range. 

3.  The authors cite similar work by Walker Soler et al. 2012 but fail to mention in the introduction how the current paper has the same goals as that earlier work.  If applied more carefully, their approach should be able to improve on the earlier work, but it is important to be clear about how similar the papers are and how the current one hopes to improve on the previous one.

- [] Clarify this point.

4.  Comparing the most heavily clustered of 15 mocks against SSA22 is unfair.
SSA22 is a ~5 sigma overdensity, which is much rarer than 1/15.

- [] Clarify this point. It is a ~6 times more dense than the average
field *as measured for LBG galaxies*. It is only ~4 times more dense
if we measure it in LAEs. We don't know the correspondence for dark
matter haloes. That's why we make a blind analysis. 

5.  There seems to be significant confusion between halo occupation fraction and Lyman alpha escape fraction.  Only in an incredibly simplified model where 100% of Ly alpha photons escape from LAEs and 0% escape from non-LAEs would the comparison in section 4.1 make sense.  In a more realistic model, it is not clear how this investigation could hope to comment on the escape fraction of Ly alpha photons.

- [] Remove the discussion on the escape fraction.

Additional suggestions:
Only allowing one LAE per dark matter halo is a simplistic model.  You
could try more complicated HOD models where f_occ(M) is not constant
to see if they make a significant difference in the results.  But
arguing against this by citing the Jose et al. 2013 paper seems
reasonable.  You could also note the lack of close-proximity LAE pairs
noted in papers led by Nick Bond as further evidence against multiple
LAEs occurring in the same halo.  

You could explore a wide range of f_occ down to properly tiny values
such as 0.001 without increasing the run time by indexing it
logarithmically. 

Be careful when turning filters into cosmological volumes - you should
not simply use the FWHM.  See Gronwall et al. 2007 for a detailed
discussion of this. 

Be careful when comparing results for halo occupation fraction for
LAEs from a survey with a 190A observed EW threshold to the more
typical 20 A EW rest-frame, which equals 80A observed at z=3. 

SSA22 is known to be an extraordinary overdensity, as you note.  The
overdense region, and possibly the entire field, will bias the
measured number density of z=3.1 LAEs upwards.  So one should be very
cautious including this field in a study like this one. 

When using a constant value for f_occ, you could use just M_min, M_max
as parameters, choose the value of f_occ that fits the observed number
density (propagating uncertainties) and then test the resulting ACF
versus data.  This would be accurate because a random sub-sample of
halos has the same clustering as the full sample. 
