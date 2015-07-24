rm03s06 is the best center
rm04s06 is the 2nd best

the richness is in field lambda_chisq

These files have the sums
    wsum - sum of weights
    dsum - sum of \Delta\Sigma*weight
    rsum - sum of radius
    npair - count of lens-source pairs

When you do a stack you should sum those fields for everything
on the stack.

For either a stack or a single cluster, you can get means and errors
in the following way.

    To get the mean delta sigma
        dsum/wsum

    To get the mean radius
        rsum/npair 

    The weights are approximately 1/err^2 where err is the error in delta
    sigma.  Thus to get an approximate error estimate on delta sigma

        sqrt(1/wsum)

    This works for the ensemble too.
