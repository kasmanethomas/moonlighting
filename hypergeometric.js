// Hypergeometric JS from http://www.math.ucla.edu/~tom/distributions/Hypergeometric.html
// See <https://gist.github.com/trevnorris/c39ac96740842e05303f>
// Annotated by Adam Novak

// Properly used with successesAvailable<popSize/2 and sampled<popSize/2

// See line 120 for usage info.

function hyp(successesObserved, sampled, successesAvailable, popSize) {
//console.log("hyp() was called with " +[successesObserved, sampled, successesAvailable, popSize]);
    // The actual hypergeometric CDF. Requires that half or less of the
    // population be successes, and that half or less of the population be
    // sampled.

    // We conceptualize the problem like this: given the population, choose the
    // successes, choose the sampled items, and look at the size of the overlap.

    // Because of this, it doesn't matter which set is the successes and which
    // is the sampled items, so we can swap those around for numerical reasons.

    // What are the sizes of the two sets we are intersecting, identifgied by
    // relative size?
    var smallerSet, largerSet;

    // best to have sampled<successesAvailable
    if (successesAvailable < sampled) {
        // The set of successes is smaller
        smallerSet = successesAvailable;
        largerSet = sampled
    } else {
        // The set of sampled items is smaller
        smallerSet = sampled;
        largerSet = successesAvailable;
    }

    // This is an intermediate value I don't really understand, which is used in
    // the middle of the cumulative CDF calculation.
    var h = 1;

    // This is the probability of having observed everything we looked at so
    // far. Except sometimes it goes above 1 and we have to fix it?
    var s = 1;

    // This is an index over which item we are at in the smaller set
    var k = 0;

    // This is an index over which item we are at in the intersection
    var i = 0;

    while (i < successesObserved) {
        // For each item in the intersection    

        while (s > 1 && k < smallerSet) {
            // Sample some items from the smaller set (?)

            // Get the probability of, after already grabbing k items for the
            // smaller set that weren't in the larger set, grabbing another
            // item for the smaller set that wasn't in the larger set.

            // Then multiply h and s by this probability

            h = h * (1 - largerSet / (popSize - k));
            s = s * (1 - largerSet / (popSize - k));

            // Advance to the next item in the smaller set.
            k = k + 1;
        }

        // Get the number of things not taken for the intersection, after taking
        // this thing: (popSize - smallerSet - largerSet + i + 1)

        // Get the number of things in the smaller set not in the intersection,
        // and multiply by the number of things in the larger set not in the
        // intersection.

        // As i goes from 0 to successesObserved, h accumulates:
        // factors of smallerSet to smallerSet - successesObserved on top
        // factors of largerSet to largerSet - successesObserved on top
        // factors of 1 to successesObserved on bottom
        // factors of (popSize - smallerSet - largerSet + 1) to 
        //    (popSize - smallerSet - largerSet + successesObserved)

        // I'm not entirely sure how this works, but it looks like we're
        // cheating a bit to calculate the sum of ratios of factorials without
        // needing to re-do lots of the multiplications.

        h = h * (smallerSet - i) * (largerSet - i) / (i + 1) / (popSize - smallerSet - largerSet + i + 1);
        s = s + h;

        // Move on to the next item in the intersection
        i = i + 1;
    }

    while (k < smallerSet) {
        // For each remaining item in the smaller set (conceptually the sampled
        // ones) that was not part of the intersection(?)

        // popSize - k is the number of items remaining to be grabbed for the
        // smaller set.

        // largerSet / (popSize - k) is the probability that the next item
        // grabbed for the smaller set would be in the larger set.

        // We take 1 - that because we know the next item grabbed for the
        // smaller set will not have been in the larger set, since it wasn't in
        // the intersection.

        // Then we multiply the probability of everything we've seen so far by
        // the probability of having successfully not added some extra thing to
        // our intersection.
        s = s * (1 - largerSet / (popSize - k));

        // Move on to grab the next item in the smaller set.
        k = k + 1;
    }

    return s;
}

// If the genome has 2600 genes, and a probe brought back 300 samples, 
// and we are looking for a theoretical maximum of 12 hits but saw 6,
// then call computeHypergeometricProbability( 2600, 300, 12, 6 )
// to get the probability of seeing 6-or-more hits.
function computeHypergeometricProbability(N, sampled, successesAvailable, successesObserved) {

    // Get the population size
    var popSize = N;

    successesObserved--;
    // we will return cumulative probability of >= successesObserved

    // This is the probability of observing this many successes or fewer, given
    // the other parameters.
    var Prob;

    if (sampled <= 0 || successesAvailable <= 0 || popSize <= 0) {
        // Nothing sampled or no successes available or nothing in the
        // population.
      //  alert("Parameters must be positive integers");
        Prob = 0;
    } else if (successesAvailable > popSize || sampled > popSize) {
        // Successes available or samples taken larger than actual population.

     //   alert("successesAvailable and sampled must be less than popSize");
        Prob = 0;
    } else if (successesObserved < 0 || successesObserved < sampled + successesAvailable - popSize) {
        // We observed negative successes, or we pulled out so much of the
        // population that we would absolutely have had to have sampled more
        // succeses than this.
        Prob = 0;
    } else if (successesObserved >= sampled || successesObserved >= successesAvailable) {
        // We saw as many successes as we had samples (or even more), or we saw
        // all (or more) of the successes available. The probability of seeing
        // that many successes *or fewer* is thus 1.
        Prob = 1;
    } else {
        // We can't short-circuit and get the answer already. We have to do real
        // math. But first, we have to precondition the input parameters to the
        // CDF function.

        if (2 * successesAvailable > popSize) {
            // More than half of the population is successes.

            if (2 * sampled > popSize) {
                // More than half the population has been sampled.

                // Change to a different problem:

                // Out of the same size population, with success and failure
                // reversed, sampling everything we didn't sample before, what's
                // the probability we sample that many or more of the things
                // that were originally unsampled failures?

                // The popSize - successesAvailable - sampled +
                // successesObserved has the +successesObserved to avoid double-
                // counting the successesAvailable that were also sampled in the
                // subtraction. So it gets the number of things that were
                // originally neither sampled nor observed.

                // P(the unsampled part contains some number of failures or
                // fewer) = P(the sampled part contains x successes or fewer)

                Prob = hyp(popSize - successesAvailable - sampled + successesObserved, popSize - sampled, popSize - successesAvailable, popSize);
            } else {
                // Half or less of the population has been sampled.

                // Change to a different different problem:

                // Out of the same size population, with success and failure
                // reversed, sampling the same number of things, what's the
                // probability of observing the same situation we actually
                // observed, but with success and failure reversed? Or something
                // with more successes than that?

                // Basically we go and get the other (smaller) tail, and do 1
                // minus that.

                Prob = 1 - hyp(sampled - successesObserved - 1, sampled, popSize - successesAvailable, popSize);
            }
        } else if (2 * sampled > popSize) {
            // Half or less of the population is successes, but we sampled more
            // than half of it.

            // Try a different problem: sample the part we didn't sample before,
            // and look at the probability of finding the corresponding number
            // of successes or more in that part, and then flip around to the
            // other tail.

            Prob = 1 - hyp(successesAvailable - successesObserved - 1, successesAvailable, popSize - sampled, popSize);
        } else {
            // We're sampling half or less of the population, and half or less
            // of it is successes. Just go straight to the actual function.
            Prob = hyp(successesObserved, sampled, successesAvailable, popSize)
        }
    }
    Prob = Math.round(Prob * 100000) / 100000;
    return 1. - Prob;
}
