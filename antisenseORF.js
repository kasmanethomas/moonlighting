// Get Fasta printout of all antisense products of moonlighters.
// Depends on utils.js.

// See usage example at the bottom.


/* getFastaRun()

Produces a listing, to the console, of translation products of putative antisense ORFs.
Call it with an array of ORF regexes. The regexes will look something like:  

                        /.{12}(..ATG)(...)+?(?=(TAA|TAG|TGA)).../g 

This particular regex says to find a 12-base leader followed by two bases, followed by ATG (start codon), 
followed by one or more codons (...) which DO NOT contain a stop codon, followed by three of any base.
The "three of any base" at the end will capture the stop codon.

The unusual construction ?(?=(TAA|TAG|TGA)) prevents capture of a stop codon while we're looking for codons.
But eventually, the search finds a stop codon and can't continue, so it will end, and then the final three bases
will be the actual stop codon itself.

The leader is obviously optional. Use it (or not) to make it easier to look for Shine Dalgarno sequences.

The expression tries to find the longest stop-free ORF. This means a hit can be short or long, depending
on where the stop codon is found. It's up to you to filter results by some length requirement, if desired.
Use the fourth argument for this.

ARGUMENTS:
moonies: An array of gene indexes (integers).
regexArray: An array of regex-literals.
offsetToStartCodon: If your regexes contain a leader section, this value is the size of all leader characters preceding the start codon.
aaMinLength: Minimum number of amino acids in the polypeptide encoded by the antisense ORF.

See usage example at end of code, below.

*/

function getFastaRun( moonies, regexArray, offsetToStartCodon, aaMinLength ) {

    let MINIMUM_POLYPEPTIDE_LENGTH = aaMinLength;
    let offsetToStart = offsetToStartCodon || 14;
    getHits = (g,rx)=>{
        let hits = reverseComplement(g).match(rx);
        if (!hits)
            return [];
        hits.sort((a,b)=>b.length - a.length);
        let results = [];
        hits.forEach(h=>{
            if (h.length > 29)
                results.push(h);
        }
        );
        return results;
    }

    // Save the world here
    let ar = [];
   
    proteins = {};

    regexArray.forEach(regex=>{
        for (var i = 0; i < genes.length; i++) {

            let gi = i;

            if (moonies.indexOf(gi) == -1)
                continue;

            let hits = getHits(genes[gi], regex);
            if (!hits)
                continue;
            hits.forEach((h,index)=>{
                let polypeptide = translate(h.slice(offsetToStart));
                let indexOfHit = reverseComplement(genes[gi]).indexOf(h.slice(offsetToStart));
                let mod = indexOfHit % 3;
                if (polypeptide.length < MINIMUM_POLYPEPTIDE_LENGTH)
                    return;

                // Persist all data in a global object 
                if (!(polypeptide in proteins)) {

                    proteins[polypeptide] = ({
                        offset: indexOfHit,
                        bases: h.slice(offsetToStart).length,
                        startCodon: h.substr(offsetToStart,3),
                        gi: gi,
                        name: getProteinName(legends[gi]),
                        mod: indexOfHit % 3,
                        ppSD: PurinePercent(h.slice(0, 15)),
                        R1: CodonPercent(h.slice(14), /[AG]../g)
                    });

                    // Create a Fasta listing
                    ar.push('>' + mod + "_" + Math.random().toString().slice(5, 10) +
                        "_" + getProteinName(legends[gi]) + '-' + index +  '\n' + polypeptide);

                }            
            } ); // hits.forEach
        }
        // for genes
    }
    );
    // regexArray.forEach

    console.log(ar.join('\n'));
    
    console.log("Total translation products: " + ar.length);

}


// EXAMPLE using a bunch of regexes reflecting different start codons that are actually used in Mtb
let asORF_regexes = [   
                        /.{12}(..TTG)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..ATG)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..GTG)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..CTG)(...)+?(?=(TAA|TAG|TGA)).../g,  
                        
                        /.{12}(TGATG)(...)+?(?=(TAA|TAG|TGA)).../g,
    
                        /.{12}(..ATGA..)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..TTGA..)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..GTGA..)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..CTGA..)(...)+?(?=(TAA|TAG|TGA)).../g,
    
                        /.{12}(..ATT)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..ATA)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..TAC)(...)+?(?=(TAA|TAG|TGA)).../g,
                        /.{12}(..ATC)(...)+?(?=(TAA|TAG|TGA)).../g                 
];


// RUN IT
function getTargets(rx) {
    let targets = [];
    genes.forEach((g,gi)=>{
        if (legends[gi].match(rx))
            targets.push(gi);
    } );
   return targets;
}
let moonies = getTargets(moonlit);
let offsetToStartCodon = 14;
let aaMinLength = 30; 
// This should produce scores (or even hundreds) of hits in any organism:
getFastaRun( moonies, asORF_regexes, offsetToStartCodon, aaMinLength );

