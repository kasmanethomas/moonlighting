/*
Load a Fasta CDS genome (e.g. a RefSeq genome from NCBI) into a tab of Chrome.
Then run the following code in the Chrome developer console:
*/

// ======= INITIALIZE GLOBAL ARRAYS ======
// The next 7 lines will create a global array of genes
// plus a global array of Fasta headers called 'legends'.

genes = document.body.textContent.toUpperCase().split(/>[^\n]+\n/);
genes.shift(); // First array item is empty, so get rid of it.
// Get rid of newlines:
genes.forEach((item,i,ar)=>{
    ar[i] = item.replace(/\n/g, "");
});
legends = document.body.textContent.match(/>[^\n]+\n/g);

/*
Moonlighting genes have no specific Gene Ontology designation. Therefore, we 
look for them individually by name (see below).
*/

// Moonlighting genes by name:
moonlitSource = ['methyltransferase Erm', 'glutamate racemase', 'elongation factor Tu', '[^l]malate synthase|glcB', 'cysteine desulfurase', 'gamma-glutamyl phosphate reductase', 'glucose-6-phosphate isomerase', '6-phosphofructokinase', '6-phosphogluconate dehydrogenase', 'enolase|phosphopyruvate hydratase', 'triose-?phosphate isomerase', 'fusA', 'pepO', 'rpoB', 'DnaK', 'GroEL', 'Ag85[ABC]', 'superoxide dismutase', '23S rRNA Methyltransferase', 'glyceraldehyde.3-phosphate dehydrogenase', 'phosphoglycerate kinase', 'fructose-bisphosphate aldolase', 'phosphoglycerate mutase', 'glutamine synthetase'];

// Convert the above list into a single regular expression:
moonlit = new RegExp(moonlitSource.join("|"),'g');

// A utility routine that gets the percent of codons that match a regex.
// The first arg is a string (representing a DNA sequence).
CodonPercent = (g,codonrx)=>{
    var t = 0;
    var codons = g.match(/.../g); // get ALL codons
    codons.forEach(c=>{
        var m = c.match(codonrx);
        if (m)
            t++;
    });
    return t / codons.length;
};

// So now we can do things like test for codons of pattern RNY:
rny = (gi)=>CodonPercent(genes[gi], /[AG].[CT]/g);

// Other common patterns we might want to look for:
nyn = (gi)=>CodonPercent( genes[gi], /.[CT]./g );
rnn = (gi)=>CodonPercent( genes[gi], /[AG]../g );
nny = (gi)=>CodonPercent( genes[gi], /..[CT]/g );

// Use this table to convert codons to Dayhoff codes:
Codons = {
    "CCA": "P",
    "CCG": "P",
    "CCT": "P",
    "CCC": "P",
    "CGA": "R",
    "CGG": "R",
    "CGT": "R",
    "CGC": "R",
    "CAA": "Q",
    "CAG": "Q",
    "CAT": "H",
    "CAC": "H",
    "CTA": "L",
    "CTG": "L",
    "CTT": "L",
    "CTC": "L",
    "GCA": "A",
    "GCG": "A",
    "GCT": "A",
    "GCC": "A",
    "GGA": "G",
    "GGG": "G",
    "GGT": "G",
    "GGC": "G",
    "GAA": "E",
    "GAG": "E",
    "GAT": "D",
    "GAC": "D",
    "GTA": "V",
    "GTG": "V",
    "GTT": "V",
    "GTC": "V",
    "ACA": "T",
    "ACG": "T",
    "ACT": "T",
    "ACC": "T",
    "AGA": "R",
    "AGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AAA": "K",
    "AAG": "K",
    "AAT": "N",
    "AAC": "N",
    "ATA": "I",
    "ATG": "M",
    "ATT": "I",
    "ATC": "I",
    "TCA": "S",
    "TCG": "S",
    "TCT": "S",
    "TCC": "S",
    "TGA": "*",
    "TGG": "W",
    "TGT": "C",
    "TGC": "C",
    "TAA": "*",
    "TAG": "*",
    "TAT": "Y",
    "TAC": "Y",
    "TTA": "L",
    "TTG": "L",
    "TTT": "F",
    "TTC": "F"
};

// Using the above table, we can translate a gene in silico.
// Note: This routine will honor stop codons. It returns only the amino acid sequence preceding the stop.
translate = (g)=>{
    let rcg = (g);
    let codons = rcg.match(/.../g);
    if (codons == null || codons.length < 1)
        return "";
    var rprotein = "";
    codons.forEach(c=>{
        rprotein += Codons[c];
    });
    
    // If a stop codon occurred, just return the first translation product
    rprotein = rprotein.split("*")[0];
    return rprotein;
};

// Sum an array.
sum = (ar)=>{
    let s = ar.reduce((s,item)=>s + item, 0);
    return s;
};
// EXAMPLE:
// sum( [1,2,3,4])
// 10

// AVERAGE
function average(data) {
    var sum = data.reduce(function(sum, value) {
        return sum + value;
    }, 0);

    var avg = sum / data.length;
    return avg;
}

// STANDARD DEVIATION
function sd(arr) {
    // Create the mean with Array.reduce
    let mean = arr.reduce((acc,curr)=>{
        return acc + curr
    }
    , 0) / arr.length;

    // Assign (value - mean) ^ 2 to every array item
    arr = arr.map((k)=>{
        return (k - mean) ** 2
    }
    )

    // Calculate the sum of array
    let sum = arr.reduce((acc,curr)=>acc + curr, 0);

    // Calculate the variance
    let variance = sum / arr.length

    // Return the Standard deviation
    return Math.sqrt(sum / arr.length)
}

// A simple Test object containing a name string and a regex.
function Test(name, rx) {

    this.name = name;
    this.rx = rx;
}

// Some Test objects that allow to find keywords in Fasta headers.
pgrs = new Test("PE-PGRS",/PE-PGRS/g);
ppe = new Test("PPE family",/PPE family/g);
membrane = new Test("membrane",/[^s]membrane/g);
// membrane not preceded by 'trans'
transmembrane = new Test("transmembrane",/transmembrane/g);
efflux = new Test("efflux",/efflux/g);
permease = new Test("permease",/permease/g);
esx = new Test("esx",/Esx/g);
sigma = new Test("Sigma factors",/Sig[A-Z]/g);
antisigma = new Test("anti-Sigma factors",/anti-sigma/g);
secretion = new Test("secretion",/secretion/g);
hypothetical = new Test("hypothetical protein",/hypothetical/g);
trna = new Test("tRNA ligase",/tRNA ligase/g);
polyketide = new Test("polyketide synthase",/polyketide/g);
rf = new Test("release factor",/release factor/g);
peptidoglycan = new Test("peptidoglycan",/peptidoglycan/g);
trna = new Test("tRNA ligase",/tRNA ligase/g);
ribosomal = new Test("ribosomal protein",/ribosomal protein/g);
// transporter not followed by the word 'permease':
mycolic = new Test("mycolate synthesis",/mycol[yia]/g);
transporter = new Test("transporter",/transporter.[^p]/g);

fattyacid = new Test("fatty-acid synthesis",/fatty/g);

// Put a bunch of tests into an array:
SpecialTestArray = [polyketide, peptidoglycan, rf, transmembrane, ribosomal, transporter, new Test("moonlighting",moonlit), pgrs, membrane, secretion, efflux, hypothetical, ppe, permease, esx, trna, fattyacid, mycolic];

/* getProteinName()
   Gets protein descriptor from Fasta legend.
*/
getProteinName = (legend)=>{
    var rx = /protein=([^\[\]])+/;
    if (legend.match(rx))
        return legend.match(rx)[0].replace("protein=", "");
    return "(unknown)";
}

/* getGeneName()
   Gets gene name from Fasta legend using gene index (gi)
*/
getGeneName = (gi)=>{
    if (isNaN(1 * gi) || typeof legends == 'undefined')
        return "(unknown)";
    let name = legends[1 * gi].match(/gene=(\w+)/);
    if (name)
        return name[1];
    return "";
}
// EXAMPLE
/*
getGeneName(0)
'dnaA'
*/

// ENRICHMENTS ROUTINE
/*
The first argument is an array of tab-delimited data items in which the first column is the gene index
(so for example, in Mtb, gene zero is dnaA) and the second column can be a gene's description from the 
"protein=" part of the Fasta header. E.g., [protein=translation initiation factor IF-2] 
See example further below.

The second arg is an array of Test objects (see definition further above).
*/
function enrichmentReport(ar, TA) {

    let getGeneIndex=(item)=>item.split('\t')[0] * 1;
    let headerHits=(gi,rx)=>(legends[gi].match(rx) || []).length; 

    let result = [];
    TA.forEach(test=>{
        let rx = test.rx;
        var possible = sum(legends.map(m=>m.match(rx) ? 1 : 0));
        var actual = sum(ar.map(m=>headerHits( getGeneIndex(m), rx )) );
        var successRatio = (actual && possible) ? actual / possible : 0;
        var coverageRatio = ar.length / genes.length;
        var enrichment = successRatio / coverageRatio;
        var hyp = computeHypergeometricProbability(genes.length, ar.length, possible, actual);
        result.push([enrichment.toPrecision(3), hyp.toFixed(3), test.name, actual + "/" + possible].join('\t'));
    });

    // sort results by fold enrichment
    result.sort((a,b)=>{
        return b.split('\t')[0] * 1 - a.split('\t')[0] * 1
    });

    // find out how many enrichments were >= 1.0
    var t = 0;
    result.forEach(r=>{
        if (r.split('\t')[0] * 1 >= 1)
            t++;
    });

    // return just enrichments that were >= 1.0
    return result.slice(0, t).join('\n');
}

// ENRICHMENT ASSAY (see example following).
// Callback (first arg) should take a gene index as the sole arg.
// 'fraction' should be 0.1, or 0.2, etc.
// This function returns no value. It writes to the console.
function enrichmentAssay(cb, fraction) {
    let genesSortedByMetric = genes.map((m,gi)=>gi + '\t' + cb(gi) + '\t' + getProteinName(legends[gi])).sort((a,b)=>{
        return b.split('\t')[1] * 1 - a.split('\t')[1] * 1
    }
    );

    // We look at enrichments in the top fraction
    let er = enrichmentReport(genesSortedByMetric.slice(0, Math.floor(genes.length * fraction)), SpecialTestArray, 1);
    console.log("Fold\tE\t\tFunction (actual/possible)");
    console.log(er);
}

// EXAMPLE
enrichmentAssay( rny, .2 );

In Mtb H37Rv, the result is:
Fold	E		Function (actual/possible)
4.67	0.000	PE-PGRS	56/60
3.00	0.058	peptidoglycan	3/5
2.29	0.000	moonlighting	16/35
1.83	0.002	PPE family	23/63
1.82	0.056	esx	8/22
1.67	0.488	efflux	1/3
1.21	0.258	ribosomal protein	14/58
1.02	0.500	transporter	17/83
*/








