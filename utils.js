// For Hypergeometric Probability computations, use the code at
// https://gist.github.com/adamnovak/f34e6cf2c08684752a9d

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



// Moonlighting genes by name, converted also to a regular expression
moonlitSource = ['methyltransferase Erm', 'glutamate racemase', 'elongation factor Tu', '[^l]malate synthase|glcB', 'cysteine desulfurase', 'gamma-glutamyl phosphate reductase', 'glucose-6-phosphate isomerase', '6-phosphofructokinase', '6-phosphogluconate dehydrogenase', 'enolase|phosphopyruvate hydratase', 'triose-?phosphate isomerase', 'fusA', 'pepO', 'rpoB', 'DnaK', 'GroEL', 'Ag85[ABC]', 'superoxide dismutase', '23S rRNA Methyltransferase', 'glyceraldehyde.3-phosphate dehydrogenase', 'phosphoglycerate kinase', 'fructose-bisphosphate aldolase', 'phosphoglycerate mutase', 'glutamine synthetase'];
moonlit = new RegExp(moonlitSource.join("|"),'g');

// A routine that gets the percent of codons that match a regex.
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


// Use this to convert codons to Dayhoff codes:
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
    // Creating the mean with Array.reduce
    let mean = arr.reduce((acc,curr)=>{
        return acc + curr
    }
    , 0) / arr.length;

    // Assigning (value - mean) ^ 2 to every array item
    arr = arr.map((k)=>{
        return (k - mean) ** 2
    }
    )

    // Calculating the sum of updated array
    let sum = arr.reduce((acc,curr)=>acc + curr, 0);

    // Calculating the variance
    let variance = sum / arr.length

    // Returning the Standered deviation
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
transporter = new Test("transporter",/transporter.[^p]/g);
// transporter not followed by 'permease'
mycolic = new Test("mycolate synthesis",/mycol[yia]/g);
fattyacid = new Test("fatty-acid synthesis",/fatty/g);

// Put a bunch of tests into an array:
SpecialTestArray = [polyketide, peptidoglycan, rf, transmembrane, ribosomal, transporter, new Test("moonlighting",moonlit), pgrs, membrane, secretion, efflux, hypothetical, ppe, permease, esx, trna, fattyacid, mycolic];




