# moonlighting

The JavaScript code in this project will allow you to do some bioinformatic exploration of protein-coding genes.

<p>All code is native plain-vanilla JavaScript. The code is meant to be run in the browser console. It can of course be adapted to run in a Node app or anywhere else, since it is generic JavaScript, but the workflow we recommend is just to load a Fasta CDS genome in a Chrome tab, then run code against it in the console.</p>

<h2>What's in this project?</h2>
<p>This project consists of standalone plain-vanilla JavaScript functions designed to be run in a browser's "developer console." 
Routines were tested in Chrome, but should work in other browsers as well. It's all plain ECMAScript2015 (browser JS) with no funky dependencies.</p>

<p>We use this code to explore moonlighting genes of bacteria (hence the name of the project), but a lot of the code is general-purpose and will let you do other kinds of explorations.</p>

<h2>How can I get started quickly?</h2>
<p>The easiest way to run the code is just copy it to the console and hit Enter, which will load functions into memory. 
You can, of course, paste the code into a Snippet (in Chrome dev IDE) and run it from there.</p>

<p>The assumption here is that you have a Fasta genome (an annotated CDS genome, perhaps a RefSeq genome downloaded from NCBI)
loaded, and displayed, in a browser window. (For example, to get <i>E. coli NCTC11775</i>, go to https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/690/815/GCF_000690815.1_ASM69081v1/ and download the file called <b>GCF_000690815.1_ASM69081v1_cds_from_genomic.fna.gz</b>, then unzip it in a folder, and drag it into Chrome.) The code operates against that genome. How? Easy. First, run the initialization code shown below (which is taken from <code>utils.js</code>). That is to say: With a genome already displayed
in your browser, go to your developer console using Command-Option-J (Mac) or Shift-Control-J (Win), then paste the following lines into
the console and execute them (with carriage return or Enter):</p>
<pre>
genes = document.body.textContent.toUpperCase().split(/>[^\n]+\n/);
genes.shift(); // First array item is empty, so get rid of it.
// Get rid of newlines:
genes.forEach((item,i,ar)=>{
    ar[i] = item.replace(/\n/g, "");
});
legends = document.body.textContent.match(/>[^\n]+\n/g);
</pre>

Now you have every gene in a global array called 'genes' and you have every Fasta header in an array called 'legends'.

If you run 'genes.length', you should see the number of genes of your genome in the console window.

If you run <code>console.log(legends[0])</code>, you should see the first Fasta header printed out on the console.

<h2>What's in the files?</h2>
<p>The <code>utils.js</code> file contains some handy functions like average(), sd() for Standard Deviation, etc., plus some interesting utilities like the CodonPercent() function, which takes a gene sequence and a regex for arguments.</p>

<p>The <code>hypergeometric.js</code> file contains a handy routine for computing cumulative hypergeometric probabilities.</p>

<p>The <code>pearson.js</code> file contains a function that lets you obtain the correlation between values in two equal-length arrays.</p>

<p>The <code>proximity.js</code> file contains routines that let you do enrichment experiments to find out what kinds of genes exist near moonlighting genes, on the genome in question. In other words, suppose you have an array of gene indexes for presumed moonlighting genes. Suppose you want to know: What kinds of genes are within plus or minus 5 genes of the moonlighters? It's easy enough to obtain the list of neighbor-genes, but the question, now, is: Do the neighbor genes over-represent certain functional categories? Are there a disproportionate number of cell-membrane genes, for example? Are there a disproportionate number of genes involved in secretion? To get a statistically meaningful answer to these sorts of questions requires careful analysis of the numbers, using hypergeometric probability analysis. To help with functional categories, we've gathered the names of genes associated with various GO (Gene Ontology) labels. See <code>proximity.js</code> for details.</p>

<p>The code>antisenseORF.js</code> file contains a function that allows you to search a group of genes in antisense for putative open reading frames, then print (to the console) results as Fasta listings. Each listing has a descriptive header and the polypeptide translation of the region in question. A usage example is given at the bottom of the file.</p>


<h2>What kinds of quick experiments can I do?</h2>
<p>Suppose you've got E. coli's genome loaded and you want to know what percentage of genes use 'TGA' as a stop codon.</p>
<pre>
// Run utils.js first, to load the sum() function. Then, run this in the console:
sum( genes.map(gene=>gene.slice( gene.length - 3 )=='TGA' ? 1:0) )/genes.length
</pre>

In <i>E. coli</i>, the answer gets printed to the console: <code>0.30580075662042877</code> Almost 31% of genes use TGA as a stop codon.

<p>If you load <code>utils.js</code>, you can make use of a function called CodonPercent(). Let's use it to determine the percentage of codons that begin with a purine:</p>
<pre>
average( genes.map(gene=> CodonPercent(gene, /[AG]../g) ))
// 0.6004191707440053
</pre>
<p>In E. coli, 60% of codons begin with A or G.</p>

<p>The <code>utils.js</code> file contains other examples of console experiments you can run against any genome. Some of the other files also contain useful examples.






 









