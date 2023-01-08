/*
This code supports an experiment that asks: What are the neighbors (on the genome) of moonlighting genes?
What kinds of genes are in close proximity to moonlighters? 

At the very bottom, two lines of code are given, to show how to run an experiment. The experiment puts a 
brief report (of enrichment results) in the console.

To obtain enrichment results, we rely on lists of gene names associated with various functional categories; categories with names
like "cell division" or "inner membrane." These are the plain-English names for broad categories of function.

To obtain the gene names associated with (for example) "inner membrane," we go to the GO site at http://amigo.geneontology.org/
and search for "inner+membrane", assigning a taxon_subset_closure_label: Bacteria, and	isa_partof_closure_label: membrane.

The various search terms and results are listed below, throughout. Gene names (e.g. 'secA2') are agglomerated into regular expressions.
The regexes are used in Test objects (see utils.js or just read below). An array of Test objects gets passed to the enrichment() function
along with a data array in which every data item is simply a gene index, followed by a tab, followed by the gene's 
'protein=' Fasta descriptor. 

*/

// Run utils.js before using these routines.

// Utility routine for obtaining the gene indexes
// of genes whose Fasta legends match regex 'rx'
function harvestGeneIndexesUsingPattern(rx) {
    let harvested = [];
    genes.forEach((g,gi)=>{
        if (legends[gi].match(rx))
            harvested.push(gi);
    }
    );
    return harvested;
}


// Filter an array of type [gi + '\t' + data]
filterPA = (pa,filterRx)=>{
    let newArray = [];
    pa.forEach(p=>{
        if (!(p.match(filterRx)))
            newArray.push(p);
    }
    );
    return newArray;
}

// neighbor analysis
neighbors = (giArray,dist)=>giArray.map(m=>(m + dist) + '\t' + getProteinName(legends[m + dist]));
function getNeighbors(ar, toDistance) {
    let allNeighbors = [];
    for (var i = 1; i <= toDistance; i++)
        allNeighbors = allNeighbors.concat(neighbors(ar, -i)).concat(neighbors(ar, i));

    // Allow no duplicates:
    let lut = {};
    let results = [];
    allNeighbors.forEach(a=>{
        if (!(a in lut)) {
            lut[a] = 1;
            results.push(a);
        }
    }
    );
    return results;
}

function enrichments(ar, TA) {

    let getActualHitsCount = (ar,rx)=>{
        let count = 0;
        ar.forEach(item=>{
            let gi = item.split('\t')[0] * 1;
            let hits = (legends[gi].match(rx) || []).length;
            if (hits != 0)
                count += hits;
        }
        );

        return count;
    }

    let result = [];
    TA.forEach(test=>{
        let rx = test.rx;
        var possible = test.name.indexOf("(GO)") != -1 ? sum(legends.map((m,gi)=>getGeneName(gi).match(rx) ? 1 : 0)) : sum(legends.map(m=>m.match(rx) ? 1 : 0));
        var actual = test.name.indexOf("(GO)") != -1 ? sum(ar.map((m)=>getGeneName(m.split('\t')[0] * 1).match(rx) ? 1 : 0)) : getActualHitsCount(ar, rx);
        var successRatio = (actual && possible) ? actual / possible : 0;
        var coverageRatio = ar.length / genes.length;
        var enrichment = successRatio / coverageRatio;
        var hyp = computeHypergeometricProbability(genes.length, ar.length, possible, actual);
        result.push([enrichment.toPrecision(3), hyp.toPrecision(4).slice(0, 5), test.name, actual + "/" + possible].join('\t'));
    }
    );
    result.sort((a,b)=>{
        return b.split('\t')[0] * 1 - a.split('\t')[0] * 1
    }
    );
    var t = 0;
    result.forEach(r=>{
        if (r.split('\t')[0] * 1 >= 1)
            t++;
    }
    );

    return result.slice(0, t).join('\n');

}
// enrichments()

/* In the code below, we use the Test() object a lot. It is defined in util.js as:

function Test(name, rx) {
    this.name = name;
    this.rx = rx; // rx means a regex
}

*/

/* cell+wall
+	taxon_subset_closure_label: Mycobacterium tuberculosis H37Rv	
+	isa_partof_closure_label: cell wall biogenesis
*/
GOcellWallOrganization = new Test("cell wall biogenesis (GO)",/wag31|pks7|Rv0519c|Rv0774c|Rv1639c|Rv1288|glf|fadD32|lppX|Rv2952|pks1|aftA|embA|fadD23|pks10|glfT1|papA5|mraY|drrB|accD4|Rv2953|fadD15|dprE1|accD3|dprE2|mmaA1|accD2|mmpL3|mmaA4|mmaA2|kasB|murE|wecA|pks11|pks13|mas|tesA|ddl|ldtA|ponA2|ldtB|Rv1433|murC|Rv2949c|glfT2|ppsE|fbpC|espA|inhA|fadD22|fadD29|fadD30|whiB3|fadD19|mmaA3|mabA|cmaA1|pcaA|cmaA2|ppsA|ppsD|ppsC|glmM|Rv2959c|ppsB|murI|alr|fadD17|fadD26|fadD28|pknB/);

/* cell+wall
taxon_subset_closure_label: Mycobacterium tuberculosis H37Rv	
+	isa_partof_closure_label: transmembrane transport
*/
GOtransmembraneTransport = new Test("transmembrane transport (GO)",/drrB|iniA|RVBD_0194|espA|esxB|espD|eccCa1|eccCb1|espC|drrA|arfA|secA1|esxA|secA2/);

GO_permease = new Test("permease (GO)",/pstA1|Rv1686c|yrbE1A|Rv2040c|irtB|RVBD_0194|proZ|proW|Rv2687c|sugB|sugA|Rv0072|Rv1283c|Rv1747|Rv0987|cycA|Rv0180c|yrbE4A|yrbE4B|yrbE2A|yrbE1B|Rv1672c|yrbE2B|yrbE3B|yrbE3A|modB|rfbD|Rv2281|Rv1272c|rocE|cysW|cysT|Rv3253c|Rv2686c|ansP1|ansP2|irtA|bacA|Rv1282c|drrB|drrC|Rv2563|pstC2|glnH/);

/* cell+division
*/
GOcellDivision = new Test("cell division (GO)",/ftsX|ftsQ|ftsZ|crgA|sepF|whiA|ftsE|wag31|fic|ftsW|rodA|fhaB|Rv3717|Rv0530|clpX|pstS3/);

GOtransmembrane = new Test("transmembrane (GO)",/Rv1226c|nirD|Rv1648|esxA|esxB|eccCa1|Rv2040c|secA2|Rv3000|dppA|ugpC|irtB|Rv3779|Rv0556|Rv3435c|Rv1069c|Rv0200|RVBD_0194|Rv1567c|proZ|proW|eccCb1|Rv2687c|espA|Rv1410c|spmT|iniA|mmr|Rv1428c|Rv1634|arsC|Rv3645|Rv3683|Rv3041c|Rv0473|sugI|Rv3239c|Rv3104c|Rv3454|Rv1672c|nanT|Rv0235c|Rv0318c|chaA|Rv3848|Rv3277|Rv3271c|narK3|narK1|kdpB|lysE|Rv0488|ftsW|kdpC|Rv2287|nuoB|nuoK|nuoA|nuoM|nuoF|nuoE|pit|Rv2281|Rv2025c|Rv1877|Rv0205|modC|Rv0073|Rv1272c|Rv2564|amt|ctpE|ctaE|crcB1|crcB2|phoP|cysW|Rv0996|mmpL7|mscL|espB|espC|espD|espI|Rv2686c|mntH|mctB|arfA|Rv2536|drrA|Rv2688c|irtA|atpE|drrB|Rv1101c|tap|phoU1|phoU2|pknF|Rv1280c|secA1|Rv0502|Rv0226c|Rv0219/);

gs = 'pilC|pilB|pilA|CBU_0056|secA|CBU_0532|tamA|gspM|yghD|epsE|BC_0044|gll0920|gll0921|glr2438|gll2495|glr3151|gll3283|PA5210|escB|c3393|escU|escC|E2348C_3947|sepQ|sctN|espA|escF|espF|CPS_3928|PFL_3778|trwE|trwI|trwK|avrA|hrpW1|VC_A0111|VC_A0108|pcrV|pcr1|sseL|srcA|outE|outC|outB|etpG|virB4-1|sscA|sifA|PSPPH_B0026|popP2|sopD|pilD|virD4|E2348C_1442|tagF|tssM|virB11|PSPTO_5426|PSPTO_5421|prtE|prtD|hlyB|SO_2234|PSPPH_2526|virD5|eccB3|eccB5|eccD1|eccD3|eccD4|eccD5|esxA|esxB|GSU2609|eccCa1|lcrH_2|eccC3|eccC4|eccC5|spa40|PFL_4660|ssaQ|sctL2|ERDMAN_4251|DVUA0106|DVUA0111|ERDMAN_4243|eccC2|ERDMAN_4245|virB6-4|A0A0N9NCT8|VC_A0284|VC_A0120|VC_A0118|VC_A0107|pYV0082|gspJ|HNE_2704|espJ|trwJ|espG5|nleA|fimU|fliS|lepB|icmQ|slrP|esaA|ECH_1042|yscO|ECH_0772|eccA5|MCA2888|virB2-3|virB2-6|GTNG_0418|CBU_1686|dotC|icmL.2|icmK|MCA2765|APH_1068|secA2|cdib|xcpQ|ppdD|ybaV|ytfM|YPO0811|YPO2044|bll3872|xcpX|xcpW|xcpV|xcpS|xcpR|gspI|gspF|gspE|gspD|tssH|Dtur_1103|Dtur_1035|Dtur_0873|tcpE|aprF|aprD|YPO3523|hofB|YPO1015|YPO0815|CBO1001|YPO2999|CBO2539|sspH2|SPO0867|ybaY|evpC|eseG|PSPPH_2180|CT_670|PSPPH_2516|PSPPH_0725|sseC|sseA|EC042_4539|EC042_4544|gspE1|exeA|PSPTO_5427|PSPTO_5423|PSPTO_5420|acrH|aopB|tatA|virB8|tssG|PSPPH_0124|PSPPH_0133|hrcN|ssaE|sscB|PSPPH_4224|pcr3|popN|outF|outG|sseJ|outH|outI|outJ|outK|outL|outM|PSPPH_0134|eibA|eibC|gldM|gldL|gldK|CPn0706 Homolog|yscN|HNE_1858|pipB2|epsL|epsD|DET1352|DET1354|pipB|DET1514|sctC|mxiC|PSPPH_2825|essB|secB|NSE_0787|NSE_0859|NSE_0860|virB4-2|NSE_0768|virB2|virB10|escJ|MCA1117|virB9|virB7|virB6|virB4|virB3|pilE|mxiM|eccCb1|eccD2|eccB1|sctN1|ECH_0989|ECH_1047|secE|ffh|yscD|sseD|espR|espK|eccE2|mshG|mshA|mshD|MCA0061|secD|ECH_0187|ECH_0252|PFL_2820|MCA2226|SPO0300|SPO3814|SPO0790|eibD|SPO0329|tamB|yueB|epsF|epsJ|epsI|epsH|CBOP19|BT_3641|HI_1008|SO_4146|SO_4147|aggC|aggA|tama|hofC|hxuB|celA|sll0547|CBO2978|FN1611|FN1818|FN2086|FN2095|FN0131|FN0292|xcsF|xcsD|secA1|pulQ|Dtur_1168|pilE1|XCC1480|ssaC|NMB1780|NMB0299|VC_0857|SCP1.169|CBO1909|PA2673|PA2463|PA2543|PA2542|PA2676|PA2677|comGA|YPO0273|lmo0708|comEA|ytfN|ppdB|PA1382|pscC|PA0692|PA0680|PA0677|pscT|YPO0272|BC_4239|BC_4324|BC_1639|nfrB|flhB|flhA|csgG|csgF|csgE|gspO|xpsD|hrcC|xcsE|xcsI|xcsJ|XCC3098|XCC4085|XCC4084|fhaC|yscC|TM_0837|ssaU|ssaT|TM_1179|TM_1117|TM_1052|Caur_1057|Caur_3127|Caur_2411|DR_0774|DR_0207|SCO2568|DR_1855|DR_1964|YPO0690|YPO2491|YPO0810|YPO0266|NMB0496|SAOUHSC_01693|NMB0057|SAOUHSC_01643|NMB1657|NMB1762|NMB2017|VC_A1084|SCO5009|VC_0199|VC_2423|SCO4526|SCO4254|VC_1446|hpmB|tapB|ppdD3|aq_1286|pilB1|THEYE_A0597|pilB2|THEYE_A0213|THEYE_A0212|rhcC2|blr3253|bll0666|blr0996|ctpG|ctpD|YPO0816|YPO2280|RB10560|sll1386|RB11294|RB11961|pilT|RB13073|RB13104|outD|RB1125|RB1147|RB3479|RB3482|RB3524|YPO3524|RB4151|tadA|SCO3556|blr7872|bll6020|bll6021|blr6022|blr6126|bll7021|rcpA|PA4540|PA4624|PA0040|PA3140|PA0685|hxcR|PA0687|pscU|xqhA|Caur_0512|sll1921|pilF|pilG|NMB1737|xpsF|xpsI|RB11595|hlyD|cpaC|YPO4005|clyB|blr3697|blr3500|blr3494|eccA2|RB9390|ssaV|mttC|gll3884|LA_1510|LA_1501|piv|vgrG3|PA0572|fliQ|rsmY|vgrG4|tse5|cbpD|plcH|pscL|pscK|pscI|pscH|pscG|pscF|exsA|exsE|exsC|pcrH|pcrG|pcrD|pcr4|pcr2|PA0099|ctpA|orfX|secG|dnaK|PA2939|exsD|pscB|pscD|pscE|secF|pscJ|lipC|mexR|fimV|xcpP|xcpT|xcpU|fliR|pscO|pscP|pscQ|vgrG4b|stk1|stp1|vgrG6|icmF2|dotU2|hsiJ2|lip2|fha2|sfa2|clpV2|hsiH2|hsiG2|hsiF2|hsiC2|hsiB2|hsiA2|PA0095|tse6|tsi6|vgrG1|tssA1|ppkA|lipA|rbsC|lon|cdrA|cdrB|plcB|ptrA|pppA|xcpY|xcpZ|aprE|tsi5|cagA|mycP1|VC_1798|virB2-1|phoP|D3G36_10545|VC_2424|sspH1|VC_2627|vgrG2|VC_A0105|VC_A0110|VC_A0113|VC_A0114|VC_A0119|VC_A0121|yadA|prtF|xpsE|VC_0402|VC_0405|VC_0406|tetX|sctN2|APH_0792|CBU_1198|APH_0924|GBAA_2190|virB2-8|virB2-2|cagY|icmB|icmD|icmE|icmS|HP_0753|virB6-2|virB6-3|APH_0256|CBU_1063|APH_0751|yggR|gspH|yghE|fliP|gspL|gspK|gspC|gspA|fliI|gspB|hofQ|gspG|psrA|ampR|vgrG4a|mexT|mexS|pdtB|pdtA|hxcW|hxcU|hxcY|hxcX|hxcT|hxcV|hxcZ|hxcQ|hxcS|lapA|lapB|hxcP|secY|hsbR|plpD|typA|plcN|vgrG2a|mep72|lasA|xphA|rsmZ|rsmA|lasB|vgrG2b|tspR|vgrG5|phoA|tle5b|epsM|mshM|mshL|impL|eccE1|eccE5|espB|espC|espD|espE|espI|espL|mshC|mshE|MCA0060|MCA0091|virB1|virB5|st20|MCA1114|MCA1118|MCA1119|MCA1121|MCA1122|MCA1123|virB9-2|MCA1277|eccA3|eccA1|VP1683|SO_0854|SO_1176|eccE3|tklG|GSU3167|tssB|tssC|tssD|fimT|SO_3523|SO_3781|SO_4159|mxdC|PFL_0092|espH|EC042_4533|EC042_4529|eccC|tcpT|espG1|hrpN|oxpG|pulG|pulP|pulE|pulF|GSU1982|GSU2221|sseB|prgO|PSPPH_4980|PSPPH_B0047|PSPPH_B0025|PSPPH_B0027|CPS_3919|sctL|BTH_II0824|CPS_4163|cag-alfa|CPS_4623|PSPTO_2543|PSPTO_2544|PSPTO_2546|PSPTO_2547|PSPTO_2550|PSPTO_2551|PSPTO_2552|PSPTO_2554|sctC1|PSPTO_2700|damX|sopD2|PSPPH_2591|steC|PSPPH_0127|gspN|PSPTO_4345|tssJ|tssK|PSPTO_4633|tssF|GSU0435|GSU0985|cagT|cagX|PSPTO_5418|PSPTO_5430|PSPTO_5432|PSPTO_5433|PFL_0701|GSU1536';
GOsecretions = new Test("secretion (GO)",new RegExp(gs));

/*
inner membrane
	taxon_subset_closure_label: Bacteria	
+	isa_partof_closure_label: membrane
*/
GOinnerMembrane = new Test("inner membrane (GO)",/arnT|YPO2181|yciB|PA5486|creD|PA4857|yadS|yccF|ybaN|ydcZ|ychE|yhaH|yiiR|marC|YPO1439|YPO1260|yoaE|STM0575|BT_1254|yejM|FN1985|Dtur_0934|XCC2498|yjdB|Caur_1084|Dtur_0933|PA5205|XCC3501|STM0356|ybbP|STM3549|yhjG|yjcH|yjiN|STM1874|yeeA|yegH|yhfK|STM0566|ybjE|STM1637|PA1439|yebZ|yfdC|yidG|yigC|yiaB|STM3170|yedI|ydcF|STM1158|STM1253|Caur_1083|Caur_0932|Caur_0931|THEYE_A0031|blr2954|Caur_3638|Caur_3637|Caur_1926|yhgN|hofC|yphA|yohC|ymgF|yhaI|yqjE|ybjO|ycfZ|ydjM|yfeZ|ycfT|yhiM|ycdZ|ylaC|ybhQ|cbrB|ybjJ|yhiD|yibH|yjgN|ytfF|igaA|yeeE|ybjM|yciC|HI_1064|eamA|eptA|HP_1343|PA5469|ygdQ|cptA|yvgT|BT_1201|BT_1813|HI_1005|yjbE|HI_0056|HI_1240|PA1037|ybiP|PA5250|PA4517|gll2754|BC_5423|opgE|SCO5481|XCC3770|DR_1187|DR_B0131|DR_2368|SCO4104|YPO0042|yhjW|NMB2010|eptB|blr3770|bll4584|RB12886|YPO2377|bll5340|PA3310|Caur_2247|PA1972|NMB1638|RB5758|BC_3709|terC|yicG|yhbX|eptC|yedA|yicL|yedE|yeiH|yidI|ymfA|pmrR|HP_0565|ygaZ|plsY1|plsY|PA0451|BC_1838|THEYE_A0592|yaaH|ybbJ|yjiG|Dtur_0511|Dtur_0336|Dtur_0766|YPO0274|CBO0870|CBO1581|Rv1283c|BT_1742|YPO0252|HI_0219|HI_1628|spr1140|LA_2073|CBO1577|CBO0866|CBO0634|cstA|spmB|CBO3310|SCO7517|FN2099|FN0221|Rv0473|BT_3284|Rv0870c|XCC2359|XCC2878|lmo1210|SAOUHSC_00621|LA_0353|Dtur_0338|CBO0538|CBO3589|yqeZ|XCC4058|TM_0865|YPO1785|lmo2062|YPO3784|lmo1211|ycnJ|PA3747|PA3235|PA3320|PA5478|PA5211|Caur_3065|PA4815|BC_5416|BC_1726|BC_2334|BC_2336|plsY2|BC_3542|gll2175|BC_1471|BC_2056|oppB|yjeT|elaB|rclC|yhcB|ygaM|XCC1193|SCO3964|SCO1797|yjiY|TM_0982|Caur_3272|SCO1826|DR_0790|DR_1112|DR_1113|DR_1457|DR_A0162|DR_A0298|DR_2142|SCO4610|SCO7226|NMB0571|SAOUHSC_00131|SAOUHSC_01677|SAOUHSC_02272|VC_0687|SCO4053|corE|THEYE_A1959|THEYE_A1590|THEYE_A1507|blr4114|bll3505|blr3521|blr3523|YPO3565|Rv1487|slr0935|RB7673|sll1678|RB2784|RB2785|RB3696|BT_4452|BT_4608|bll7749|bll4878|bll4952|bll6527|bll6908|PA4441|PA4606|PA3631|PA5430|BT_4556|BT_1541|nfeD|blr2899|blr2888|BC_3451|BC_3450|yqgA|satP|yqjD|yhhL|yejB|btsT|ypjD|lolC|yqiJ|Dtur_0119|GSU0394|GSU0830|GSU2135|Dtur_1765|Caur_0121|ybhN|yjfL|yijD|ygdD|ydgK|yidH|yiaA|yneE|ctrB|Caur_2575|yeaL|ycjF|yobD|yecN|yaiY|lolE|gfcA|yfjS|yebE|yjeO|yjdF|yfjD|yafY|yqiK|yeeR|ygbE|ybcI|ygiZ|CBU_0933|BC_3075|PA0659|PA0563|Rv3000|YPO2511|HI_1701|sanA|HI_1643|LA_3345|yubF|BT_1208|BT_0567|spr0034|LA_2113|SO_3801|SO_4029|CBO0651|CBO0363|CBO2900|CBO2822|yadH|XCC3978|YPO0112|yjkA|XCC0380|XCC1165|XCC1692|SAOUHSC_02638|FN0533|DR_0210|DR_0187|PA2811|ybbM|lmo2147|lmo2442|lmo0668|lmo0672|lmo1004|YPO1307|yajR|ycbC|PA1036|glr0923|BC_1617|BC_2231|BC_3183|BC_4795|BC_4941|BC_0349|BC_2001|BT_1919|TM_0193|SCO4629|STM0098|ygjQ|TM_1674|SCO1700|DR_1094|YPO0087|SAOUHSC_00322|NMB0004|NMB1979|SAOUHSC_02620|SAOUHSC_02753|SCO2472|VC_0590|aq_1986|aq_1210|THEYE_A0462|THEYE_A0244|THEYE_A1593|sll1606|blr3189|slr1647|RB9488|bsl7612|blr8071|blr8105|bll4980|bll5861|blr6578|PA3985|PA4233|PA5383|PA2066|BT_1337|BT_0879|BT_4609|blr1392|RB1526|BC_5174|gll4063|LA_0970|elyC|yabI|dedA|yjeM|fetB|yqjA|yghB|sfbC|Dtur_0563|GSU2664|Dtur_1052|ttpC|exbB|GSU1480|GSU1332|cysW|Caur_3239|Caur_3238|Caur_2538|Caur_0426|Caur_0324|ycaI|yagU|cysU|secA|potH|secA2|secA1|potB|ydcU|ydcV|GSU2782|GSU3400|GSU1611|exbD|oppC|glr3481|wzt|nagE|secG|secY|secE/);

/* outer membrane
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOouterMembrane = new Test("outer membrane (GO)",/ydeT|yqiG|ybgQ|slp|yehB|yraJ|yhcD|yfcU|yeaY|elfC|tolC|fimD|htrE|sfmD|chiP|pldA|ompX|ompW|slyB|btuB|qseG|ompN|yiaT|ompF|phoE|lolB|yaiO|ybhC|yncD|fadL|bamE|pal|ompA|yjfN|rclB|yjfY|yqhH|yhcN|yahO|bglH|bhsA|yehP|yehM|mdtP|lamB|cirA|ybiJ|mcbA|fepA|lolA|bsmA|nanC|lpp|cusC|ompL|fecA|nfrA|bcsC|yhjY|yaiW|yfeN|fiu|bamA|wza|bamC|bamD|bamB|lpoA|lpoB|fhuA|uidC|rzoR|dolP|mlaA|fhuE|pgaA|ypjA|ynfB|yfgH|yiaD|yfaZ|appX|flu|ychO|yzcX|eaeH|pqiC|yceK|mliC|yceI|ymgD|yobA|ecnB|ecnA|osmB|yeeJ|yghG|blc|ygiB|ygeR|yddL|yddB|tamA|lolC|lptE|lptD|yliI|marB|yncE|osmE|csgG|lptA|gltF|ycfJ|malM|yphF|nlpI|rzoD|acrZ|mepS|rsxG|digH|rcsF|mltA|csgF|csgE|mltB|pgaB|gspD|ompT|flgH|gfcE|ompG|flgG|yjbF|nlpE|ybaE|potF|sfmC|dsbG|potD|osmY|ybgP|hisJ|ycbF|livJ|yraI|asr|pqiB|btuF|yhcA|efeO|spy|pdeC|pspE|pdeG|yiaO|cydX|pdeN|elfD|sanA|nanM|ytfQ|yadV|yfcS|alsB|pdeB|ydhX|ybcL|ydeN|nlpD|rlpA|tsx|mipA|lptC|surA|csgB|emtA|mltF|pagP|fimC|letB|pdeD|cpxP|mlaC|aphA|hofC|ivy|ycgV|livK|ppdB|ynfF|zraP|rbsB|ugpB|sbp|yqiH|mreC|mltC|ppk|ydcS|flgI|gspM|yghD|lnt|ppiA|ptrA|dcp|hdeB|hdeA|glnH|gspF|xylF|mglB|gspH|fliF|nrfB|rcnB|bglX|mepK|yejA|bisC|araF|oppA|thiB|lolE|mdtN|yjcS|artI|ssuA|iap|rseB|tamB|ampC|gspJ|gspI|gspC|sapA|pqqL|ftsP|fkpA|ccmG|gspG|cpoB|skp|loiP|ybaV|bepA|lolD|amiD|ompC|yfaL|hofQ|tesA|dmsA|ygiS|cusF|proX|msrP|gcd|fepB|dsbC|dsbA|degS|torZ|ygiW|cueO|nikA|ddpA|hyaB|amiB|yghE|pstS|malS|nrfA|osmF|opgB|metQ|prc|chiA|envC|degP|slt|torA|opgG|gspL|fhuD|opgD|phoA|gltI|fecB|tauA|degQ|ynfE|mppA|gsiB|agp|gspE|amiC|pgpB|cusB|zinT|ggt|fdnG|malE|tcyJ|envZ|tolB|modA|cpxA|acrA|dacC|amiA|cpdB|ldtB|appA|cysP|dppA|eco|efeB|rna|argT|napA|lpxT|sodC|artJ|evgS|emrA|mepA|treA|torC|kdpD|ansB|ushA|mltG|narX|mrdA|acrB|dcuS|dacA|glrK|lysP|tonB|dpiB|cusS|sgrR|mrcB|rcsC/);

/* peptidoglycan
+   taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOpeptidoglycan = new Test("peptidoglycan (GO)",/glmM|mraY|murA|rlpA|pbpG|oppC|oppB|metN|nlpI|nagE|mltA|pgpB|pal|ompA|mipA|mltD|mrdA|mltC|murJ|motB|dacA|ampG|ldtE|mepH|rrrQ|dacC|mpaA|ampH|ddlA|ftsI|torY|ynbB|oppD|ldcA|glmU|amiB|yghE|amiA|rfe|ddlB|ftsX|ampD|mrdB|ftsW|murI|murD|bacA|murF|murB|ldtB|murC|mepK|elyC|mltB|emtA|oppF|mpl|mtgA|lpoA|mrcB|mrcA|ycjG|chiA|pbpC|murQ|envC|yiaT|dedD|nagK|ldtC|nagZ|flgJ|ldtD|lpoB|yieL|yegX|mltF|slt|amiD|ybjG|mepS|mepM|dacB|lpxM|mepA|torC|ytfB|anmK|damX|gspL|lpp|ldtA|mltG|ftsN|murE|ispU|mppA|dacD|murG|ldtF|amiC|rrrD/);

/* plasma membrane
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOplasmaMembrane = new Test("plasma membrane (GO)",/ysaB|hemX|yajG|yghQ|essQ|yahG|sieB|yahC|yicS|ybhN|ybhI|chiQ|yghJ|marC|yohP|yjfL|yihY|yjeT|yjcH|yijD|yiiR|yoaI|yedD|yjiG|ynjI|ynjF|ynjC|ydjZ|ydjX|yniB|yhhN|yhfL|hokA|ygdD|ydcF|ykgR|yifL|ydgK|yidX|yidG|ynfA|yiaW|yeaL|ydgC|yphA|yohK|ppdA|yfjW|ycjF|yfeH|ypdK|azuC|yoaJ|yoaK|rzoQ|ythA|ymgE|ygdR|ygdI|hemY|fxsA|yqgA|yhdV|yhdU|ychE|ynaM|ymiC|ykiC|yahV|ynfS|ydjK|yobD|ppdC|yhaI|yhaH|creD|yqjE|ybjO|yniD|ypfN|yecN|ybbJ|ybaN|ycfZ|yceB|yaiY|yccF|ybhM|ybhL|ymcE|yebO|yeaQ|ydhI|yraQ|ykgH|ychQ|yeaY|yeiH|yibI|dmsD|yfeZ|yfdY|yfhR|yfdI|ymfR|ymfE|bcsB|yhiI|yhhM|yhhL|ycdZ|ycdU|yfiB|yaiZ|ylaC|ybfB|ybhQ|yicG|cbrB|atoE|yidI|yicL|yicH|ytiC|gfcA|gfcB|ytjB|yfjS|ybjT|dedA|yncO|yhiD|ydiK|yadS|ynfC|ydgU|yncL|ybjP|mokC|yebE|yjiN|yjiK|yjiJ|yjiH|yjgN|ytfF|yhdP|yciQ|yhgE|yjeO|yjdF|gatC|yohJ|ybdJ|ybbP|rarD|yigG|yifK|ydbH|yeeE|yqhA|yqjF|ygjV|gspB|yeiB|yigF|yhhH|yohO|ygjQ|ygjI|ygcG|djlB|yafY|ygdQ|yqeJ|yhgN|ytjA|yeeR|yafT|ygbE|ymiA|yedI|yajI|ybcI|yhfT|ygiZ|ydcL|ynaJ|pmrR|ybjM|ypjD|ypfJ|yhbE|yciC|cysW|djlA|hycD|yhcB|yedE|hokD|hflD|ybaY|yihG|hokB|yqfA|ldrD|yfdH|agaD|agaC|mdtL|sgcC|srlA|hybO|hybC|mcbA|manZ|manY|ghxQ|mdtG|ghxP|agaS|emrD|ydhC|gdx|ldrB|hypA|ibsD|pqiA|yfgO|ybiR|nuoI|syd|hyfE|hycF|hdeD|ygiM|yqaE|ampE|yidH|ibsE|ibsC|shoB|ibsB|ibsA|hokC|yibN|ghoT|ldrC|ldrA|ybfA|ycjP|yagU|ygeQ|ydcU|yiaV|yhhT|ycjO|ybhS|ybhR|ydhU|nuoH|rseC|ortT|yedA|yfdV|lsrB|potH|hisM|yadH|yfgM|phnD|yhhJ|ycjV|yqiK|artQ|artM|ptsG|dmsB|chbC|glvC|fryC|srlE|glvB|yddG|ccmG|safA|osmE|rsxA|cvpA|yccM|yaaU|yiaB|yqiJ|araJ|ydjE|arnC|nuoJ|nuoE|nuoA|abrB|uspB|yfcA|yebZ|yebQ|nimT|yoaE|ybeX|citT|paeA|ygaM|ydiN|ydiM|ydhK|ydhJ|rsxD|yiaA|nuoL|yegH|yohC|acrF|mdtD|ymgF|yeaV|trkA|yqjD|gltF|yohD|ynfM|ybhG|yajR|ttdT|ydjM|yihN|yiaD|ycfT|yhjE|yhjD|yhhS|tap|ycfJ|mdtH|hokE|yqcE|yfbL|hsrA|napH|ccmA|yfaV|eamB|yccS|ftsE|ybjJ|yibH|appX|lsrD|pdeA|lsrC|ansP|flk|yfcC|yfcJ|rsxE|yhdX|yeeA|hyfH|yfjD|ydhP|ygcS|yphD|narH|ycaD|ygaP|yhfK|yqjA|qmcA|dinQ|lptG|potI|yadI|lptF|gltK|gltJ|frwD|hflK|yjfF|sgcA|pldB|yejE|yejB|ytfT|cmtB|yhdY|znuB|fruB|manX|crr|hisQ|frdD|frdC|livM|frdB|yjcE|frdA|sdhA|fruA|zntB|dmsA|treB|murP|emrY|bglF|ccmH|alaE|ccmE|emrK|uraA|ybhF|macA|ascF|ulaA|cmtA|hyaA|narI|mtlA|kbaZ|livH|lolC|btuC|srlB|modB|mlaE|ppiD|yggE|araH|tatE|ftsA|exbD|sgcB|lgt|xylH|yehY|yehW|yejA|alsC|ptsA|lolE|mqo|rbsC|sapC|agaV|chbB|chbA|sapA|ddpC|ddpB|fetB|eamA|ydeA|ydcZ|narG|yjjP|ecnB|ecnA|osmB|ftsB|yciB|era|ygfX|pspA|hemG|ybgC|nuoN|nuoM|nuoK|nuoB|nfrB|acrE|yfbV|ftsL|elaB|rclC|yagG|yajC|yqaA|dld|napG|fsr|nuoG|yneE|nuoC|yejM|letA|ynbB|ynbA|kch|gndA|nlpC|pspD|ftsX|ycaI|dcuA|pheP|mreD|hyfC|yaeF|yfgG|ftp|yidE|ybjL|ydgI|drpB|nlpA|nrfD|pdeC|arnE|yqeG|yabI|yihP|pdeG|yiaN|yiaM|yhjV|dedD|cydX|pdeK|yhiM|pgaC|xanP|adeQ|rsxC|pdeN|yicJ|ycaM|lfhA|yafL|sanA|secF|yfbS|ydcV|sohB|lgoT|pspG|tsgA|xanQ|aaeA|aaeB|uidB|lpxH|bsmA|ymfA|nuoF|damX|gspK|pdeB|hyfF|alx|rsxB|dinF|wcaD|ybbY|yceO|yghB|livG|feoB|livF|artP|cydC|mdtM|cydD|bcr|ygiS|proX|potC|potB|gsiD|gsiC|dppC|dppB|degS|sgrT|nikA|gatB|ptsP|ddpA|lptB|mglC|hflC|lapB|malG|osmF|ydcT|fdoH|fdnH|sapF|metQ|malF|cysU|tatA|tatB|tatC|frwB|thiP|proW|yejF|sapB|ssuC|bcsQ|ddpF|ulaC|fryB|ugpE|ugpA|fryA|ddpD|gutM|fecE|fecC|fecD|tauA|degQ|mppA|gsiB|ydeE|yjjB|dmsC|dsbD|pppA|pnuC|trg|yodB|proP|glpG|oppC|oppB|nlpI|lapA|yidD|letB|pdeD|yahN|rhaT|ftsQ|gspO|psiE|ynfH|nepI|pspC|hofC|hybB|mlaD|pitB|lysP|napC|atpI|ppdB|yggT|proY|tyrP|mdtO|kgtP|zipA|pgaD|htpX|yhjG|rodZ|yceJ|acrZ|sppA|rmuC|tqsA|cydH|yhjX|yjhF|nanX|lepB|igaA|fliM|yjeM|aer|gspA|psuT|nupX|ybaT|cynX|fliN|mreC|djlC|asmA|ygbN|rsxG|fdhF|nagE|mdfA|emrB|mdtK|panF|frvB|emrE|malX|artJ|macB|emrA|mdtI|frwC|mngA|mdtJ|gatZ|hisP|dacA|tonB|rutG|proQ|ydcO|gspM|yghD|uacT|lnt|ygaZ|pepN|msrQ|ybiO|zitB|ybbW|yadG|flhB|flhA|gspF|yfdC|zapA|araE|bcsF|dgcE|mdtC|mdtB|hofB|dsbB|torY|ynbD|ybaL|gspH|tolR|fliF|yijE|sdhD|secD|ycaL|ygaH|fliP|yidC|cydB|cydA|cdsA|elyC|hyaC|nrfC|fliO|yccA|nrfG|nrfF|wbbK|wbbH|yihO|yjbF|bcsG|ghrA|dgcT|gadC|nrfE|mdtN|ubiA|yegT|ypdI|tolQ|fliL|mgtS|plsY|tamB|ytfB|yfeW|efeU|lspA|mscM|tehA|gspJ|gspI|wzyE|gspC|lplT|mzrA|hyfD|ftsN|fliR|amtB|kefB|gspG|potG|ytfR|ftsY|glnP|btuD|cusA|parC|gsiA|nikC|minD|malE|dppF|rbsA|nikB|dhaM|modF|exbB|cyoA|malK|ssuB|sapD|mglA|dppD|araG|actP|frvA|fepG|fepD|dppA|yfhM|thiQ|lolD|cysA|proV|yojI|znuC|tcyL|btsT|msbA|fhuB|potA|ulaB|ugpC|lpxB|modC|alsA|cstA|tauC|fetA|narZ|nlpD|mscL|minE|rlpA|mdlA|ydfJ|psd|wecF|ppx|yeeO|shiA|nupC|motA|dtpD|nhaB|mltA|cdh|mltD|plsC|acrD|hycC|zupT|yjbB|gcd|narV|glpF|galP|motB|mntP|leuE|dgcJ|cyoE|csgG|lptC|ampG|bioP|plsB|hyfB|glcE|dgcF|mdtA|pspB|arsB|appC|appB|glcF|clsA|pgpC|ccp|hyaB|yfeO|cobS|fliG|yghE|gntU|hofN|nirC|focA|dcuC|dcuB|ccmD|ccmC|ccmB|corA|pgsA|crcB|mgrB|ybdG|mdlB|oppF|focB|aqpZ|sdaC|opgB|wzxE|narU|rfbX|lacY|mlaF|clsB|menA|prc|pssA|gltP|mtr|envC|wecH|bcsA|mdtF|mdtE|opgC|sdhB|phnC|fdrA|degP|adeP|ccmF|lit|pstC|rcnA|xylE|uhpT|baeS|essD|fliQ|hcaT|sbmA|aroP|rhtC|opgH|mscS|argO|mhpT|pgpA|gspL|hofO|fliJ|narK|yddA|chaA|lsrA|waaL|uhpB|cheR|wcaJ|wzxC|ynfE|yqeI|hlyE|gudP|mscK|yhbX|eptC|xapB|gspE|garP|exuT|pstA|plaP|tisB|dnaA|nikE|nikD|fhuC|afuC|yehX|tcyN|mrcA|fepC|metI|rbbA|srlR|tauB|metN|phoU|atpB|ptsN|gatA|fadD|aas|htpG|mraY|mntH|kdgT|yhdZ|umuC|glcA|atoS|cusB|setA|nhaA|wzzB|pitA|nupG|glpC|dtpB|lpxL|gltS|murJ|opgE|fdoI|fdnI|envZ|dgcP|cdgI|xylG|rcsD|lhgD|entS|cpxA|ribB|acrA|clcB|qseC|wzzE|glpA|mreB|brnQ|dacC|lpxP|ampH|dadA|ftsI|btsS|ubiB|puuP|cvrA|sdhC|rhtB|satP|gntP|creC|rfe|waaA|secG|dgkA|glpD|glpB|cyoD|cyoC|cyoB|mrdB|ftsW|cybB|bacA|betT|mgtA|clcA|tnaB|hybA|nanT|ftsH|pstB|torS|potE|cycA|tdcC|ubiD|fhuF|mtgA|melB|dcuD|ygeY|putP|pdeF|trkG|pbpC|caiT|eutH|ghrB|yrbG|eptB|yhhQ|zntA|poxB|pntA|ldtD|fieF|pntB|arnF|arnT|yidK|lpxT|csrD|fadK|setC|sstT|evgS|eptA|basS|secY|kup|lysO|ybjG|dgcI|lpxK|secA|trkH|kdpF|secE|rng|cheZ|mepM|dauA|lpxM|mnmE|lldD|lldP|torC|yfdG|idnT|yaaJ|rstB|kdpC|kdpB|glpT|kdpA|kdpD|dtpC|tcyP|fucP|kefC|mmuP|setB|cyuP|ndh|dsdX|mltG|uhpC|gntT|cheB|adiC|dtpA|dgcN|abgT|tolA|entD|betA|murG|frlA|zraS|pyrS|dosC|codB|rhtA|atpG|putA|ptsH|entB|dgcZ|kefF|copA|hprS|cheW|narX|glnQ|pgpB|mrdA|rseP|acrB|dcuS|barA|arcB|dctA|pyrD|glrK|ppk|wzc|cysZ|uup|fepE|gltL|cadB|dgcQ|tsr|phoQ|tar|etk|rseA|dacB|gabP|hda|tufB|tufA|yjeH|dpiB|narQ|cusS|cysQ|entF|dgoT|ftsK|nadR|fecR|csgD|atpC|phoR|atpH|atpF|dgcC|mrcB|cadC|rpoH|lpd|gpt|cheA|atpE|lepA|rcsC|dgcM|dnaK|atpD|atpA/);

/* fatty acid
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOfattyacid = new Test("fatty acid (GO)",/atoE|cfa|fabY|aas|menE|caiC|fadD|fadL|fabF|gnsB|ydiO|yqeF|fabB|ygcQ|caiA|ydiR|yciA|scpB|ydiF|fadM|prpB|paaG|atoA|atoD|atoB|fadK|fadH|fadB|fadA|fixB|fadI|fadE|fabH|fabD|oxc|paaH|paaF|caiD|paaJ|prpE|fadJ|fabR|prpD|prpC|acpS|accB|scpA|acpH|acnB|fabI|fabG|fabZ|fabA|fadR|accC|ackA|accA|lpxL|scpC|lpxP|acpP|waaA|acpT|tdcF|tesB|glnB|lpxM|lipB|lipA|tdcE|tdcG|tesA|plsB|accD|prpR|tdcD|tdcB|tdcA|pta|tcyJ|gcvH|acrA|putP|tcyL|hmp|acrB|tcyN|sucB|aceF/);

/* lipoprotein
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOlipoprotein = new Test("lipoprotein (GO)",/ysaB|lolC|lptE|ecnB|ecnA|osmE|osmB|yajG|rlpA|lnt|nlpI|chiQ|gfcE|yghJ|pal|slyB|yedD|yhfL|acrA|yifL|yidX|qseG|rzoQ|slp|yqhH|ygdR|ygdI|yfgH|nlpC|yhdV|rzoR|rzoD|lgt|dolP|ftp|yceB|pqiC|yehP|yehM|nlpA|ybaY|metQ|lpoA|ygeQ|yeaY|rcsF|yjbF|mlaA|yiaD|lolE|lolD|yfiB|phnL|lpoB|gfcB|ybbA|yfjS|ypdI|ynfC|ybjP|lolB|lolA|lspA|yaiW|lpp|cusC|nlpE|yafY|borD|yafT|mliC|ybhC|yajI|ydcL/);

/* antiporter
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOantiporter = new Test("antiporter (GO)",/nhaA|nhaB|mdfA|mdtK|citT|clcB|ybaL|cvrA|dcuC|dcuB|dcuA|clcA|potE|cadB|ydgI|ttdT|yjcE|emrE|caiT|yrbG|gadC|uhpT|yfbS|mdtI|mdtM|glpT|yfcC|kefC|narK|chaA|dinF|uhpC|adiC|mdtJ|bcr|kefB/);

/* exporter
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOtransporter = new Test("transporter (GO)",/lsrB|ybaE|potH|potG|potF|ydeT|hisQ|ydeE|eamA|ydeA|yncD|ydcS|ydcO|lptG|lolC|lptE|mntH|kefG|kefF|mscL|ytfR|treB|kdgT|copA|ygiS|murP|glnP|yhdZ|btuD|mdlA|btuC|yghQ|glcA|cusB|uacT|ghxQ|ydfJ|emrY|yaaU|setA|nhaA|yqiG|araJ|proX|sbcC|potI|potD|potC|potB|pnuC|srlB|ydjE|cusA|yeeO|shiA|pitA|yfgO|proP|oppC|oppB|gsiD|gsiC|gsiA|nupG|nupC|nuoN|ybiR|ybiO|fiu|nuoM|nuoK|nuoJ|nuoI|fadL|nuoE|ybhI|motA|nuoB|nuoA|metN|zitB|ybgQ|dtpD|nagE|nikC|nhaB|chiP|glnQ|gfcE|wza|ybbW|tsx|yadG|ompA|mdfA|yfcA|bglF|dtpB|acrE|malE|lptF|livH|acrD|hycF|hycE|hycC|hisM|hisJ|zupT|yjbB|gltS|gltK|gltJ|glpF|glnH|murJ|galP|motB|modB|tcyJ|fepB|yebQ|mntP|emrB|leuE|nimT|acrB|ynjC|dppC|yahN|dppB|xylG|xylF|yagG|mglB|znuB|mlaE|cysW|mdtK|csgG|modA|yfdC|rhaT|dppF|dctA|btuB|csgF|csgE|entS|lptC|citT|lptA|artQ|artM|araH|araE|ampG|bioP|acrA|ydiN|ydiM|feoB|ydhK|ydhJ|clcB|tatE|nepI|hyfB|rbsA|fsr|nuoL|nuoG|mdtC|mdtB|mdtA|yadH|arsB|appC|appB|panF|brnQ|livJ|sgrT|cysZ|nuoC|sgrR|nikE|nikD|nikB|nikA|atpC|yohK|ybaL|yedE|modF|ompG|ycjN|puuP|oppD|ycgV|gatB|cvrA|ddpA|livK|acrF|livF|livM|yfeO|frwD|mdtD|rhtB|satP|gntU|gntP|tolR|exbD|exbB|yijE|mlaD|mlaB|pstS|secD|secG|pitB|nirC|focA|fhuC|mdtG|dcuC|lysP|dcuB|dcuA|ccmC|ccmB|fecA|sgcB|cyoD|cyoC|cyoA|cyoB|lptB|corA|mglC|mrdB|ftsW|ydjK|yeaV|betT|xylH|trkA|mgtA|atpD|atpA|atpG|atpH|atpF|clcA|afuC|crcB|fepE|ghxP|alaE|bglH|ybdG|malK|malG|tnaB|osmF|yehY|yehX|yehW|yehB|ynfM|ydcT|nanT|ssuB|yejA|sapF|sapD|pstB|yidE|mglA|mdlB|gltL|dppD|ybjL|artP|araG|potE|ompN|cadB|ydgI|proY|cycA|oppF|focB|aqpZ|tdcC|sdaC|tyrP|yraJ|alsC|yajR|ttdT|mdtO|mdtP|ycjP|actP|yjcE|wzxE|metQ|kgtP|tcyN|yjfF|narU|lamB|rfbX|tonB|rbsB|araF|melB|lacY|mlaF|malF|chbC|arnE|cirA|dcuD|yhcD|yqeG|putP|hisP|emrK|yfcU|cysU|cysP|pgaA|frvA|frvB|tatA|tatB|ydcU|tatC|emrE|cydC|yihP|yihO|yihN|fepC|fepG|fepD|frwB|trkG|dppA|oppA|lptD|caiT|thiB|eutH|thiP|thiQ|metI|gltP|mtr|yiaN|yiaM|argT|yrbG|yhjV|lolE|lolD|yhjE|yhjD|mdtF|mdtE|rbbA|yhhS|yhhQ|zntA|phnD|gadC|phnC|cysA|fieF|proW|proV|mdtH|acrZ|mdtN|yqcE|yojI|xanP|uraA|yhhT|fepA|hsrA|ravA|ccmA|arnF|yfaV|adeP|ccmF|mdtL|atoE|yejF|yejE|glvC|yidK|emrD|adeQ|rbsC|yicL|setC|yicJ|pstC|sapC|sapB|yhhJ|agaV|rcnA|ydhC|rutG|yegT|xylE|eamB|sstT|znuC|yccS|malX|artJ|artI|uhpT|elfC|ssuA|ssuC|ybhF|secY|ybbA|kup|ycaM|macB|macA|livG|lysO|ftsE|ybjJ|sthA|emrA|secA|trkH|hcaT|sbmA|kdpF|secE|ugpB|phoU|sbp|yfbS|yejB|aroP|ascF|tcyL|ydiK|tqsA|lsrD|ydcV|rhtC|ycjO|dauA|atpE|mdtI|ybhS|ybhR|mscS|gabP|tolC|ompF|phoE|yhjX|lldP|lgoT|btsT|mdtM|yjiJ|frwC|lsrC|gdx|hyfI|tsgA|yjhF|nanX|idnT|msbA|xanQ|ytfT|ytfQ|ddpF|ansP|ulaA|yaaJ|kdpC|kdpB|atpB|fimD|glpT|kdpA|aaeA|aaeB|argO|uidB|ompC|efeU|fhuB|fhuA|htrE|potA|uidC|mscM|yjeM|mhpT|nanC|yjeH|dtpC|fryC|alsB|yfcC|ptsN|gatA|cmtA|ulaB|ulaC|yfcJ|fruB|fryB|manZ|manY|tcyP|yddL|cydD|yphE|ybbP|nuoF|srlE|fucP|kefC|ynjD|tehA|rarD|yhdY|yhdX|yifK|mmuP|fliI|ugpC|ugpE|ugpA|narK|manX|chbB|chbA|glvB|ptsG|sapA|setB|psuT|nupX|lplT|sufC|yeeA|ycjV|znuA|sfmD|ddpC|cyuP|fryA|hyfH|ndh|fhuD|dsdX|ybaT|nuoH|cusC|yddA|chaA|lsrA|dinF|ddpD|yphF|cynX|uhpC|modC|alsA|cstA|gntT|gltI|fecE|fecC|fecD|fecB|tauC|tauB|tauA|amtB|ygjI|mngA|ydhP|wzxC|adiC|mppA|mdtJ|hyfG|gudP|ygcS|yphD|ddpB|fetB|dtpA|ompL|gsiB|fetA|bcr|ybbY|mscK|yddG|abgT|ycaD|tolA|fruA|xapB|mtlA|frlA|yhfK|kefB|yfaL|fhuE|ygbN|codB|garP|exuT|dgoT|yhbE|yedA|zntB|rhtA|yqjA|yghB|fdhF|pstA|plaP|tisB/);

/* transporter
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOexporter = new Test("exporter (GO)",/ydeE|eamA|ydeA|yjjP|yjjB|ygaZ|zitB|mntP|leuE|nimT|yahN|entS|paeA|ydiM|nepI|rhtB|yijE|ygaH|ccmD|ccmC|ccmB|alaE|sapF|sapD|cydC|fieF|setC|sapC|sapB|rcnA|eamB|ybhF|lysO|rhtC|ybhS|gdx|argO|yjeH|cydD|kefC|fetB|fetA|yddG|kefB|rhtA/);

/* efflux
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOefflux = new Test("efflux (GO)",/ydeE|eamA|ydeA|kefG|kefF|cusB|emrY|setA|araJ|zitB|mdfA|acrE|acrD|mntP|emrB|acrB|entS|acrA|ydhK|ydhJ|nepI|fsr|mdtC|mdtB|mdtA|acrF|mdtD|rhtB|mdtG|alaE|ybhG|tdcC|mdtO|mdtP|narU|arnE|emrK|emrE|yiaV|mdtF|mdtE|fieF|mdtH|acrZ|mdtN|arnF|mdtL|emrD|setC|rcnA|eamB|yccS|ybhF|macB|macA|lysO|emrA|yibH|tqsA|rhtC|mdtI|ybhS|ybhR|tolC|mdtM|gdx|aaeA|aaeB|kefC|tehA|narK|setB|lplT|yeeA|cusC|ydhP|mdtJ|bcr|yhfK|kefB|zntB|rhtA|yfdV/);

/* tRNA ligase
+	taxon_subset_closure_label: Bacteria	
+	taxon_subset_closure_label: Escherichia coli K-12
*/
GOtrnaLigase = new Test("tRNA ligase (GO)",/epmA|lysU|lysS|thrS|serS|ileS|gltX|leuS|hisS|argS|pheS|valS|tyrS|gluQ|cysS|tcdA|glnS|glyS|metG|alaS|proS|aspS|trpS|asnS|pheT|glyQ|tilS|rtcB/);





// ========= EXPERIMENT STARTS HERE 

function getProximityReport(geneIndexes, proximity, filtered) {

    let PROXIMITY = proximity;
    let FILTERED = filtered;
    let moonies = geneIndexes;

    // FILTERED means remove hypotheticals, as follows
    if (FILTERED) {
        finalArray = filterPA(getNeighbors(moonies, PROXIMITY), /hypothetical/);

    } else
        finalArray = getNeighbors(moonies, PROXIMITY);

    // don't count ourselves
    finalArray = filterPA(finalArray, moonlit);

    console.log(finalArray.length + " genes");

    console.log("PROXIMITY: " + PROXIMITY);
    console.log("FILTERED: " + FILTERED);
    console.log("Fold\tE\t\tFunction");
    console.log(enrichments(finalArray, [].concat(hypothetical).concat(GOcellWallOrganization).concat(GOtransmembrane).concat(GOcellDivision).concat(GO_permease).concat(GOsecretions).concat(GOinnerMembrane).concat(GOplasmaMembrane).concat(GOpeptidoglycan).concat(GOfattyacid).concat(GOlipoprotein).concat(GOtransporter).concat(GOantiporter).concat(GOexporter).concat(GOefflux).concat(GOtrnaLigase).concat(GOtransmembraneTransport).concat(GOouterMembrane)))
}

// RUN IT:

// Get gene indexes of all moonlighting genes as 'moonies' array.
moonies = harvestGeneIndexesUsingPattern(moonlit);

getProximityReport( moonies, 5, false );
