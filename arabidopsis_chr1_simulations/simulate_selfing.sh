#! /bin/bash

model="model0"

mkdir ${model}
cd ${model}

h="selfing"
for selfingpercent in "0"  "25" "50" "75"
do
mkdir h_${h}; cd h_${h}
s=1000

scalingfactor=100

if [ "${h}" == "0.0" ]
    then fitnessblock=""
    else 
read -r -d '' fitnessblock << endmsg 
1:13000 fitness(m1) {
    h = 1/((1/0.987) - 39547 * mut.selectionCoeff/${scalingfactor});
    //h = mut.mutationType.dominanceCoeff;
    if (homozygous) {
        return ((1.0 + 0.5*mut.selectionCoeff)*(1.0 + 0.5*mut.selectionCoeff));
    } else {
        return (1.0 + mut.selectionCoeff * h);
    }
}
endmsg
fi

cat > slim_h${h}_full_split${s}_self${selfingpercent}.slim << EOM
// use page 76 to randomly generate exons
initialize() {
    initializeMutationRate(0.7e-8*(2.31/3.31)*${scalingfactor}); // no synonymous, lowered mutation rate by factor 2.31/3.31
    
    //nonsynonymous drawn from a DFE Huber arabidopsis paper
    // importantly, SLiM computes the fitness of the heterozygote and homozygote as 1+sh and 1+s
    // dadi and others compute it as 1+2sh and 1+2s
    // to convert from dadi to SLiM, E[s] = -shape*scale
    // scale up by a factor of ${scalingfactor}
    // E[s], shape
    initializeMutationType("m1", 0.5, "g", -0.00048655*${scalingfactor}, 0.185);
    // neutral
    // initializeMutationType("m2", ${h}, "f", 0.0);
    
    // p1 marker mutations
    initializeMutationType("m10", 0.5, "f", 0.0);
    m10.convertToSubstitution = F; //marker mutations are always tracked
    
    //genomic element: exon and uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
    //initializeGenomicElementType("g1", c(m2,m1), c(1.0,2.31));
    //genomic element: exon and uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
    initializeGenomicElementType("g1", c(m1), c(1.0)); //with no synonymous sites
    //genomic element: neutral
    //initializeGenomicElementType("g2", c(m2), c(1.0));
    
    //read in exon and recomb info
    info_lines = readFile("/u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_arabidopsis/sim_seq_info.txt");
    
    //recombination
    rec_ends = NULL;
    rec_rates = NULL;
    for (line in info_lines[substr(info_lines, 0, 2) == "rec"])
    {
        components = strsplit(line, " ");
        rec_ends = c(rec_ends, asInteger(components[1]));
        rec_rates = c(rec_rates, asFloat(components[2]));
    }
    initializeRecombinationRate(rec_rates*${scalingfactor}, rec_ends); //multiply rec rates by scaling factor
    
    //neutral stuff
    //for (line in info_lines[substr(info_lines, 0, 2) == "neu"])
    //{
    //    components = strsplit(line, " ");
    //    neutral_starts = asInteger(components[1]);
    //    neutral_ends = asInteger(components[2]);
    //    initializeGenomicElement(g2, neutral_starts, neutral_ends);
    //}
    
    //exons
    for (line in info_lines[substr(info_lines, 0, 2) == "exo"]) 
    {
        components = strsplit(line, " ");
        exon_starts = asInteger(components[1]);
        exon_ends = asInteger(components[2]);
        initializeGenomicElement(g1, exon_starts, exon_ends);
    }
}

${fitnessblock}

//read burn in population
1 early() {
    defineConstant("simnum", getSeed());
    setSeed(getSeed() + $RANDOM);
    sim.addSubpop("p1",1000);
    // when set to "F" does not remove fixed variants
    // this isn't necessary because we're keeping track of all populations
    m1.convertToSubstitution = T;
}

10000 early() { // after burn-in, split populations into two: p1 and p2
    sim.addSubpopSplit("p2", 1000, p1);
    p1.setSubpopulationSize(1000);
    
    // this also isnt necessary, but sets the migration rates to 0
    p1.setMigrationRates(c(p2),c(0.0)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.0)); //migration rate INTO p2
    
    // make p2 a partially selfing population
    p2.setSelfingRate(0.${selfingpercent});
}

//add marker mutations every 1000 bp
$((10000 + ${s} - 1)) late() {
    for (pos in asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/1000))*1000)-1)) {
        p1.genomes.addNewDrawnMutation(m10, pos);
    }
}

$((10000 + ${s})) early() {
    p1.setMigrationRates(c(p2),c(0.05)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.05)); //migration rate INTO p2
}

$((10000 + ${s})) late() {
    p1.setMigrationRates(c(p2),c(0.0)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.0)); //migration rate INTO p2
    p2.setSubpopulationSize(1000);
    p1.setSubpopulationSize(1000);
}

$((10000 + ${s} + 1000)) { 
    p1.outputVCFSample(1000, replace=F, filePath="/u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_arabidopsis/burnin_${model}/h_${h}/h${h}_self${selfingpercent}_" + simnum + "_p1.vcf");
    p2.outputVCFSample(1000, replace=F, filePath="/u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_arabidopsis/burnin_${model}/h_${h}/h${h}_self${selfingpercent}_" + simnum + "_p2.vcf");
}

EOM

cd ..
done


cd ..