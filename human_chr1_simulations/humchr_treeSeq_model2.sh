#! /bin/bash

model="model2"
Nanc="10000"
Ns="1000"
Nr1="10000"
Nr2="10000"
splittime="20000"

scalingfactor="5"

mkdir burnin_${model}
cd burnin_${model}

for h in "0.0" "0.5" "s"
do
mkdir h_${h}; cd h_${h}

if [ "${h}" == "0.0" ]; then 
    fitnessblock=""
    hh="0.0"
elif [ "${h}" == "0.5" ]; then
read -r -d '' fitnessblock << endmsg 
1:$( echo "( ${Nanc} * 10 + ${splittime} + ${Nanc} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f") fitness(m1) {
    h = mut.mutationType.dominanceCoeff;
    if (homozygous) {
        return ((1.0 + 0.5*mut.selectionCoeff)*(1.0 + 0.5*mut.selectionCoeff));
    } else {
        return (1.0 + mut.selectionCoeff * h);
    }
}
endmsg
hh="0.5"
else 
read -r -d '' fitnessblock << endmsg 
1:$( echo "( ${Nanc} * 10 + ${splittime} + ${Nanc} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f") fitness(m1) {
    h = (0.5)/(1 - 7071.07*(mut.selectionCoeff/${scalingfactor}));
    //h = mut.mutationType.dominanceCoeff;
    if (homozygous) {
        return ((1.0 + 0.5*mut.selectionCoeff)*(1.0 + 0.5*mut.selectionCoeff));
    } else {
        return (1.0 + mut.selectionCoeff * h);
    }
}
endmsg
hh="0.5"
fi

cat > slim_h${h}_full_split${s}.slim << EOM
// use page 76 to randomly generate exons
initialize() {
    initializeMutationRate(1.5e-8*(2.31/3.31)*${scalingfactor});
    initializeTreeSeq();
    
    //nonsynonymous drawn from a DFE from Kim et al.
    // importantly, SLiM computes the fitness of the heterozygote and homozygote as 1+sh and 1+s
    // dadi and others compute it as 1+2sh and 1+2s
    // to convert from dadi to SLiM, E[s] = -shape*scale*2/(2Na) = -shape*scale*2/(2*12378)
    // scale up by a factor of ${scalingfactor}
    // E[s], shape
    initializeMutationType("m1", ${hh}, "g", -0.01026*${scalingfactor}, 0.186);
    //synonymous -- assumed neutral here
    initializeMutationType("m2", ${hh}, "f", 0.0);
    //noncoding -- assumed neutral here
    //initializeMutationType("m3", ${hh}, "f", 0.0);
    //beneficial -- for the time being, is left out
    //initializeMutationType("m4", ${hh}, "f", 0.0);
    
    //genomic element: exon and uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
    //initializeGenomicElementType("g1", c(m2,m1), c(1.0,2.31));
    initializeGenomicElementType("g1", c(m1), c(1.0)); // no synonymous mutations, mutation rate is scaled
    //genomic element: intron
    //initializeGenomicElementType("g2", c(m3), c(1.0));
    //genomic element: intergenic
    //initializeGenomicElementType("g3", c(m3), c(1.0));	
    
    //read in exon and recomb info
    info_lines = readFile("/u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/sim_seq_info.txt");
    
    //recombination
    rec_ends = NULL;
    rec_rates = NULL;
    for (line in info_lines[substr(info_lines, 0, 2) == "rec"])
    {
        components = strsplit(line, " ");
        rec_ends = c(rec_ends, asInteger(components[1]));
        rec_rates = c(rec_rates, asFloat(components[2]));
    }
    initializeRecombinationRate(0.5*(1-(1-2*rec_rates)^${scalingfactor}), rec_ends); //multiply rec rates by scaling factor
    
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
    
    //conserved non-coding is neutral
    //for (line in info_lines[substr(info_lines, 0, 2) == "cnc"])
    //{
    //    components = strsplit(line, " ");
    //    cnc_starts = asInteger(components[1]);
    //    cnc_ends = asInteger(components[2]);
    //    initializeGenomicElement(g3, cnc_starts, cnc_ends);
    //}
}

${fitnessblock}

//read burn in population
1 early() {
    defineConstant("simnum", getSeed());
    setSeed(getSeed() + $RANDOM);
    sim.addSubpop("p1",$( echo "${Nanc} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
}

$( echo "( ${Nanc} * 10 ) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) early() { // after burn-in, split populations into two: p1 and p2
    sim.addSubpopSplit("p2", $( echo "${Nr1} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ), p1);
    p1.setSubpopulationSize($( echo "${Ns} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
}

$( echo "( ${Nanc} * 10 + ${splittime} ) / ${scalingfactor} - 1" | bc -l | xargs printf "%.0f" ) late() {
    sim.treeSeqRememberIndividuals(sim.subpopulations.individuals);
}

$( echo "( ${Nanc} * 10 + ${splittime} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) early() {
    sim.addSubpop("p3", $( echo "${Nr2} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
    p3.setMigrationRates(c(p1,p2), c(0.05, 0.95));
}

$( echo "( ${Nanc} * 10 + ${splittime} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) late() {
    p3.setMigrationRates(c(p1,p2), c(0, 0));
    
    p1.setSubpopulationSize(0);
    p2.setSubpopulationSize(0);
}

$( echo "( ${Nanc} * 10 + ${splittime} + ${Nanc} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f") {
    p3.setSubpopulationSize(1000);
}

$( echo "( ${Nanc} * 10 + ${splittime} + ${Nanc} ) / ${scalingfactor}" | bc -l | xargs printf "%.0f") late() {
    sim.treeSeqOutput("/u/home/b/bkim331/project-klohmueldata/bernard_data/admixture_simulations_fullsims/burnin_${model}/h_${h}/h${h}_full_" + simnum + "_p1.trees");
}

EOM

cd ..
done


cd ..