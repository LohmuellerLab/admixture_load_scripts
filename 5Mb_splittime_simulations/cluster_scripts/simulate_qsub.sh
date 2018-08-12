#! /bin/bash

scalingfactor="5"

for splittime in "100" "250" "500" "1000" "2500" "5000" "10000" "15000" "20000" "25000" "30000" "35000" "40000"; do

#demographic parameters before scaling
#time in generations from 0. remember that 10N generations are used for burn-in
# Nanc: size of ancestral population before divergence
# Nsource: size of p1 from time 10N to 13N
# Nrecipient1: size of p2 from time 10N to (12N-50)
# NrecipientB: size of p2 from time (12N-50) to (12N)
# Nrecipient2: size of p2 from time 12N to 13N
model="model4"
Nanc="10000"
Nsource="10000"
Nrecipient1="1000"
NrecipientB="1000"
Nrecipient2="10000"

cwd=$( pwd )

mkdir sims
cd sims

mkdir ${model}
cd ${model}

for h in "0.0"; do
mkdir split_${splittime}; cd split_${splittime}

if [ ${h} = "0.5" ]; then
read -r -d '' fitnessblock << endmsg 
1:$( echo "130000 / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) fitness(m1) {
    sadj = mut.selectionCoeff * mut.mutationType.dominanceCoeff;
    if (homozygous) {
        return ((1.0 + sadj)*(1.0 + sadj));
    } else {
        return (1.0 + sadj);
    }
}
endmsg
else 
    fitnessblock=""
fi

for r in "1e-9"; do

cat > slim_h${h}_r${r}_${scalingfactor}x.slim << EOM
function (f$)calcFST(o<Subpopulation>$ subpop1, o<Subpopulation>$ subpop2)
{
    mutlist = c(sim.mutationsOfType(m1), sim.mutationsOfType(m2), sim.mutationsOfType(m3));
    p1_p = sim.mutationFrequencies(subpop1, mutlist);
    p2_p = sim.mutationFrequencies(subpop2, mutlist);
    n1 = size(subpop1.individuals);
    n2 = size(subpop2.individuals);
    fst_num = (p1_p-p2_p)^2 - ((p1_p*(1-p1_p))/(n1-1)) - ((p2_p*(1-p2_p))/(n2-1));
    fst_denom = p1_p*(1-p2_p) + p2_p*(1-p1_p);
    fst_num = mean(fst_num);
    fst_denom = mean(fst_denom);
    fst = fst_num/fst_denom;
    return fst;
}

function (f$)calcFST2(o<Subpopulation>$ subpop1, o<Subpopulation>$ subpop2)
{
    mutlist = c(sim.mutationsOfType(m1), sim.mutationsOfType(m2), sim.mutationsOfType(m3));
    p1_p = sim.mutationFrequencies(subpop1, mutlist);
    p2_p = sim.mutationFrequencies(subpop2, mutlist);
    mean_p = (p1_p + p2_p) / 2.0;
    H_t = 2.0 * mean_p * (1.0 - mean_p);
    H_s = p1_p * (1.0 - p1_p) + p2_p * (1.0 - p2_p);
    fst_num = H_t - H_s;
    fst_denom = H_t;
    return mean(fst_num)/mean(fst_denom);
}
  
function (f$)calcFST3(o<Subpopulation>$ subpop1, o<Subpopulation>$ subpop2){
    mutlist = c(sim.mutationsOfType(m1), sim.mutationsOfType(m2), sim.mutationsOfType(m3));
    p1_p = sim.mutationFrequencies(subpop1, mutlist);
    p2_p = sim.mutationFrequencies(subpop2, mutlist);
    mean_p = (p1_p + p2_p) / 2.0;
    H_t = 2.0 * mean_p * (1.0 - mean_p);
    H_s = p1_p * (1.0 - p1_p) + p2_p * (1.0 - p2_p);
    fst = 1.0 - H_s/H_t;
    return mean(fst);
}

// use page 76 to randomly generate exons
initialize() {
    initializeMutationRate(1.5e-8*${scalingfactor});
	
	//nonsynonymous drawn from a DFE from Kim et al.
	// importantly, SLiM computes the fitness of the heterozygote and homozygote as 1+sh and 1+s
	// dadi and others compute it as 1+2sh and 1+2s
        initializeMutationType("m1", ${h}, "g", -0.01314833*${scalingfactor}, 0.186);
	// initializeMutationType("m1", ${h}, "f", 0.0);
	//synonymous -- assumed neutral here
	initializeMutationType("m2", ${h}, "f", 0.0);
	//noncoding -- assumed neutral here
	initializeMutationType("m3", ${h}, "f", 0.0);
	//beneficial -- for the time being, is left out
	//initializeMutationType("m4", ${h}, "f", 0.0);
	
	// p1 marker mutations
	// m10 is within genes, m11 is outside of genes
	initializeMutationType("m10", 0.5, "f", 0.0);
	m10.convertToSubstitution = F; //marker mutations are always tracked
	initializeMutationType("m11", 0.5, "f", 0.0);
	m11.convertToSubstitution = F; //marker mutations are always tracked
	
	//genomic element: exon and uses a mixture of syn and nonsyn at a 1:2.31 ratio (Huber et al.)
	initializeGenomicElementType("g1", c(m2,m1), c(1.0,2.31));
	//genomic element: intron
	initializeGenomicElementType("g2", c(m3), c(1.0));
	//genomic element: intergenic
	initializeGenomicElementType("g3", c(m3), c(1.0));	
	
	// Generate random genes along approximately 100kb
	base = 0;
	
	while (base < 5000000) {
		//make first noncoding
		nc_length = asInteger(runif(1, 100, 5000));
		initializeGenomicElement(g3, base, base + nc_length - 1);
		base = base + nc_length;
		
		//make first exon
		ex_length = asInteger(rlnorm(1, log(50), log(2))) + 1;
		initializeGenomicElement(g1, base, base + ex_length - 1);
		base = base + ex_length;
		
		//make additional intron-exon pairs
		do
		{
			in_length = asInteger(rlnorm(1, log(100), log(1.5))) + 10;
			initializeGenomicElement(g2, base, base + in_length - 1);
			base = base + in_length;
			
			ex_length = asInteger(rlnorm(1, log(50), log(2))) + 1;
			initializeGenomicElement(g1, base, base + ex_length - 1);
			base = base + ex_length;
		}
		while (runif(1) < 0.8); //20% probability of stopping
	}
	
	nc_length = asInteger(runif(1, 100, 5000));
	initializeGenomicElement(g3, base, base + nc_length - 1);
	
	// one recombination rate
	initializeRecombinationRate(${r}*${scalingfactor});
}

${fitnessblock}

//burn in
1 early() {
    defineConstant("simnum", getSeed());
    setSeed(getSeed() + $RANDOM);
    sim.addSubpop("p1", $( echo "${Nanc} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
    // when set to "F" does not remove fixed variants
    // this isn't necessary because we're keeping track of all populations
    m2.convertToSubstitution = T;
}

$( echo "100000 / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) early() { // after burn-in, split populations into two: p1 and p2
    sim.addSubpopSplit("p2", $( echo "${Nrecipient1} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ), p1);
    p1.setSubpopulationSize($( echo "${Nsource} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
    
    // this also isnt necessary, but sets the migration rates to 0
	p1.setMigrationRates(c(p2),c(0.0)); //migration rate INTO p1
	p2.setMigrationRates(c(p1),c(0.0)); //migration rate INTO p2
	
    //get start and stop vectors for noncoding
    starts = c();
    ends = c();
    for (ge in sim.chromosome.genomicElements)
        if (ge.genomicElementType == g3) {
            starts=c(starts,ge.startPosition);
            ends=c(ends,ge.endPosition);
            }
  	//add marker mutations every 500 bp
  	withingenes = 0;
  	outsidegenes = 0;
	for (pos in asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/500))*500)-1)) {
	    if (pos <= ends[sum(pos >= starts) - 1]) { // if outside genes
		    p1.genomes.addNewDrawnMutation(m11, pos); // mark with m11
		    outsidegenes = outsidegenes + 1;
		}
		else { // if in genes
		    p1.genomes.addNewDrawnMutation(m10, pos); // mark with m10
		    withingenes = withingenes + 1;
		}
	}
	
	defineConstant("within", withingenes);
	defineConstant("outside", outsidegenes);
	
    cat("#generation,p1Fraction,p2Fraction,FST_Bhatia,FST_Hudson2,FST_mean\n");}

$( echo "(100000 + ${splittime}) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) early() {
    p1.setMigrationRates(c(p2),c(0.00)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.05)); //migration rate INTO p2
}

$( echo "(100000 + ${splittime}) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) late() {
    p1.setMigrationRates(c(p2),c(0.0)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.0)); //migration rate INTO p2
    p2.setSubpopulationSize($( echo "${Nrecipient2} / ${scalingfactor}" | bc -l | xargs printf "%.0f" ));
}

$( echo "(110000 + ${splittime}) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) { }

$( echo "100000 / ${scalingfactor}" | bc -l | xargs printf "%.0f" ):$( echo "(110000 + ${splittime}) / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) late() {
  if (sim.generation % $( echo "50 / ${scalingfactor}" | bc -l | xargs printf "%.0f" ) == 0){
    
    p1g = p1.genomes;
    p2g = p2.genomes;
    
    //quantify introgressed ancestry
    p1Count = sum(p1g.countOfMutationsOfType(m10)) + sum(p1g.countOfMutationsOfType(m11));
    maxp1Count = p1g.size() * size(asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/500))*500)-1));
    p1Fraction = 1 - (p1Count / maxp1Count); // one minus here because im looking at the fraction of p2 ancestry in p1
        
    //quantify introgressed ancestry
    // p2g = p2.genomes; //defined earlier up
    p2Count = sum(p2g.countOfMutationsOfType(m10)) + sum(p2g.countOfMutationsOfType(m11));
    maxp2Count = p2g.size() * size(asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/500))*500)-1));
    p2Fraction = p2Count / maxp2Count;
    
    fst1 = calcFST(p1,p2);
    fst2 = calcFST2(p1,p2);
    fst3 = calcFST3(p1,p2);
    
    cat(sim.generation + "," + p1Fraction + "," + p2Fraction + "," + fst1 + "," + fst2 + "," + fst3 + "\n");
    }
}

EOM

queue="h_rt=24:00:00,h_data=8G,arch=intel*"

cat > submitjob.sh << EOM
#! /bin/bash
#$ -cwd
#$ -A bkim331
#$ -l ${queue}
#$ -N slim_${model}_s${splittime}
#$ -o /dev/null
#$ -j y
#$ -M bkim331@gmail.com
#$ -m a
#$ -t 1-200:1
. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
cd ${cwd}/sims/${model}/split_${splittime}
/u/home/b/bkim331/bin/slim -seed \${SGE_TASK_ID} slim_h${h}_r${r}_${scalingfactor}x.slim > slim_h${h}_r${r}_\${SGE_TASK_ID}.csv
EOM

#qsub submitjob.sh

done

cd ..

done

cd ..
cd ..

done