#! /bin/bash

#read -r -d '' fitnessblock << endmsg 
#1:13000 fitness(m1) {
#    h = 1/((1/0.987) * 39547 * mut.selectionCoeff)
#    if (homozygous) {
#        return ((1.0 + 0.5*mut.selectionCoeff)*(1.0 + 0.5*mut.selectionCoeff));
#    } else {
#        return (1.0 + mut.selectionCoeff * h);
#    }
#}
#endmsg

#scaling factor should be 1, 2, 5, or 10
scalingfactor="5"

Nanc="10000"
Nsource="10000"
Nrecipient1="10000"
NrecipientB="1000"
Nrecipient2="10000"

for h in 0.5; do
mkdir h_${h}_model1; cd h_${h}_model1

read -r -d '' fitnessblock << endmsg 
1:$( expr 130000 / ${scalingfactor} ) fitness(m1) {
    sadj = mut.selectionCoeff * mut.mutationType.dominanceCoeff;
    if (homozygous) {
        return ((1.0 + sadj)*(1.0 + sadj));
    } else {
        return (1.0 + sadj);
    }
}
endmsg

for r in "1e-6" "1e-7" "1e-8" "1e-9"
do

cat > slim_h${h}_r${r}_model1_${scalingfactor}x.slim << EOM
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
    sim.addSubpop("p1",$( expr ${Nanc} / ${scalingfactor} ));
    // when set to "F" does not remove fixed variants
    // this isn't necessary because we're keeping track of all populations
    m2.convertToSubstitution = T;
}

$( expr 100000 / ${scalingfactor} ) early() { // after burn-in, split populations into two: p1 and p2
    sim.addSubpopSplit("p2", $( expr ${Nrecipient1} / ${scalingfactor} ), p1);
    p1.setSubpopulationSize($(expr ${Nsource} / ${scalingfactor} ));
    
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
	
    cat("#generation,FitnessLarge,FitnessSmall,p1Fraction,p1FractionWithin,p1FractionOutside,p2Fraction,p2FractionWithin,p2FractionOutside,p1MeanDel,p1MeanNeu,p1MeanMut,p1IndivDel,p1MeanHomNeu,p1MeanHomDel,p1load_s,p1load_w,p2MeanDel,p2MeanNeu,p2MeanMut,p2IndivDel,p2MeanHomNeu,p2MeanHomDel,p2load_s,p2load_w\n");}

$( expr 119950 / ${scalingfactor} ) early() {
    p2.setSubpopulationSize($(expr ${NrecipientB} / ${scalingfactor} ));
}

$( expr 120000 / ${scalingfactor} ) early() {
    p1.setMigrationRates(c(p2),c(0.00)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.05)); //migration rate INTO p2
}

$( expr 120000 / ${scalingfactor} ) late() {
    p1.setMigrationRates(c(p2),c(0.0)); //migration rate INTO p1
    p2.setMigrationRates(c(p1),c(0.0)); //migration rate INTO p2
    p2.setSubpopulationSize($(expr ${Nrecipient2} / ${scalingfactor} ));
}

$( expr 130000 / ${scalingfactor} ) { }

$( expr 100000 / ${scalingfactor} ):$( expr 130000 / ${scalingfactor} ) late() {
  if (sim.generation % $(expr 50 / ${scalingfactor}) == 0){
    
    p1g = p1.genomes;
    p2g = p2.genomes;
    
    //quantify introgressed ancestry
    p1Count = sum(p1g.countOfMutationsOfType(m10)) + sum(p1g.countOfMutationsOfType(m11));
    maxp1Count = p1g.size() * size(asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/500))*500)-1));
    p1Fraction = 1 - (p1Count / maxp1Count); // one minus here because im looking at the fraction of p2 ancestry in p1
    
    //within genes
    p1CountWithin = sum(p1g.countOfMutationsOfType(m10));
    maxp1CountWithin = p1g.size() * within;
    p1FractionWithin = 1 - (p1CountWithin / maxp1CountWithin); // one minus here because im looking at the fraction of p2 ancestry in p1
    
    //outside of genes
    p1CountOutside = sum(p1g.countOfMutationsOfType(m11));
    maxp1CountOutside = p1g.size() * outside;
    p1FractionOutside = 1 - (p1CountOutside / maxp1CountOutside); // one minus here because im looking at the fraction of p2 ancestry in p1
    
    //quantify average number of deleterious mutations per haplotype
    p1MeanDel = mean(p1g.countOfMutationsOfType(m1));
    p1MeanNeu = mean(p1g.countOfMutationsOfType(m2)+p1g.countOfMutationsOfType(m3));
    p1MeanMut = mean(p1g.countOfMutationsOfType(m1)+p1g.countOfMutationsOfType(m2)+p1g.countOfMutationsOfType(m3));
    
    //quantify average number of deleterious mutations per individual
    p1IndivDel = mean(p1.individuals.countOfMutationsOfType(m1));
    p1HomsDel = sum(p1.individuals.countOfMutationsOfType(m1)) - size(p1.individuals.uniqueMutationsOfType(m1));
    p1HomsNeu = sum(p1.individuals.countOfMutationsOfType(m2)) - size(p1.individuals.uniqueMutationsOfType(m2)) + sum(p1.individuals.countOfMutationsOfType(m3)) - size(p1.individuals.uniqueMutationsOfType(m3));
    p1MeanHomDel = p1HomsDel/size(p1.individuals);
    p1MeanHomNeu = p1HomsNeu/size(p1.individuals);
    
    //quantify reduction in fitness due to...
    fitnessvec_s = c();
    fitnessvec_w = c();
    for (indiv in p1.individuals){
        allmuts = c(indiv.genomes[0].mutationsOfType(m1),indiv.genomes[1].mutationsOfType(m1));
        uniqmuts = indiv.uniqueMutationsOfType(m1);
        
        if (size(uniqmuts) > 0){
            indivfitnessvec_s = c();
            indivfitnessvec_w = c();
        
            for (u in uniqmuts){
                places = (allmuts.id == u.id);
                uu = allmuts[places];
                reductionFitness = 0;
                if ((m1.dominanceCoeff == 0.0) & (size(uu) == 2)) {
                    reductionFitness = sum(uu.selectionCoeff)/2;
                } else if ((m1.dominanceCoeff == 0.0) & (size(uu) == 1)) {
                    reductionFitness = 0;
                } else if (m1.dominanceCoeff == 0.5) {
                    reductionFitness = product(1 + uu.selectionCoeff * m1.dominanceCoeff) - 1; //add/mult
                    // reductionFitness = sum(prod(1-uu.selectionCoeff))/2; //strict additive
                }
                
                if (2*max(c(p1g.size(),p2g.size()))*abs(u[0].selectionCoeff) <= 1){
                    indivfitnessvec_s = c(indivfitnessvec_s, 1);
                    indivfitnessvec_w = c(indivfitnessvec_w, 1+reductionFitness);
                } else {
                    indivfitnessvec_s = c(indivfitnessvec_s, 1+reductionFitness);
                    indivfitnessvec_w = c(indivfitnessvec_w, 1);
                }
            }
            
            fitnessvec_s = c(fitnessvec_s, product(indivfitnessvec_s));
            fitnessvec_w = c(fitnessvec_w, product(indivfitnessvec_w));
            
        } else {
            fitnessvec_s = c(fitnessvec_s, 1);
            fitnessvec_w = c(fitnessvec_w, 1);
        }
    }
    
    meanFitnessP1 = mean(fitnessvec_s*fitnessvec_w);
    p1load_s=mean(fitnessvec_s);
    p1load_w=mean(fitnessvec_w);
    
    //quantify introgressed ancestry
    // p2g = p2.genomes; //defined earlier up
    p2Count = sum(p2g.countOfMutationsOfType(m10)) + sum(p2g.countOfMutationsOfType(m11));
    maxp2Count = p2g.size() * size(asInteger(c(1,(1:((sim.chromosome.lastPosition+1)/500))*500)-1));
    p2Fraction = p2Count / maxp2Count;
    
    //within genes
    p2CountWithin = sum(p2g.countOfMutationsOfType(m10));
    maxp2CountWithin = p2g.size() * within;
    p2FractionWithin = (p2CountWithin / maxp2CountWithin);
    
    //outside of genes
    p2CountOutside = sum(p2g.countOfMutationsOfType(m11));
    maxp2CountOutside = p2g.size() * outside;
    p2FractionOutside = (p2CountOutside / maxp2CountOutside);
    
    //quantify average number of deleterious mutations per haplotype
    p2MeanDel = mean(p2g.countOfMutationsOfType(m1));
    p2MeanNeu = mean(p2g.countOfMutationsOfType(m2)+p2g.countOfMutationsOfType(m3));
    p2MeanMut = mean(p2g.countOfMutationsOfType(m1)+
    p2g.countOfMutationsOfType(m2)+p2g.countOfMutationsOfType(m3));
    
    //quantify average number of deleterious mutations per individual
    p2IndivDel = mean(p2.individuals.countOfMutationsOfType(m1));
    p2HomsDel = sum(p2.individuals.countOfMutationsOfType(m1)) - size(p2.individuals.uniqueMutationsOfType(m1));
    p2HomsNeu = sum(p2.individuals.countOfMutationsOfType(m2)) - size(p2.individuals.uniqueMutationsOfType(m2)) + sum(p2.individuals.countOfMutationsOfType(m3)) - size(p2.individuals.uniqueMutationsOfType(m3));
    p2MeanHomDel = p2HomsDel/size(p2.individuals);
    p2MeanHomNeu = p2HomsNeu/size(p2.individuals);
    
    //quantify reduction in fitness due to...
    fitnessvec_s = c();
    fitnessvec_w = c();
    for (indiv in p2.individuals){
        allmuts = c(indiv.genomes[0].mutationsOfType(m1),indiv.genomes[1].mutationsOfType(m1));
        uniqmuts = indiv.uniqueMutationsOfType(m1);
        
        if (size(uniqmuts) > 0){
            indivfitnessvec_s = c();
            indivfitnessvec_w = c();
        
            for (u in uniqmuts){
                places = (allmuts.id == u.id);
                uu = allmuts[places];
                reductionFitness = 0;
                if ((m1.dominanceCoeff == 0.0) & (size(uu) == 2)) {
                    reductionFitness = sum(uu.selectionCoeff)/2;
                } else if ((m1.dominanceCoeff == 0.0) & (size(uu) == 1)) {
                    reductionFitness = 0;
                } else if (m1.dominanceCoeff == 0.5) {
                    reductionFitness = product(1 + uu.selectionCoeff * m1.dominanceCoeff) - 1;
                    // reductionFitness = sum(uu.selectionCoeff)/2;
                }
                
                if (2*max(c(p1g.size(),p2g.size()))*abs(u[0].selectionCoeff) <= 1){
                    indivfitnessvec_s = c(indivfitnessvec_s, 1);
                    indivfitnessvec_w = c(indivfitnessvec_w, 1+reductionFitness);
                } else {
                    indivfitnessvec_s = c(indivfitnessvec_s, 1+reductionFitness);
                    indivfitnessvec_w = c(indivfitnessvec_w, 1);
                }
            }
            
            fitnessvec_s = c(fitnessvec_s, product(indivfitnessvec_s));
            fitnessvec_w = c(fitnessvec_w, product(indivfitnessvec_w));
        
        } else {
            fitnessvec_s = c(fitnessvec_s, 1);
            fitnessvec_w = c(fitnessvec_w, 1);
        }
    }
    
    meanFitnessP2 = mean(fitnessvec_s*fitnessvec_w);
    p2load_s=mean(fitnessvec_s);
    p2load_w=mean(fitnessvec_w);
    
    cat(sim.generation + "," + meanFitnessP1 + "," + meanFitnessP2 + "," + p1Fraction + "," + p1FractionWithin + "," + p1FractionOutside + "," + p2Fraction + ","+ p2FractionWithin + "," + p2FractionOutside + ","  +  p1MeanDel + "," + p1MeanNeu + "," + p1MeanMut + "," + p1IndivDel + "," + p1MeanHomNeu + "," + p1MeanHomDel + "," + p1load_s + "," + p1load_w + "," + p2MeanDel + "," + p2MeanNeu + "," + p2MeanMut + "," + p2IndivDel + "," + p2MeanHomNeu + "," + p2MeanHomDel + "," + p2load_s + "," + p2load_w + "\n");
    }
}

EOM

done

cd ..
done