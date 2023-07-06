
process CCDSFilter {
  
   input:
   tuple val(chrom), val(pos), val(id), val(ref), val(alt), val(qual), val(filter), val(info) 

   output:
   path('CCDS.vcf')
   
   shell:
   '''
   loc=`pwd`
   echo -e "!{chrom}\t!{pos}\t!{id}\t!{ref}\t!{alt}\t!{qual}\t!{filter}\t!{info}" > one.vcf
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/one.vcf" --output_file "${loc}/CCDS.vcf" -plugin CCDSFilter
   '''
}

process CSN {

   input:
   path(in_vcf)

   output:
   path('CSN.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/CSN.vcf" -plugin CSN
   '''
}

process Carol {

   input:
   path(in_vcf)

   output:
   path('Carol.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/Carol.vcf" -plugin Carol
   '''
}

process Downstream {

   input:
   path(in_vcf)

   output:
   path('Downstream.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/Downstream.vcf" -plugin Downstream
   '''
}
 
process Draw {
   
   input:
   path(in_vcf)

   output:
   path('Draw.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/Draw.vcf" -plugin Draw
   '''
}

process GXA {

   input:
   path(in_vcf)

   output:
   path('GXA.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/GXA.vcf" -plugin GXA
   '''
}

process LOVD { 

   input:
   path(in_vcf)

   output:
   path('LOVD.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/LOVD.vcf" -plugin LOVD
   '''
}

process NMD {

   input:
   path(in_vcf)

   output:
   path('NMD.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/NMD.vcf" -plugin NMD
   '''
}

process NearestGene {

   input:
   path(in_vcf)

   output:
   path('NearestGene.vcf')

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/NearestGene.vcf" -plugin NearestGene
   '''
}

process AncestralAllele {
   
   input:
   path(in_vcf)

   output:
   path('AncestralAllele.vcf'), emit: vcf

   shell:
   '''
   loc=`pwd`
   vep --cache --offline --format vcf --vcf --force_overwrite --input_file "${loc}/!{in_vcf}" --output_file "${loc}/AncestralAllele.vcf" -plugin AncestralAllele,homo_sapiens_ancestor_GRCh38.fa.gz 
   '''
}

process makeHeader {
	
   input:
   each item

   shell:
   '''
   echo !{item} >> header.txt
   '''
}

process CompileVCF {

   publishDir "${params.outdir}", mode: "copy"

   input:
   path '*.vcf'
   path 'input.vcf.gz'

   output:
   path 'final.vcf.gz'

   //recreate VCF using original header plus any extra info fields created by plugins
   //some plugins only update the "latest command" header line, I omitted that due to no additional info provided and is complex
   shell:
   '''
   zcat input.vcf.gz | grep -e "^##" > final.vcf
   grep -e "^##INFO=<ID=CSQ" *.vcf | uniq | cut -f 2 -d : | uniq >> final.vcf
   grep -e "^##CAROL" *.vcf | uniq | cut -f 2 -d : | uniq >> final.vcf
   grep -e "^##LOVD" *.vcf | uniq | cut -f 2 -d : | uniq >> final.vcf
   grep -e "^##NMD" *.vcf | uniq | cut -f 2 -d : | uniq >> final.vcf
   grep -e "^##NearestGene" *.vcf | uniq | cut -f 2 -d : | uniq >> final.vcf
   zcat input.vcf.gz | grep -e "^#CHROM" >> final.vcf
   tail -n 1 *.vcf | grep -v "==>" | grep -e "^." >> final.vcf
   gzip final.vcf
   '''
}

workflow {
   input_vcf = Channel.fromPath(params.input)
   //splits input vcf into a string of each line, split converts this to a map of each line then to toList converts it to the expected list type for future operators
   input_vcf.splitText( decompress:true ) { it.split().toList() } 
   //using the logic of how vcf files are built, put all header information into one channel and all genomic regions into another
	    .branch {
		header: it.first() =~ '#' 
		body: it.first() !=~ '#'
	     }
	    .set { processed_vcf }
   CCDSFilter(processed_vcf.body) | CSN | Carol | Downstream | Draw | GXA | LOVD | NMD | NearestGene | AncestralAllele
   CompileVCF(AncestralAllele.out.vcf.collect(), input_vcf)
}
