#!/usr/bin/perl
use Cwd;
use File::Spec::Functions qw(rel2abs); 
use File::Basename; 
use Scalar::Util qw(looks_like_number);

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# PEGASUS Version 1.1
# Priyanka Nakka, Sohini Ramachandran Lab, Brown University 2015
# modified from VEGAS (Liu et al AJHG 2010) 

$dir = dirname(rel2abs($0));

print "Options in effect: \n";
print "-input $ARGV[0]\n";

# Options
for($i = 1; $i <= $#ARGV; ++$i){
	$field = $ARGV[$i];
	if(length($field)>1){
		if($field eq "-pop"){ # Specify reference pop
			$chip = "$dir/$ARGV[$i+1]";
			print "-pop $ARGV[$i+1]\n";
		}
		if($field eq "-chr"){ # Do test for specific chromosome
			$dochr = "$ARGV[$i+1]";
			print "-chr $dochr\n";
			if ($dochr > 23 || $dochr < 1){
				die "Error: -chr must be a number between 1 and 23\n";
			}
		}
		if($field eq "-out"){ # Specify outfile
			print "-out $ARGV[$i+1]\n";
			$outfile = "$ARGV[$i+1]";
		}
		if($field eq "-keeptimestamp"){ # Do not delete time_stamp folder
			print "-keeptimestamp\n";
			$keeptimestamp = 1;
		}
		if($field eq "-custom"){ # Use custom individual genotypes 
			print "-custom $ARGV[$i+1]\n";
			$custom = "$ARGV[$i+1]";
			$chip = "$dir/hapmapCEU";
			#$chip = "$dir/1kgEUR";
		}
		if($field eq "-upper"){ # Upper boundary (use with -custom or -ld-file only)
			print "-upper $ARGV[$i+1]\n";
			$upper = "$ARGV[$i+1]";
		}
		if($field eq "-lower"){ # Lower boundary (use with -custom or -ld-file only)
			print "-lower $ARGV[$i+1]\n";
			$lower = "$ARGV[$i+1]";
		}
		if($field eq "-glist"){ # Lower boundary (use with -custom or -ld-file only)
			print "-glist $ARGV[$i+1]\n";
			$glist = "$ARGV[$i+1]";
		}
		if($field eq "-ld-file"){ # Specify user-generated LD file in default PLINK format
			print "-ld-file $ARGV[$i+1]\n";
			$chip = "$dir/hapmapCEU";
			#$chip = "$dir/1kgEUR";
			$ldfile = "$ARGV[$i+1]";	
		}
		if($field eq "-fisher"){ # calculate gene scores for males and females separately and combine using fisher's method
			$fisher = 1;	
		}
	}
}

# Make working directory
my ($Second, $Minute, $Hour, $Day, $Month, $Year, $WeekDay, $DayOfYear, $IsDST) = localtime(time) ;
$time_stamp = $DayOfYear.$Hour.$Minute.$Second.$dochr;
mkdir "$time_stamp", 0777 or die "Error: cannot make working directory. Check permission settings: $!";

print "\nCreating working directory: $time_stamp\n";

# Check that p-values file exists
if(!-e "$ARGV[0]"){
	die "Error: $ARGV[0] not found\n";
}

open(PVALS, "$ARGV[0]"); # Location of GWAS results file
@pvals = <PVALS>;

# Check that custom genotypes exist - and that length of $custom.bim is the same as @pvals
if(defined($custom)){
	if(!-e "$custom.bed"){
		die "Error: $custom.bed not found\n";
	}
	if(!-e "$custom.bim"){
		die "Error: $custom.bim not found\n";
	}
	if(!-e "$custom.fam"){
		die "Error: $custom.fam not found\n";
	}
}

# Check upper and lower are being used correctly
if(defined($upper) || defined($lower)){
	if(!defined($custom)){
		if(!defined($ld-file)){
			die "Error: -upper and -lower can only be used for custom individual genotypes or ld files\n";
		}
	}
}

chdir "$time_stamp" or die "Error: cannot change to working directory: $!";

# Check that R and plink exists
system("echo -n \$PATH > path");

open(PATH, "path");
$path = <PATH>;
@path = split(/:/, $path);

$Rtimes=0;
$plinktimes=0;

foreach $i(@path){
	if(-e "$i/R"){
		$Rtimes=$Rtimes+1;
	}
	if(-e "$i/plink"){
		$plinktimes=$plinktimes+1;
	}
}
if($Rtimes==0){
	chdir "../" or die "cannot change directory: $!";
	system("rm -R $time_stamp");
	die "Error: Cannot find R - make sure it can be accessed through \$PATH \n$!\n";
}
if($plinktimes==0){
	chdir "../" or die "cannot change directory: $!";
	system("rm -R $time_stamp");
	die "Error: Cannot find plink - make sure it can be accessed through \$PATH \n$!\n";
}

close PATH;
system("rm path");

# Check that corpcor and compquadform exist
&checkpackages;

# Write files to server
open(TEMPPVALS, ">temppvals");
print TEMPPVALS "@pvals";

close PVALS;
close TEMPPVALS;

# Validations

# Check pvalues in file are valid - in decimal or standard scientific notation (0.0034, 1e-6 etc...)
print "\nValidating input file...";
&numbercheck;
print "done.\n";

# More validations
if(!defined($custom)){
	if (!defined($ldfile)) {
		if (!defined($chip)){
			die "Error: Invalid SNP-set selection\n";
		}
		if (!-d "$chip/geneset"){
			die "Error: $chip/geneset does not exist\n";
		}
		if (!-d "$chip/genebeds"){
			die "Error: $chip/genebeds does not exist\n";
		}
	}
}

# Do test
system("echo \"Chr Gene nSNPs Start Stop Test-Stat Pvalue Error Best-SNP SNP-pvalue\" > gene-basedoutput.out");

	
if(defined($custom)){
	if(defined($dochr)){
		&docustomchr;
	}
	else{
		&customtest;
	}
}
elsif(defined($ldfile)){
	if(defined($dochr)){
		&doldchr;
	}
	else{
		&ldfiletest
	}
}
else{
	if(defined($dochr)){
		&dochrtest;
	}
	else{
		&defaulttest;
	}
}
# Clean up files

if(defined($outfile)){
	system("cp gene-basedoutput.out ../$outfile.out");
}
else{
	system("cp gene-basedoutput.out ../gene-basedoutput.out");
}

chdir "../" or die "Error: Cannot change directory: $!";
unless(defined($keeptimestamp)){
	system("rm -R $time_stamp");
}

print "\nGene-based test complete!\n";


#Subroutines:
sub defaulttest{ 

	print "\nAnnotating GWAS results...";
	&makemerger;
	print "done\n\n";

	for ($chr = 1; $chr <= 23; $chr++) {
		if ($chr==23) {$chr="X";}
		open(MERGEFILE, "merged$chr"); # Check if there are genes on chromosome
		my @mergefile = <MERGEFILE>;
		close(MERGEFILE);
		if(scalar(@mergefile) == 1){
			next;
		}
		
		print "Starting chromosome $chr\n";
		open(ALLGENES, "$chip/geneset/allgene$chr");

		@allgenes = <ALLGENES>;
		close ALLGENES;
	
		&dogene;
	}
}

sub dochrtest{ # Chromosome test
	print "\nAnnotating GWAS results...";
		&makemergerdochr;
	print "done\n\n";

	for ($chr = $dochr) {

		print "Starting chromosome $chr\n";
	
		open(ALLGENES, "$chip/geneset/allgene$chr");

		@allgenes = <ALLGENES>;
		close ALLGENES;
	
		&dogene;
	}
}

sub makemerger{ # Makemerge - genome-wide
system("
R --vanilla --slave  <<EOF

options(warn=-1)

#read in GWAS results
pvals <- read.table('temppvals',header=F,colClasses=c('character','numeric'))
names(pvals) <- c('SNP','PVALUE')

#read in start/stop positions for each gene
#read.table('$chip/glist-hg18',colClasses=c('character',rep('integer',2),'character')) -> startstop
read.table('$chip/glist-hg19',colClasses=c('character',rep('integer',2),'character')) -> startstop
names(startstop) <- c('chr','start','stop','gene')

snpgenevec <- seq(1,23)

for(chr in snpgenevec){
	if (chr == 23) {chr = 'X'}
	filename <- paste('$chip/geneset/snpgene',chr,sep='')
	#read in which SNPs map to which gene
	read.table(file=filename,colClasses=c('character','character')) -> wp
	names(wp) <- c('gene','rs')

	merge(startstop,wp,by='gene',all.y=T) -> wp

	#convert p-values to chi2-df1 stats
	qchisq(pvals[,2],1,lower.tail=F) -> chi2 
	data.frame(pvals[,1],chi2) -> chi2withrs
	names(chi2withrs) <- c('SNP','chisq')

	merge(wp,chi2withrs,by.x='rs',by.y='SNP') -> mer

	writename <- paste('merged',chr,sep='')
	write.table(mer,file=writename,row.names=F,quote=F)
}

EOF");
}

sub makemergerdochr{ # Makemerge - single chromosome

system("
R --vanilla --slave  <<EOF

options(warn=-1)

#read in GWAS results
pvals <- read.table('temppvals',header=F,colClasses=c('character','numeric'))
names(pvals) <- c('SNP','PVALUE')

#read in start/stop positions for each gene
#read.table('$chip/glist-hg18',colClasses=c('character',rep('integer',2),'character')) -> startstop
read.table('$chip/glist-hg19',colClasses=c('character',rep('integer',2),'character')) -> startstop
names(startstop) <- c('chr','start','stop','gene')

#read in which SNPs map to which gene
read.table('$chip/geneset/snpgene$dochr',colClasses=c('character','character')) -> wp
names(wp) <- c('gene','rs')

merge(startstop,wp,by='gene',all.y=T) -> wp

#convert p-values to chi2-df1 stats
qchisq(pvals[,2],1,lower.tail=F) -> chi2 
data.frame(pvals[,1],chi2) -> chi2withrs
names(chi2withrs) <- c('SNP','chisq')

merge(wp,chi2withrs,by.x='rs',by.y='SNP') -> mer

write.table(mer,file='merged$dochr',row.names=F,quote=F)

EOF");
}

sub dogene{ # Main gene-based test outline
			
	foreach $gene (@allgenes){
		chomp($gene);
		
		open(MERGEDFILE,"merged$chr");
		@tempgene = grep /\b $gene \b/i, <MERGEDFILE>;
		close(MERGEDFILE);
		
		open(PRINTTEMPGENE,">tempgene");
		print PRINTTEMPGENE "@tempgene";
		close(PRINTTEMPGENE);
		
		if(-z "tempgene"){next;} # If file is empty, move on
		
		system("awk  '{print \$1;}' tempgene > tempgene.snp");
		
		open(PREDEFGENE, "$chip/genebeds/$gene.bim");
		@predefgene = <PREDEFGENE>;
		close PREDEFGENE;

		#check that gene exists in geneset
		#if (scalar(@predefgene) == 0){
		#	print "Gene doesn't exist...\n";
		#}
		#else{	
			&hapsims;
		#}
	}
	
}

sub hapsims {	
	&plink;
				system("

R --vanilla --slave  <<EOF

options(warn=-1)

rm(list = ls(all = TRUE))

numsnps <- $numsnps
matrix(scan('ld.ld',quiet=T),nc=numsnps) -> co

#check that co is positive definite. Make diagonals 1.0001
library(corpcor)

if(is.positive.definite(co)==F){
	co <- make.positive.definite(co)
}
if(is.positive.definite(co)==F){
	matrix(scan('ld.ld',quiet=T),nc=numsnps) -> co
	for(i in 1:numsnps){
		co[i,i] <- 1.0001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.01
	}
}

library(CompQuadForm)

## now do the svd and CDF calculation
decomp <- svd(co)[1]
lambda <- decomp$d
lambda = as.numeric(unlist(lambda))
read.table('tempgene',header=F) -> gene
ass <- gene[,6]
sum(ass,na.rm=T) -> sumass
sum(!is.na(ass)) -> len # get number of snps with non-missing test stat
length(ass) -> asslen
ans = imhof(sumass[1], lambda, h = rep(1, length(lambda)), delta = rep(0,length(lambda)), epsabs = 10^(-16), epsrel = 10^(-16), limit = 1000)
empp <- ans$Qq
if (empp[1] < 0) {
	empp[1] = 2.22e-16 ## to correct for underflow errors, set to machine precision of R
}
snpp <- pchisq(max(ass),1,lower.tail=F)
snpname <- gene[which(ass==max(ass)),1][1]

write.table(data.frame('$chr','$gene',numsnps,unique(gene[4]),unique(gene[5]),sumass,empp,snpname,snpp),'alwayswritep',row.names=F,col.names=F, quote = F)

EOF
");
		open(CORRECTEDP, "alwayswritep");
		$correctedp = <CORRECTEDP>;
		open(PRINTCORRECTEDP, ">>gene-basedoutput.out");
		print PRINTCORRECTEDP "$correctedp";
}

sub sort_by_number{ $b <=> $a }

sub plink{ # Use plink to generate ld matrix, then check for errors
	
	if(defined($custom)){
		system("plink --bfile custom$chr --extract tempgene.snp --matrix --r --noweb --silent --out ld > /dev/null");
	}
	elsif(!defined($ldfile)){
		system("plink --bfile $chip/genebeds/$gene --extract tempgene.snp --matrix --r --noweb --silent --out ld > /dev/null");
	}
	

	open(PLINKLD, "ld.ld");
	@plinkld = <PLINKLD>;
	foreach $plinkld(@plinkld){
		$plinkld =~ s/nan/0/g; #replace nan's with 0 in ld matrix
		$plinkld =~ s/NaN/0/g;
	}
	if(scalar(@plinkld == 1)){
		@plinkld = 1;
	}
	
	open(PRINTPLINKLD, ">ld.ld");
	print PRINTPLINKLD "@plinkld";
	$numsnps = scalar(@plinkld);
	close PRINTPLINKLD;
	close PLINKLD;
	
}

sub customtest{ # using custom set of SNPs - default
	print "\nReading custom genotypes...done\n\n";
	#open(HGLIST,"../glist-hg18") or die "Error: Cannot find glist-hg18\n";
	open(HGLIST,"../glist-hg19") or die "Error: Cannot find glist-hg19\n";
	@hglist = <HGLIST>;
	close(HGLIST);

	for ($chr = 1; $chr <= 23; $chr++) {
		system("plink --bfile ../$custom --chr $chr --make-bed --out custom$chr --noweb --silent > /dev/null");
		if ($chr==23) {$chr ='X';}
		open(ALLGENES, "$chip/geneset/allgene$chr") or die "Error: Cannot find $chip/geneset/allgene$chr\n";
		@allgenes = <ALLGENES>;
		#close(ALLGENES);
		if ($chr=='X') {$chr =23;}
		
		foreach $gene (@allgenes){

			if(-e "tempgene.snp"){system("rm tempgene.snp");}
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			
			chomp($gene);
			#system("grep -w \"$gene\" ../glist-hg18 > gene.position");
			system("grep -w \"$gene\" ../glist-hg19 > gene.position");
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/\s+/,@glistpos[0]);
		
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
	
			system("plink --bfile custom$chr --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null"); 

			if(-e "plink.snplist"){
				system("mv plink.snplist tempgene.snp");
			}
			
			if (-z "tempgene.snp"){next;}

			system("grep -F -w -f tempgene.snp temppvals > tempgene.pvalue");	
			
			&maketempgener;
			&hapsims;
		}
	}
}

sub docustomchr { # using custom set of SNPs - chromosome
	print "\nReading custom genotypes...done\n\n";
	
	#open(HGLIST,"../glist-hg18") or die "Error: Cannot find glist-hg18\n";
	open(HGLIST,"../glist-hg19") or die "Error: Cannot find glist-hg19\n";
	@hglist = <HGLIST>;
	#close(HGLIST);

	for ($chr=$dochr){
		system("plink --bfile ../$custom --chr $chr --make-bed --out custom$chr --noweb --silent > /dev/null");
		if ($chr==23) {$chr ='X';}
		open(ALLGENES, "$chip/geneset/allgene$chr");
		@allgenes = <ALLGENES>;
	#	close(ALLGENES);
		if ($chr=='X') {$chr =23;}

		foreach $gene (@allgenes){
			if(-e "tempgene.snp"){system("rm tempgene.snp");}
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			
			chomp($gene);
			#system("grep -w \"$gene\" ../glist-hg18 > gene.position");
			system("grep -w \"$gene\" ../glist-hg19 > gene.position");
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/\s+/,@glistpos[0]);
		
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
	
			system("plink --bfile custom$chr --chr $chr --from-bp $start --to-bp $stop --write-snplist --noweb --silent > /dev/null"); 
			if(-e "plink.snplist"){system("mv plink.snplist tempgene.snp");}
			
			if (-z "tempgene.snp"){next;}

			system("grep -F -w -f tempgene.snp temppvals > tempgene.pvalue");	
			
			&maketempgener;
			&hapsims;
		}
	}
}

sub ldfiletest{ # using custom LD file - default
	
	#open(HGLIST,"../glist-hg18") or die "Error: Cannot find glist-hg18\n";
	open(HGLIST,"../glist-hg19") or die "Error: Cannot find glist-hg19\n";
	@hglist = <HGLIST>;
	close(HGLIST);

	for ($chr = 1; $chr <= 23; $chr++) {	
		system("awk '\$1==$chr {print \$0}' ../glist-hg19 > allgene$chr\n"); ## find all genes on chromosome
		if ($chr==23) {$chr =='X';}
		open(ALLGENES, "allgene$chr") or die "Error: Cannot find allgene$chr\n";
		#open(ALLGENES, "$chip/geneset/allgene$chr") or die "Error: Cannot find $chip/geneset/allgene$chr\n";
		@allgenes = <ALLGENES>;
		close(ALLGENES);
		
		system("awk '\$1==$chr {print \$0}' $ldfile > temptempld"); ## subset ld file by chromosome

		foreach $gene (@allgenes){

			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			if(-e "templd"){system("rm templd");}
			if(-e "ld.ld"){system("rm ld.ld");}
			
			chomp($gene);
			@items = split(/\s+/,$gene);
			$gene = $items[3];
			#system("grep -w \"$gene\" ../glist-hg18 > gene.position");
			system("grep -w \"$gene\" ../glist-hg19 > gene.position");
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/\s+/,@glistpos[0]);
			
			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
			## find all snps within gene using temptempld file only
			system("

R --vanilla --slave  <<EOF

options(warn=-1)
rm(list = ls(all = TRUE))

ld = read.table('temptempld', header = F, colClasses=c('integer','integer', 'character', 'integer', 'integer', 'character', 'numeric'))
names(ld) <- c('CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R')
ld = ld[which(ld[,2]>=$start & ld[,2]<=$stop),]
ld = ld[which(ld[,5]>=$start & ld[,5]<=$stop),]

if (nrow(ld) >=1) {
	snps = unique(c(ld[,3], ld[,6]))
	temppvals = read.table('temppvals', header = F, colClasses = c('character', 'numeric'))
	snpps = temppvals[which(temppvals[,1] %in% snps),]
	snps <- snpps[,1]
	if (nrow(snpps)==1) {
		write.table(snpps, 'tempgene.pvalue', row.names=F,col.names=F, quote = F)
		system('echo 1 > ld.ld')
	}
	if (nrow(snpps) > 1) {
		write.table(snpps, 'tempgene.pvalue', row.names=F,col.names=F, quote = F)
		ld = ld[,c(3,6,7)]
		dummy = as.data.frame(cbind(snpps[,1], snpps[,1], rep(1.0, nrow(snpps))))
		names(dummy)<- names(ld)
		ld = rbind(ld, dummy)
		ldmat = reshape(ld, direction='wide', timevar = 'SNP_B', idvar = 'SNP_A')
		ldmat = ldmat[order(match(ldmat[,1], snpps[,1])),] #change row order
		rownames(ldmat) = ldmat[,1]
		ldmat = ldmat[,-1]
		ldmat = ldmat[,order(match(names(ldmat), paste('R.', snpps[,1], sep = '')))]# change column order to make (upper) triangular matrix
		ldmat[is.na(ldmat)] <- 0.0
		ldmat[ldmat == 'NaN'] <- 0.0
		ldmat = apply(ldmat, 1, as.numeric)
		if (isSymmetric(unname(as.matrix(ldmat))) == FALSE) {
  			ldmat = ldmat+ t(ldmat)
  			diag(ldmat) <- 1.0
		}
		write.table(ldmat, 'ld.ld', row.names=F,col.names=F, quote = F)
	}
}
EOF
");

			if (-e "tempgene.pvalue"){
				if (-e "ld.ld") {
					&maketempgener;
					&hapsims;
				}
			}
			else {next;}
		}
	}
}

sub doldchr{ # using custom LD file - for one chromosome
	
	
	#open(HGLIST,"../glist-hg18") or die "Error: Cannot find glist-hg18\n";
	open(HGLIST,"../glist-hg19") or die "Error: Cannot find glist-hg19\n";
	@hglist = <HGLIST>;
	close(HGLIST);

	for ($chr = $dochr) {
		system("awk '\$1==$chr {print \$0}' ../glist-hg19 > allgene$chr\n"); ## find all genes on chromosome
		if ($chr==23) {$chr =='X';}
		open(ALLGENES, "allgene$chr") or die "Error: Cannot find allgene$chr\n";
		#open(ALLGENES, "$chip/geneset/allgene$chr") or die "Error: Cannot find $chip/geneset/allgene$chr\n";
		@allgenes = <ALLGENES>;
		close(ALLGENES);
		
		system("awk '\$1==$chr {print \$0}' $ldfile > temptempld"); ## subset ld file by chromosome

		foreach $gene (@allgenes){
			if(-e "tempgene.pvalue"){system("rm tempgene.pvalue");}
			if(-e "gene.position"){system("rm gene.position");}
			if(-e "templd"){system("rm templd");}
			if(-e "ld.ld"){system("rm ld.ld");}

			chomp($gene);
			@items = split(/\s+/,$gene);
			$gene = $items[3];
			#system("grep -w \"$gene\" ../glist-hg18 > gene.position");
			system("grep -w \"$gene\" ../glist-hg19 > gene.position");
			@glistpos = grep(/\b$gene\b/, @hglist);	
			chomp(@glistpos);
			@glistpos = split(/\s+/,@glistpos[0]);

			if(!defined($lower)){
				$start = @glistpos[1]-50000;
			}
			else{
				$start = @glistpos[1]-$lower;
			}
			if(!defined($upper)){
				$stop = @glistpos[2]+50000;
			}
			else{
				$stop = @glistpos[2]+$upper;
			}
			## find all snps within gene using temptempld file only
			system("

R --vanilla --slave  <<EOF

options(warn=-1)
rm(list = ls(all = TRUE))

ld = read.table('temptempld', header = F, colClasses=c('integer','integer', 'character', 'integer', 'integer', 'character', 'numeric'))
names(ld) <- c('CHR_A','BP_A','SNP_A','CHR_B','BP_B','SNP_B','R')
ld = ld[which(ld[,2]>=$start & ld[,2]<=$stop),]
ld = ld[which(ld[,5]>=$start & ld[,5]<=$stop),]

if (nrow(ld) >=1) {
	snps = unique(c(ld[,3], ld[,6]))
	temppvals = read.table('temppvals', header = F, colClasses = c('character', 'numeric'))
	snpps = temppvals[which(temppvals[,1] %in% snps),]
	snps <- snpps[,1]
	if (nrow(snpps)==1) {
		write.table(snpps, 'tempgene.pvalue', row.names=F,col.names=F, quote = F)
		system('echo 1 > ld.ld')
	}
	if (nrow(snpps) > 1) {
		write.table(snpps, 'tempgene.pvalue', row.names=F,col.names=F, quote = F)
		ld = ld[,c(3,6,7)]
		dummy = as.data.frame(cbind(snpps[,1], snpps[,1], rep(1.0, nrow(snpps))))
		names(dummy)<- names(ld)
		ld = rbind(ld, dummy)
		ldmat = reshape(ld, direction='wide', timevar = 'SNP_B', idvar = 'SNP_A')
		ldmat = ldmat[order(match(ldmat[,1], snpps[,1])),] #change row order
		rownames(ldmat) = ldmat[,1]
		ldmat = ldmat[,-1]
		ldmat = ldmat[,order(match(names(ldmat), paste('R.', snpps[,1], sep = '')))]# change column order to make (upper) triangular matrix
		
		ldmat[is.na(ldmat)] <- 0.0
		ldmat[ldmat == 'NaN'] <- 0.0
		ldmat = apply(ldmat, 1, as.numeric)
		if (isSymmetric(unname(as.matrix(ldmat))) == FALSE) {
  			ldmat = ldmat+ t(ldmat)
  			diag(ldmat) <- 1.0
		}
		write.table(ldmat, 'ld.ld', row.names=F,col.names=F, quote = F)
	}
}
EOF
");

			if (-e "tempgene.pvalue"){
				if (-e "ld.ld") {
					&maketempgener;
					&hapsims;
				}
			}
			else {next;}
		}	
	}
}


sub maketempgener{
unless(-z "tempgene.pvalue"){
	system("

R --vanilla --slave <<EOF
options(warn=-1)
pvals <- read.table('tempgene.pvalue',header=F, colClasses = c('character', 'numeric'))
pvals = pvals[order(pvals[,1]),]
positions <- read.table('gene.position',header=F)
tempgene <- data.frame(pvals[,1],rep(positions[,4],length(pvals[,1])),rep(positions[,1],length(pvals[,1])),rep(positions[,2],length(pvals[,1])),rep(positions[,3],length(pvals[,1])),qchisq(pvals[,2],1,lower.tail=F))
write.table(tempgene,'tempgene',row.names=F,col.names=F,quote=F)

EOF");

system("awk '{print \$1;}' tempgene > tempgene.snp");
}
}
	
sub checkpackages{
	system("
R --vanilla --slave <<EOF
options(warn=-1)
write.table(.packages(all.available=T),'R-packagelist',quote=F,col.names=F,row.names=F)
EOF");

	system("grep -w \"CompQuadForm\" R-packagelist >> testpackagelist");
	system("grep -w \"corpcor\" R-packagelist >> testpackagelist");
	open(TESTPACKAGELIST,"testpackagelist");
	my @testpackagelist = <TESTPACKAGELIST>;
	if(scalar(@testpackagelist ne 2)){
		die "Error: Missing R packages: PEGASUS requires corpcor and CompQuadForm to run\n";
	}
}

sub numbercheck{
	system("
R --vanilla --slave <<EOF
options(warn=-1)
isnumeric <- function(s) !is.na(as.numeric(s))
pvals <- read.table('temppvals',header=F,colClasses=c('character','numeric'))
ifelse(isnumeric(pvals[,2]),0,1) -> isnumber
if(sum(isnumber) > 0){
	write.table(pvals[which(isnumber==1),],'temp',row.names=F,col.names=F,quote=F)
}
ifelse(pvals[,2]<0 | pvals[,2] > 1,1,0) -> notpvalue
if(sum(notpvalue) > 0 | is.na(sum(notpvalue))){
	write.table(pvals[which(notpvalue==1),],'invalid-pvalues',row.names=F,col.names=F,quote=F)
}
EOF");

	if(-e "temp"){
		system("cat temp >> invalid-pvalues");
	}

	if(-e "invalid-pvalues"){
		system("mv invalid-pvalues ../invalid-pvalues");
		die "Error: Invalid p-values detected. See \"invalid-pvalues\" for list of problem SNPs\n";
	}
}