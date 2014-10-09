#!/usr/bin/env perl
package CAP;
use strict;
use warnings;
#use Config::Simple; # for kbase config


use JSON;

use Data::Dumper;


use AWE::Workflow;
use AWE::Task;
use AWE::TaskInput;
use AWE::TaskOutput;





sub new {
	my ($class, %h) = @_;
	if (defined($h{'shocktoken'}) && $h{'shocktoken'} eq '') {
		$h{'shocktoken'} = undef;
	}
	my $self = {
		aweserverurl	=> $ENV{'AWE_SERVER_URL'} ,
		shockurl	=> $ENV{'SHOCK_SERVER_URL'},
		clientgroup	=> $ENV{'AWE_CLIENT_GROUP'},
		shocktoken	=> $h{'shocktoken'}
	};
	bless $self, $class;
#	$self->readConfig();
	return $self;
}



sub aweserverurl {
	my ($self) = @_;
	return $self->{'aweserverurl'};
}
sub shockurl {
	my ($self) = @_;
	return $self->{'shockurl'};
}
sub clientgroup {
	my ($self) = @_;
	return $self->{'clientgroup'};
}
sub shocktoken {
	my ($self, $value) = @_;
	if (defined $value) {
		$self->{'shocktoken'} = $value;
	}
	return $self->{'shocktoken'};
}
sub readConfig {
	my ($self) = @_;
	my $conf_file = $ENV{'KB_TOP'}.'/deployment.cfg';
	unless (-e $conf_file) {
		die "error: deployment.cfg not found ($conf_file)";
	}
	my $cfg_full = Config::Simple->new($conf_file );
	my $cfg = $cfg_full->param(-block=>'AmethstService');
	unless (defined $self->{'aweserverurl'} && $self->{'aweserverurl'} ne '') {
		$self->{'aweserverurl'} = $cfg->{'awe-server'};
		unless (defined($self->{'aweserverurl'}) && $self->{'aweserverurl'} ne "") {
			die "awe-server not found in config";
		}
	}
	unless (defined $self->{'shockurl'} && $self->{'shockurl'} ne '') {
		$self->{'shockurl'} = $cfg->{'shock-server'};
		unless (defined(defined $self->{'shockurl'}) && $self->{'shockurl'} ne "") {
			die "shock-server not found in config";
		}
	}
	unless (defined $self->{'clientgroup'} && $self->{'clientgroup'} ne '') {
		$self->{'clientgroup'} = $cfg->{'clientgroup'};
		unless (defined($self->{'clientgroup'}) && $self->{'clientgroup'} ne "") {
			die "clientgroup not found in config";
		}
	}
}

sub create_cap_workflow {
	my ($self, $input_ref) = @_;

	my $workflow = new AWE::Workflow(
		"pipeline"=> "cap",
		"name"=> "cap",
		"project"=> "cap",
		"user"=> "kbase-user",
		"clientgroups"=> $self->clientgroup,
		"noretry"=> JSON::true
	);

	
	
	
	
	#task 1  (khmer)
	# input: *renamed
	#generates simulated metagenomes combined into one file
	#cat *renamed > metagenome-cumulative.fa
	#requires khmer and dependencies
	#python ~/khmer/scripts/normalize-by-median.py -k 20 -C 10 -N 4 -x 3e8 -s norm10k20.kh metagenome-cumulative.fa
	#rm norm10k20.kh #do not need
	#rm metagenome-cumulative.fa #do not need
	# output: ????
	
	
	#task 2 (velvet)
	#input: metagenome-cumulative.fa.keep
	#requires velvet install
	#velveth assembly 21 metagenome-cumulative.fa.keep
	#velvetg assembly
	# output: assembly/contigs.fa
	
	#task 3 (cap)
	#input: contigs.fa
	#can be done anytime after assembly but prior to bedtools
	#python coverage-bed-reference.py contigs.fa #produces contigs.fa.bed
	#output: contigs.fa.bed
	
	
	#task4 (bowtie)
	#input: contigs.fa
	#bowtie2-build contigs.fa assembly
	#output: ????
	
	
	
	
	#taskgroup 5 (bowtie)
	#can be done in parallel
	#input: *renamed
	#for x in *renamed; do bowtie2 -x assembly -f $x -S $x.sam; done
	#output: *renamed.sam
	
	
	#taskgroup 6 (cap)
	#requires samtools, can be done in parallel
	#input: *sam
	#for x in *sam; do samtools view -b -t contigs.fa $x > $x.bam; done
	#output: *sam.bam
	
	
	#taskgroup 7 (cap)
	#requires bedtools, can be done in parallel
	#input: *bam
	#for x in *bam; do bedtools bamtobed -i $x > $x.bed; done
	#output: *.bed
	
	
	#taskgroup 8 (cap)
	#requires preceding completed
	#input:  meta*bed?
	#for x in meta*bed; do coverageBed -a $x -b contigs.fa.bed > $x.reads.mapped; done
	#output: .reads.mapped
	
	#????alternate for preceeding, user dependent
	#????for x in meta*bed; do coverageBed -a $x -b contigs.fa.bed -d > $x.bed.coverage.perbase; done
	#requires all bedfiles completed, can be run in parallel but on the same computer, ask andreas about dict structures
	
	
	#taskgroup 9 (cap)
	#input: *mapped
	#for x in *mapped; do python get-rpkm.py $x; done
	#output: *rpkm
	
	
	
	#task 10 (cap)
	#requires all rpkm calculated
	#input: *rpkm
	#python merge.py *rpkm
	
	#requires curl (in cap)
	#curl "http://api.metagenomics.anl.gov/1/annotation/similarity/mgm4566339.3?type=ontology&source=Subsystems" > annotations.txt
	
	#python best-hit.py annotations.txt
	#requires R and dependencies phyloseq, plyr, ggplot, saves output as RData
	#also requires a meta.txt file
	#R < core.R --vanilla
	
	#output: ???
	
	my $json = JSON->new;
	
	print "AWE workflow:\n".$json->pretty->encode( $workflow->getHash() )."\n";
	
}


my $wf = new CAP();
$wf->create_cap_workflow();


