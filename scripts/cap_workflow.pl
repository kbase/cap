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


sub other_stuff {
	
	
	#task 1  (khmer)
	# input: *renamed
	#generates simulated metagenomes combined into one file
	#cat *renamed > metagenome-cumulative.fa
	#requires khmer and dependencies
	#python ~/khmer/scripts/normalize-by-median.py -k 20 -C 10 -N 4 -x 3e8 -s norm10k20.kh metagenome-cumulative.fa
	#rm norm10k20.kh #do not need
	#rm metagenome-cumulative.fa #do not need
	# output: metagenome-cumulative.fa.keep

	
	
	#task 2 (velvet)
	#input: metagenome-cumulative.fa.keep
	#requires velvet install
	#velveth assembly 21 metagenome-cumulative.fa.keep
	#velvetg assembly
	# output: assembly/contigs.fa
	

	
}


sub shock_resource {
	my ($host, $node, $filename) = @_;
	
	unless (defined $node) {
		#http://shock.metagenomics.anl.gov/node/80cce328-f8ce-4f2c-b189-803ac12f9e44
		my ($h, $n) =$host =~ /^(.*)\/node\/(.*)(\?download)?/;
		
		unless (defined $n) {
			die "could not parse shock url";
		}
		$host = $h;
		$node = $n;
		
	}
	
	
	my $res = {"resource" => "shock",
				"host" => "$host",
				"node" => "$node"};
	
	if (defined $filename) {
		$res->{"filename"} = $filename ;
	}
	return $res;
}

sub task_resource {
	my ($task, $pos) = @_;
	
	my $res = {"resource" => "task",
		"task" => $task, #string
		"position" => $pos # number;
	};
	
	return $res;
}



sub create_cap_workflow {
	my ($self, $assembly, $metatxt, $mgmid, $input_ref) = @_;

	my $workflow = new AWE::Workflow(
		"pipeline"=> "cap",
		"name"=> "cap",
		"project"=> "cap",
		"user"=> "kbase-user",
		"clientgroups"=> $self->clientgroup,
		"noretry"=> JSON::true
	);

	
	
	my $h = $self->shockurl;
	
	##############################
	# files, assemblu and mgm id , meta.txt
	
	#task 1 (cap)
	#input: contigs.fa
	#can be done anytime after assembly but prior to bedtools
	#python coverage-bed-reference.py contigs.fa #produces contigs.fa.bed
	#output: contigs.fa.bed
	
	my $newtask = $workflow->newTask('app:CAP.coverage-bed-reference.default');
	$newtask->{'cmd'}->{'app_args'} = [shock_resource($assembly)];
	my $t1 = $newtask->taskid();
	
	
	#task2 (bowtie)
	#input: contigs.fa
	#bowtie2-build contigs.fa assembly
	#output: assembly.*
	
	$newtask = $workflow->newTask('app:Bowtie2.bowtie2-build.default');
	$newtask->{'cmd'}->{'app_args'} = [task_resource($t1, 0)];
	my $t2 = $newtask->taskid();
	
	#taskgroup 3 (bowtie)
	#can be done in parallel
	#input: *renamed, assembly.*
	#for x in *renamed; do bowtie2 -x assembly -f $x -S $x.sam; done
	#output: *renamed.sam
	
	my @taskgroup3 = ();
	for (my $i = 0 ; $i < @{$input_ref} ; $i++) {
		$taskgroup3[$i] = $workflow->newTask('app:Bowtie2.bowtie2.default');
		$taskgroup3[$i]->{'cmd'}->{'app_args'} = [	shock_resource($input_ref->[$i]),
											task_resource($t2, 0), task_resource($t2, 1), task_resource($t2, 2), task_resource($t2, 3), task_resource($t2, 4), task_resource($t2, 5)];# this line is bowtie database
	}
	
	
	#taskgroup 4 (cap)
	#requires samtools, can be done in parallel
	#input: *sam
	#for x in *sam; do samtools view -b -t contigs.fa $x -o $x.bam; done
	#output: *sam.bam
	my @taskgroup4 = ();
	for (my $i = 0 ; $i < @{$input_ref} ; $i++) {
		$taskgroup4[$i] = $workflow->newTask('app:Samtools.samtools.default');
		$taskgroup4[$i]->{'cmd'}->{'app_args'} = [ task_resource($taskgroup3[$i]->taskid(), 0) , shock_resource($assembly)];
	}
	
	
	
	#taskgroup 5 (cap)
	#requires bedtools, can be done in parallel
	#input: *bam
	#for x in *bam; do bedtools bamtobed -i $x > $x.bed; done
	#output: *.bed
	my @taskgroup5 = ();
	for (my $i = 0 ; $i < @{$input_ref} ; $i++) {
		$taskgroup5[$i] = $workflow->newTask('app:Bedtools.bedtools.bamtobed');
		$taskgroup5[$i]->{'cmd'}->{'app_args'} = [ task_resource($taskgroup4[$i]->taskid(), 0)];
	}
	
	
	#taskgroup 6 (cap)
	#requires preceding completed
	#input:  meta*bed, contigs.fa.bed(3)
	#for x in meta*bed; do coverageBed -a $x -b contigs.fa.bed > $x.reads.mapped; done
	#output: .reads.mapped
	my @taskgroup6 = ();
	for (my $i = 0 ; $i < @{$input_ref} ; $i++) {
		$taskgroup6[$i] = $workflow->newTask('app:Bedtools.coverageBed.default');
		$taskgroup6[$i]->{'cmd'}->{'app_args'} = [ task_resource($taskgroup5[$i]->taskid(), 0)];
	}
	
	
	#????alternate for preceeding, user dependent
	#????for x in meta*bed; do coverageBed -a $x -b contigs.fa.bed -d > $x.bed.coverage.perbase; done
	#requires all bedfiles completed, can be run in parallel but on the same computer, ask andreas about dict structures
	
	
	#taskgroup 7 (cap)
	#input: *mapped
	#for x in *mapped; do python get-rpkm.py $x; done
	#output: *rpkm
	my @taskgroup7 = ();
	for (my $i = 0 ; $i < @{$input_ref} ; $i++) {
		$taskgroup7[$i] = $workflow->newTask('app:CAP.get-rpkm.default');
		$taskgroup7[$i]->{'cmd'}->{'app_args'} = [ task_resource($taskgroup6[$i]->taskid(), 0)];
	}

	
	
	#task 8 (cap)
	#requires all rpkm calculated
	#input: *rpkm, meta.txt
	#python merge.py *rpkm
	
	#requires curl (in cap)
	#curl "http://api.metagenomics.anl.gov/1/annotation/similarity/mgm4566339.3?type=ontology&source=Subsystems" > annotations.txt
	
	#python best-hit.py annotations.txt
	#requires R and dependencies phyloseq, plyr, ggplot, saves output as RData
	#also requires a meta.txt file
	#R < core.R --vanilla
	
	#output: metag.RData
	
	$newtask = $workflow->newTask('app:CAP.final.default');
	$newtask->{'cmd'}->{'app_args'} = [shock_resource($metatxt), $mgmid, 'RPKM files here, from taskgroup7'];

	
	
	
	my $json = JSON->new;
	
	print "AWE workflow:\n".$json->pretty->encode( $workflow->getHash() )."\n";
	
}


my $test = 'http://shock.metagenomics.anl.gov/node/80cce328-f8ce-4f2c-b189-803ac12f9e44';

my $wf = new CAP();
$wf->create_cap_workflow($test, $test, "some_mgm_id", [$test, $test]);


