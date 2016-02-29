#!/usr/bin/perl -w
#
# Script for setting up a MySQL database
# structure which will hold resequencing
# data.
#
# By: Magnus Isaksson 2008
#

use strict;

use Getopt::Long;
use File::Spec;
use DBI;

# Include custom libs in the path
use FindBin;
#use lib $FindBin::Bin."/../lib";

#use Bio::OlinkGenomics::Utils;
#use Log::Log4perl qw/:easy/;

# Subroutine prototypes
sub usage;
sub tableExist;
sub createNewTable;

# Globals
my ($sample_ID, $dsn, $db_user, $db_passwd, $log4perl_file, $output_dir);

# Parse the command line
GetOptions(	
	"s=s" => \$sample_ID,
	"o=s" => \$output_dir,
	"d=s" => \$dsn,
	"u=s" => \$db_user,	
	"p=s" => \$db_passwd
) || &usage();

# Missing options
if ( !($sample_ID) || !($dsn) || !($db_user) || !($db_passwd)) { usage(); }

# For the logfile
my $DEFAULT_LOG4PERL_FILE = $FindBin::Bin . '/../config/log4perl.properties';

# initialize log for perl
#$log4perl_file = File::Spec->rel2abs($DEFAULT_LOG4PERL_FILE);
#Log::Log4perl->init($log4perl_file);
#my $logger = get_logger();


# Open connection to database.
my $dbh = DBI->connect($dsn, $db_user, $db_passwd);

	# Check if table exist.
	if (&tableExist() eq "Yes") {
		$dbh->disconnect();
		#$logger->error( "Table $sample_ID exists in database! We are using settings: ".$dsn.".");
		print "Table $sample_ID exists in database! We are using settings: ".$dsn.".";
		exit(1);
	} else {		
		# Make new table.
		&createNewTable();
	}
	
	print "Successfully created table ".$sample_ID." with database settings: ".$dsn.".\n";
	#$logger->debug("Successfully created table ".$sample_ID." with database settings: ".$dsn.".");
	
# Close connection to database.
$dbh->disconnect();


# Sub creates new table in database.
sub createNewTable {
	
	my $db_settings;
	
	if ( defined($output_dir) ) { 
		$db_settings = "ENGINE=MyISAM DEFAULT CHARSET=latin1 DATA DIRECTORY='".$output_dir."'"; 
	} else {
		$db_settings = "ENGINE=MyISAM DEFAULT CHARSET=latin1";
	}
	my $sth = $dbh->prepare(qq{ 
		CREATE TABLE `$sample_ID` (
			`db_id` int(10) unsigned NOT NULL auto_increment,
			`id` char(20) NOT NULL,
			`start` int(10) unsigned NOT NULL,
			`end` int(10) unsigned NOT NULL,
			`ref_id` char(24) NOT NULL,
			`ref_start` int(10) unsigned NOT NULL,
			`ref_end` int(10) unsigned NOT NULL,
			`ref_strand` enum('+','-') NOT NULL default '+',
			`score` tinyint(3) unsigned NOT NULL,
			`hits` smallint(5) unsigned NOT NULL,
			PRIMARY KEY  (`db_id`)
		) $db_settings; 
	});
	
	if (!($sth->execute())) {
		#$logger->error("Unable to create table ".$sample_ID." with database settings: ".$dsn."!");
		$dbh->disconnect();
		exit(1);
	}
	
	$sth->finish();
}

# Sub checks if table exists in database.
sub tableExist {
	
	my $sth = $dbh->prepare(qq{ 
		SHOW TABLES LIKE '$sample_ID'; 
	});
	$sth->execute(); 
	my ($exists) = $sth->fetchrow_array();
	$sth->finish();

	if ($exists) {
		$exists = "Yes";
	} else {
		$exists = "No";
	}

	return $exists;
}

# Sub to tell you how to use this script.
sub usage {
  print "\nUsage: $0 -s <sample_id> -o <db_data_directory> -d <Perl_DBI_string> -u <db_user> -p <db_password>\n 
 -s Sample id name, will be name of database table.
 -o Absolut directory where database data will be located (default is MySQL settings).
 -d Perl DBI string example DBI:mysql:SOLiD:127.0.0.1.
 -u Database user name.
 -p Database password.\n\n";
  exit(1);
}
