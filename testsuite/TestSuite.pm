package TestSuite;

##
## Methods for the TiberCAD test suite
##


use strict;
use warnings;


BEGIN {

  use Exporter ();

  our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);

  $VERSION = 1.0;

  @ISA = qw(Exporter);
  @EXPORT = qw();
  %EXPORT_TAGS = ();

  @EXPORT_OK = qw(run_test
                  check_compilation
                  read_configuration
                  get_option
                  set_option
                  send_email
                  cleanup_test_dir);

}
our @EXPORT_OK;


use Cwd;
use Text::ParseWords;
use Term::ANSIColor;
use File::Basename;
use Net::SMTP;


##
## some module-global variables
##

our $makefile;
our $make;
our $make_run;
our $config;
our $outdir;
our $refdir;
our $to_dev_null;
our $to_log_file;
our $textwidth;
our $padchar;
our $default_precision;
our @known_suffixes;
our %options;



%options = (
    'verbose' => 0,
    'testonly' => 0,
    'smtp_server' => "localhost",
    'smtp_port' => 25,
    'smtp_user' => '',
    'smtp_passwd' => '',
    'mail_sender' => "$ENV{USER}",
    'package' => "unknown",
    'logfile' => "run.log"
    );

$make_run = "make run";
$make = "make";
$makefile = "Makefile";
$config = "test.conf";
$outdir = "output";
$refdir = "reference";
$to_dev_null = " 2> /dev/null > /dev/null";
$to_log_file = " 2> $options{'logfile'} > $options{'logfile'}";
$textwidth = 50;
$padchar = ".";
$default_precision = 1e-4;
@known_suffixes = ('dat', 'gmv', 'vtu', 'xyz');



##
## Read the configuration file
##
## The argument is the name of the config file
##
sub read_configuration($);


##
## Do a compilation and check for success
##
## Takes directory as argument
## Returns 1 in case of error, 0 otherwise
##
sub check_compilation($);


##
## Run the tests in a given directory
##
## Takes directory as argument
## Returns 1 in case of error, 0 otherwise
##
sub run_test($);


##
## Get the named option
##
sub get_option($);


##
## Set an option
## 
## First argument is the option name, second is value
##
sub set_option($$);


##
## Send an email with status report
##
## Argument is a hash reference containing the mail content:
##   "message" => the message,
##   "subject" => subject
##
sub send_email($);

##
## Run the executable
##
## Returns 1 in case of error, 0 otherwise
##
sub run_executable($);


##
## Cleans up the directory
##
sub cleanup_test_dir($);


##
## Test the output data against the reference data
##
## Returns 1 in case of error, 0 otherwise
##
sub compare_results_with_reference($$$);


##
## Extract data fields from a .dat data file
## (usually created in sweeps and for grace)
##
## Returns the data in a list
##
sub extract_data_from_characteristic($$);


##
## Extract data fields from a .gmv file
##
## Returns the data in a list
##
sub extract_data_from_gmv($$);


##
## Extract data fields from a .vtu file
##
## Returns the data in a list
##
sub extract_data_from_vtu($$);

##
## Extract data fields from a .xyz file
##
## Returns the data in a list
##
sub extract_data_from_xyz($$);


##
## Compare numeric data in two lists with a certain precision
##
## returns 0 if the data lists are equal, 1 if there is a difference
##
sub compare_data_vectors($$$);


##
## Compare string data in two lists
##
## returns 0 if the string lists are equal, 1 if there is a difference
##
sub compare_string_vectors($$);


##
## Reads the check sections from the config file
##
sub read_checks_from_config($$);


##
## Trim leading and trailing white space
##
sub trim($);


##
## Print 'ok' or 'failure', depending on argument
##
sub print_result($);


##
## Return 'ok' or 'failure', depending on argument
##
sub sprint_result($);


##
## Pad with character
##
## First argument is the string to pad
## Second argument is the padding character
## Returns padded string
##
sub pad_to_textwidth($$);


##
## Split variable name and precision given a string of the form
## variable_name(precision)
##
sub split_variable($); 


##
## Tells level of verbosity
##
sub verbose;





###################################################################
## Implementations
###################################################################


sub read_configuration($) {

  my $config = $_[0];

  open(CF, "<", $config) || die "Configuration file invalid!";

  while (<CF>) {

    SWITCH: {
            
       /\s*#/ &&
         last SWITCH;

       /\s*topdir\s*=(.+)/ && do {
         my $topdir = trim($1);
         if ($topdir =~ /^\s*\/+/) {
           $options{'topdir'} = $topdir;
         }
         else {
           $options{'topdir'} = getcwd() . "/" . $topdir;
         }
         last SWITCH;
       };
       
       /\s*mail_recipient\s*=(.+)/ && do {
         $options{'mail_rcpt'} = trim($1);
         last SWITCH;
       };
       
       /\s*mail_sender\s*=(.+)/ && do {
         $options{'mail_sender'} = trim($1);
         last SWITCH;
       };
       
       /\s*smtp_server\s*=(.+)/ && do {
         $options{'smtp_server'} = trim($1);
         last SWITCH;
       };
      
       /\s*smtp_port\s*=(.+)/ && do {
         $options{'smtp_port'} = trim($1);
         last SWITCH;
       };
      
       /\s*smtp_user\s*=(.+)/ && do {
         $options{'smtp_user'} = trim($1);
         last SWITCH;
       };
      
       /\s*smtp_passwd\s*=(.+)/ && do {
         $options{'smtp_passwd'} = trim($1);
         last SWITCH;
       };

       /\s*verbose\s*=(.+)/ && do {
         $options{'verbose'} = trim($1);
         last SWITCH;
       };

       /\s*svn\s*=(.+)/ && do {
         $options{'svn'} = trim($1);
         last SWITCH;
       };

       /\s*package_name\s*=(.+)/ && do {
         $options{'package'} = trim($1);
         last SWITCH;
       };
    }
  }

  close(CF);

}



sub check_compilation($) {

  my $fail = 1;

  my $dir = getcwd;
  
  if (not -d  $_[0]) { return 1 };
  
  chdir($_[0]);

  print("\nchdir to ", getcwd(), "\n") if verbose();
  print(pad_to_textwidth("Check compilation", $padchar)) if verbose();

  $fail = system("$make $to_dev_null") if -e $makefile;
  
  print_result($fail) if verbose();
  
  chdir($dir);

  return $fail;
}


sub run_test($) {

  my $dir = $_[0];
  
  if (not -d $dir) {
    return 1;
  }

  if (not -e "$dir/$config") {
    # this is not assumed to be an error
    return 0;
  }

  my $olddir = getcwd;
  chdir($dir);
  print("\nRunning test in $dir ...\n") if verbose();

  my %checks = ();
  my %opts = ();
  # read the configuration
  read_checks_from_config(\%checks, \%opts) && return 0;

  my $command = "";
  (defined $opts{"commandline"}) && ($command = $opts{"commandline"});

  my $infile = "";
  (defined $opts{"inputfile"}) && ($infile = $opts{"inputfile"});

  my $out = $outdir;
  (defined $opts{"outputdir"}) && ($out = $opts{"outputdir"});

  my $ref = $refdir;
  (defined $opts{"referencedir"}) && ($ref = $opts{"referencedir"});
 
  my $failure = 0;
  if (not $options{'testonly'}) {
    if ($command ne "") {
      $failure = system("$command > /dev/null 2>&1");
    }
    else {
      $failure = run_executable($infile);
    }
  }
  $failure = compare_results_with_reference(\%checks, $out, $ref) unless $failure;

  chdir($olddir);

  return $failure;
}



sub run_executable($) {

  my $infile = $_[0];
  my $fail = 0;

  if (($infile ne "") && -f $infile) {
    print(pad_to_textwidth("Run", $padchar)) if verbose();
    $fail = system("$options{'topdir'}/bin/tibercad $infile $to_log_file");
  } 
  else {

    if (not -e $makefile) {
      return 1;
    } 

    print(pad_to_textwidth("Run", $padchar)) if verbose();

    $fail = system("$make_run $to_log_file");
  }

  print_result($fail) if verbose();

  return $fail;
}



sub cleanup_test_dir($) {

  system("$make -C $_[0] clean $to_dev_null");
}



sub compare_results_with_reference($$$) {

  my $result = 0;
  
  my $checks = $_[0];
  my $outputdir = $_[1];
  my $referencedir = $_[2];

  # the number of checks
  my $num_checks = scalar(keys(%$checks));
  
  print("Check data:\n") if verbose();
  
  while ( my ($file, $vars) = each(%$checks)) {
    my $datafile = "$outputdir/$file";
    my $reffile = "$referencedir/$file";

    my ($filename, $path, $suffix) = fileparse($file, @known_suffixes);
    if (! length($suffix)) {
      print "Filetype unknown: $file";
      return 1;
    }
    if (! -e $reffile ) {
      print "No such reference file: $reffile";
      return 1;
    }

    if ($#{$vars} < 0) {
      print(pad_to_textwidth("    $file (no variables to check)",
                                                 $padchar)) if verbose();

      print_result(not -e $datafile) if verbose();
      next;
    }
    print(pad_to_textwidth("    $file", " ")) if verbose();

    if (not -e $datafile) {
      $result = 1;
      print_result(1) if verbose();
      next;
    }
    print("\n") if ($#{$vars} >= 0 && verbose());


    foreach my $var (@{$vars}) {

      my ($variable, $prec, $type) = split_variable($var);

      print(pad_to_textwidth("      $variable", $padchar)) if verbose();
      my @data1 = ();
      my @data2 = ();

      # assign the right data extraction method
      my $extract_func = 0;

      SWITCH: for ($suffix) {
          /dat/   && do {
            $extract_func = \&extract_data_from_characteristic;
            last SWITCH;
          };

          /gmv/   && do {
            $extract_func = \&extract_data_from_gmv;
            last SWITCH;
          };

          /vtu/   && do {
            $extract_func = \&extract_data_from_vtu;
            last SWITCH;
          };

          /xyz/   && do {
            $extract_func = \&extract_data_from_xyz;
            last SWITCH;
          };

      }

      @data2 = &$extract_func($variable, $reffile);

      if ($#data2 < 0) {
        print "Variable $variable seems to not exist in $reffile\n";
        return 1;
      }

      @data1 = &$extract_func($variable, $datafile);

      my $fail = 0;
      if ($type eq 's')
      {
        $fail = compare_string_vectors(\@data1, \@data2);
      }
      else
      {
        $fail = compare_data_vectors(\@data1, \@data2, $prec);
      }

      if ($fail != 0) { $result = 1 };
      print_result($fail) if verbose();
    }
  }

  return $result;
}




sub extract_data_from_characteristic($$) {

  my ($variable, $file) = @_;

  open(SF, '<', $file);

  eof(SF) && print "Data file does not exist!\n";

  my @data = ();

  my $pos = 0;
  my $oldpos = 0;
  while (not eof) {
    $oldpos = $pos;
    $pos = tell;
    my $line = <SF>;
    trim($line);
    ($line !~ /#.*/) && last;
  }
  $pos = $oldpos;

  seek SF, $pos, 0;

  my @vars = (split(/\s+/, trim(<SF>)));
  if ($#vars < 0) { print "Data file has no data!\n" };

  my $index = 0;
  foreach my $item (@vars)
  {
    ($item eq $variable) && last;

    $index++;
  }

  # the first one is the hash
  $index--;

  if ($index == $#vars) { return @data };

  while (not eof)  {
    my @l = (split(/\s+/, trim(<SF>)));
    if ($#l > 0) { push @data, ($l[$index]) }; 
  }

  return @data;
}


sub extract_data_from_xyz($$) {

  my ($variable, $file) = @_;

  open(SF, '<', $file);

  eof(SF) && print "Data file does not exist!\n";

  my @data = ();

  my $pos = 0;
  my $oldpos = 0;
  while (not eof) {
    $oldpos = $pos;
    $pos = tell;
    my @line = (split(/\s+/, trim(<SF>)));
    ($#line >= 3) && last;
  }
  $pos = $oldpos;

  seek SF, $pos, 0;

  my $index = 0;
  foreach my $item ('species', 'x', 'y', 'z')
  {
    ($item eq $variable) && last;

    $index++;
  }

  if ($index == 4) { return @data };

  while (not eof)  {
    my @l = (split(/\s+/, trim(<SF>)));
    if ($#l > 0) { push @data, ($l[$index]) }; 
  }

  return @data;
}





sub extract_data_from_gmv($$) {

  my ($variable, $file) = @_;

  open(SF, '<', $file);

  eof(SF) && print "Data file does not exist!\n";


  # first, goto variable section
  while (not eof) {
    my $line = <SF>;
    ($line =~ /^variable/) && last;
  }

  # now, find the variable
  my @data = ();
  while (not eof) {
    my @line = (split(/\s+/, <SF>));
    if ($#line < 0) { next };

    if ($line[0] eq $variable) {
      @data = split(/\s+/, <SF>);
      last;
    }
  }

  return @data;
}






sub extract_data_from_vtu($$) {

  use XML::Parser;
  use MIME::Base64;

  my ($variable, $file) = @_;

  # variables visible to all XML parser handlers
  our $currvar = $variable;
  our $append = 0;
  our $base64 = 0;
  our $data = "";
  our @dataarray = ();


  # handler for starting tag
  sub startElement {
    my ($parseinst, $element, %attrs) = @_;

    if ($element =~ 'DataArray') {
      if ((exists $attrs{'Name'}) && ($attrs{'Name'} =~ $currvar)) {
        ($attrs{'format'} =~ 'binary') && ($base64 = 1);
        $append = 1;
      }
    }
  }


  # handler for ending tag
  sub endElement {
    my ($parseinst, $element) = @_;

    if ($append) {
      if ($base64) {
        # we have to decode and unpack
        my $decoded = decode_base64($data);
        my @unpacked = unpack("If*", $decoded);
        shift @unpacked;
        @dataarray = (@dataarray, @unpacked);
      }
      else {
        # we just split it up
        @dataarray = (@dataarray, split(/\s+/, trim($data)));
      }
    }

    # reset data container
    $data = "";
    $append = 0;
  }


  # handler for data
  sub characterData {
    my ($parseinst, $chardata) = @_;

    # the data might be read in pieces, so we put them together
    # as strings and extract data afterwards
    if ($append) {
      $data = $data . $chardata;
    }
  }


  # default handler
  sub default {
    #my ($parseinst, $data) = @_;
  }


  #
  # the function body
  #

  my $parser = new XML::Parser;
  $parser->setHandlers( Start => \&startElement,
                        End => \&endElement,
                        Char => \&characterData,
                        Default => \&default
                      );

  $parser->parsefile($file);

  return @dataarray;
}




sub compare_data_vectors($$$) {

  my @data1 = @{$_[0]};
  my @data2 = @{$_[1]};
  my $prec = $_[2];
  my $eps = 1e-8;

  my $n = $#data1;

  if ($n != $#data2) { return 1 };
  my $maxdiff = 0;
  for (my $i = 0; $i <= $n; $i++) {
    my $diff = 0;
    my $comp = 0;
    #if (abs($data1[$i]) != 0)
    #{
    #  $diff = abs($data1[$i] - $data2[$i])/abs($data1[$i]);
    #  $comp = $prec;
    #} 
    #else 
    {
      $diff = abs($data1[$i] - $data2[$i]);
      $comp = ($prec * abs($data1[$i]) + $eps);
    }
    if ($diff > $comp) {
      if ($diff > $maxdiff) {
         $maxdiff = $diff;
      }
    }
  }

  return $maxdiff;
}


sub compare_string_vectors($$) {

  my @data1 = @{$_[0]};
  my @data2 = @{$_[1]};

  my $n = $#data1;

  if ($n != $#data2) { return 1 };
  my $difference = 0;
  for (my $i = 0; $i <= $n; $i++) {
    if (trim($data1[$i]) ne trim($data2[$i])) {
      $difference = 1;
      last;
    }
  }

  return $difference;
}




sub read_checks_from_config($$) {

  my $mine = $_[0];
  my $opts = $_[1];

  open(CF, $config) || return 1;
  
  while (<CF>) {
     
     /\s*#/ && next;
     
     chomp $_;

     # format is:
     #  check = filename : var1 var2 var3
     if (/\s*check\s*=([^=:]+):?([^=:]*)/) {
       my @words = parse_line('\s+', 1, trim($2));
       $mine->{trim($1)} = [@words];
     }
     elsif (/\s*inputfile\s*=([^=]+)/) {
       $opts->{"inputfile"} = trim($1);
     }
     elsif (/\s*outputdir\s*=([^=]+)/) {
       $opts->{"outputdir"} = trim($1);
     }
     elsif (/\s*referencedir\s*=([^=]+)/) {
       $opts->{"referencedir"} = trim($1);
     }
     elsif (/\s*commandline\s*=([^=]+)/) {
       $opts->{"commandline"} = trim($1);
     }
  }
  close CF;

  return 0;
}



sub trim($) {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}



sub get_option($) {
  return $options{$_[0]} if defined($options{$_[0]});

  return "";
}




sub set_option($$) {
  $options{$_[0]} = $_[1];
}




sub send_email($)
{
  my $rcpt = get_option("mail_rcpt");
  my $from = get_option("mail_sender");
  my $port = get_option("smtp_port");
  my $user = get_option("smtp_user");
  my $passwd = get_option("smtp_passwd");
  my $subject = "";
  my $msg = "";

  $subject = $_[0]->{'subject'} if defined($_[0]->{'subject'});
  $msg = $_[0]->{'message'} if defined($_[0]->{'message'});

  if (length($rcpt) > 0) {
    print("Sending report to $rcpt\n") if verbose();

    my $smtp;
    if (length($passwd) > 0) { 
      $smtp = Net::SMTP::SSL->new(get_option("smtp_server"),
                              Port => $port,
                              Debug => 0);

      $smtp->auth($user, $passwd) or die "Could not authenticate $!";
    }
    else {
      $smtp = Net::SMTP->new(get_option("smtp_server"),
                             Debug => 0);
    }
    $smtp->mail($from);
    my @addr = split(',', $rcpt);
    $smtp->to(@addr);

    $smtp->data();
    $smtp->datasend("To: $rcpt\n");
    $smtp->datasend("Subject: $subject\n");
    $smtp->datasend($msg);
    $smtp->dataend();

    $smtp->quit;
  }
}




sub print_result($) {

  if ($_[0] == 0) {
    print(color("green"), "ok\n", color("reset"));
  } else {
    print(color("red"),  "failed: max|x - x_ref| = $_[0]\n", color("reset"));
  }
}




sub sprint_result($) {

  if ($_[0] == 0) {
    return (color("green"), "ok\n", color("reset"));
  } else {
    return(color("red"),  "failed\n", color("reset"));
  }
}




sub pad_to_textwidth($$) {
  my $char = $_[1];
  $char = " " if (length($char) != 1); 
  my $string = $_[0];
  my $pad = $textwidth - length($_[0]);
  for (my $i = 0; $i < $pad; $i++) {
    $string = "$string$char";
  }
  return $string;
}
  

sub split_variable($) {

  my $var = "";
  my $prec = $default_precision;
  my $type = 'n'; # for numeric data

  $_[0] =~ /([^\s()<>]+)(<[sn]>)?(\(([\s\deE\+-\.]+)\))?/;
  $var = $1 if defined($1);
  $prec = 1.0 * $4 if defined($4);
  $type = 's' if (defined($2) && $2 eq '<s>');

  return ($var, $prec, $type);
}


sub verbose {
  return $options{"verbose"};
}




END { }

1;
