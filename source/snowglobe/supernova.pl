#!/usr/bin/perl
# arguments:  channel file name, experiment configuration name, noweight
# e.g.  ./supernova.pl water wc100kt30prct

# flux will fixed at ./fluxes/pinched_0.dat
$fluxname = "pinched_0";
$channame = $ARGV[0];
$expt_config = $ARGV[1];


$chanfilename = "channels/channels_".$channame.".dat";

# If we are using a non-standard binning, this binning must be defined at the top of
# the channel file.  The format is:
# % bins emin emax sampling_points sampling_min sampling_max
# and all energy values must be in GeV
# the default values are...
$ewins{'bins'} = "200";
$ewins{'emin'} = "0.0005";   #GeV
$ewins{'emax'} = "0.100";    #GeV
$ewins{'sampling_points'} = "200";
$ewins{'sampling_min'} = "0.0005"; #GeV
$ewins{'sampling_max'} = "0.100"; #GeV

$energywindowset = "standard"; # or _he for high energy Ev

open(CHANFILE,$chanfilename);# open $chanfilename as filehandler CHANFILE
$firstchanline = <CHANFILE>; # read a line of CHANFILE
$firstchanline =~ s/^s+//;   # remove leading whitespace
if (index($firstchanline, "%") != -1){ # if % not found, return -1
    @binningarray = split /\s+/, $firstchanline; # array(@,start with 0), \s+ as delimiter  
    $ewins{'bins'} = @binningarray[1];
    $ewins{'emin'} = @binningarray[2];
    $ewins{'emax'} = @binningarray[3];
    $ewins{'sampling_points'} = @binningarray[4];
    $ewins{'sampling_min'} = @binningarray[5];
    $ewins{'sampling_max'} = @binningarray[6];
}
close(CHANFILE);
# these energy window/binning values are now ready for the preamble

# now that we've read in the channel file, we should know what energy windows we have
if ($ewins{'emax'} > 0.199 and $ewins{'emax'} < 0.201) {
    $energywindowset = "_he";
    print "using high energy window!\n";
}

# Create the globes file

$globesfilename = "supernova.glb";

open(GLOBESFILE,">$globesfilename"); # open for writing


# Here we add the globes preamble, with any modifications needed for rebinning taken care of
open(PREAMBLE,"glb/preamble.glb");
$ready_to_modify = 0;
while(<PREAMBLE>) {
    if (index($_, "Energy window") != -1){$ready_to_modify = 1;} # $_ is the line feed in
    $ourline = $_;
    if($ready_to_modify) {
        keys %ewins;  # reset the internal iterator so a prior each() doesn't affect the loop
        while(my($k, $v) = each %ewins){ # my define two local enclosing variables
                                         # each returns (key,value) pair of hash ewins iteratively
           if (index($ourline, $k) != -1){
               @vallist = split /\s+/, $ourline;
               $ourline =~ s/$vallist[2]/$v/;
               last;   # we have done the replacement, no need to continue
                       # last <=> break in C
           }
       }   # end while loop over ewin replacements
    }  # any necessary modification is done
    print GLOBESFILE $ourline; # print $ourline into GLOBESFILE
}
close(PREAMBLE);

# Put in the flux info

$fluxfilename = "fluxes/".$fluxname.".dat";
unless(-f $fluxfilename) {
    print "Flux file name ",$fluxfilename," not found\n";
    exit;
}


open(FLUX,"glb/flux.glb");

while(<FLUX>) {
# Replace the flux file name with the input argument
    if (/flux_file/) {
	$_ = "        \@flux_file=  \"snowglobe/".$fluxfilename."\"\n";
    }
    print GLOBESFILE $_;
}
close(FLUX);

# Now go through the channels and put in the relevant lines

# First, smearing for each channel

unless (-f $chanfilename) {
    print "Channel file name ",$chanfilename," not found\n";
    exit;
}
open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# skip energy window info line
    if (/%/) {next;}

# Grab the channel name
    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];

    $output_line = "include \"snowglobe/smear/smear_".$chan_name."_".$expt_config.".dat\"\n";


    print GLOBESFILE $output_line;
}

close(CHANFILE); # end of section ``Flux data`` in creating supernova.glb
                 # if no background, otherwise end will hit at line 202

#  Detector info

$detfilename = "detector_configurations.dat";

unless (-f $detfilename) {
    print "Detector file name ",$detfilename," not found\n";
    exit;
}

open(DETFILENAME,$detfilename);

while(<DETFILENAME>) {
    chop($_); # remove the last character from $_
    
# Skip comments
    if (/\#/) {next;}
    $_=~s/\s*//; # remove leading blanks, leading non-blank will be matched beacuse of the *
    $_=~s/\s+/ /g; # make all(g for all, no g for first match) \s to be one char empty

    @stuff = split(/\ /,$_);

    $detname = $stuff[0];
    $masses{$detname} = $stuff[1];
    $normfact{$detname} = $stuff[2];

    if ($detname eq "" || $masses{$detname} eq "" || $normfact{$detname} eq ""){ next;}
}
close(DETFILENAME);

# Mass and target normalization by species

#%masses = ("sk1",22.5,"sk2",22.5,"100kt15pc",100, "100kt30pc",100,"ar17",17,"scint50kt",50,"halo",1.0);
#%normfact = ("sk1",2/18,"sk2",2/18,"100kt15pc",2/18, "100kt30pc",2/18,"ar17kt",1/40,"scint50kt",2/14,"halo",1/208);

if ($masses{$expt_config} == 0){
    print "Error: please enter a valid experiment configuration\n";
    exit;

}

# reference target mass
$target_mass = sprintf("%15.8f",$masses{$expt_config}*$normfact{$expt_config});  # This is ktons of free protons

print "Experiment config: ",$expt_config,"  Mass: ",$masses{$expt_config}," kton \n";
#print " Target mass: ",$normfact{$expt_config}," ",$target_mass,"\n";

# Add the background smearing here, for the given detector configuration
#  (Not yet implemented: for multiple background channels, 
# read them from a file labeled by detector configuration)
$do_bg = 0;
$bg_chan_name = "bg_chan";

$bg_filename = "backgrounds/".$bg_chan_name."_".$expt_config.".dat";
if (-e $bg_filename) {
    $do_bg = 1;
    print "Using background file ",$bg_filename,"\n";
} else {
    print "No background file for this configuration\n";
}

if ($do_bg == 1) {
    $output_line = "include \"snowglobe/smear/smear_".$bg_chan_name."_".$expt_config.".dat\"\n";
    print GLOBESFILE $output_line; 
}

open(DETECTOR,"glb/detector.glb");
while(<DETECTOR>) {
# Replace the flux file name with the input argument
    if (/mass/) {
	$_ = "\$target_mass=  ".$target_mass."\n";
    }
    print GLOBESFILE $_;
}
close(DETECTOR);


print GLOBESFILE "\n /******** Cross-sections *********/\n \n";

# Now the cross-sections.  Note that some of these are repeated 
# even thougn it is not necessary (xscns for several flavors can be in the 
# same file).
#  This is just to make a consistent loop over channels.

open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# skip energy window info line
    if (/%/) {next;}

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];

    $output_line = "cross(\#".$chan_name.")<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@cross_file= \"snowglobe/xscns/xs_".$chan_name.".dat\"\n";
#    print $output_line;
    print GLOBESFILE $output_line;

    $output_line = "\>\n";
#    print $output_line;
    print GLOBESFILE $output_line;

}

close(CHANFILE);

# Now the fake bg channel cross section, if it exists

if ($do_bg == 1) {
    $output_line = "cross(\#".$bg_chan_name.")<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@cross_file= \"snowglobe/xscns/xs_zero.dat\"\n";
#    print $output_line;
    print GLOBESFILE $output_line;
    
    $output_line = "\>\n";
#    print $output_line;
    print GLOBESFILE $output_line;
}

print GLOBESFILE "\n \/******** Channels *********\/\n \n";

# Now, the channel definitions

unless(-f $chanfilename) {
    print "Channel file name ",$chanfilename," not found\n";
    exit;
}

open(CHANFILE,$chanfilename);

while(<CHANFILE>) {

# skip energy window info line
    if (/%/) {next;}

# Grab the channel name

    $_=~s/\s*//;
    $_=~s/\s+/ /g;

    @stuff = split(/\ /,$_);

    $chan_name = $stuff[0];
    $cpstate = $stuff[2];
    $inflav = $stuff[3];

    $output_line = "channel(\#".$chan_name."_signal)<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@channel= \#supernova_flux:  ".$cpstate.":    ".$inflav.":     ".$inflav.":    \#".$chan_name.":    \#".$chan_name."_smear\n"; 
#    print $output_line;
    print GLOBESFILE $output_line;

# Get the post-smearing efficiencies by channel

    $eff_file = "effic/effic_".$chan_name."_".$expt_config.".dat";
#    print $eff_file,"\n";
    open(EFF_FILE,$eff_file);
    while(<EFF_FILE>) {
	$output_line = "       \@post_smearing_efficiencies = ".$_;
	print GLOBESFILE $output_line;
    }
    close(EFF_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

    $output_line = "\n\>\n\n";
    print GLOBESFILE $output_line;

}

close(CHANFILE);

if ($do_bg == 1) {
# Now make a fake channel for the background (just one possible bg file for now)

# This is dummy info
    $cpstate = "-"; $inflav = "e";

    $output_line = "channel(\#".$bg_chan_name."_signal)<\n";
    print GLOBESFILE $output_line;

    $output_line = "      \@channel= \#supernova_flux:  ".$cpstate.":    ".$inflav.":     ".$inflav.":    \#".$bg_chan_name.":    \#".$bg_chan_name."_smear\n"; 
#    print $output_line;
    print GLOBESFILE $output_line;

# Get the pre-smearing backgrounds by channel

    $bg_file = "backgrounds/".$bg_chan_name."_".$expt_config.".dat";
    print $bg_file,"\n";
    open(BG_FILE,$bg_file);
    while(<BG_FILE>) {
	$output_line = "       \@pre_smearing_background = ".$_;
	print GLOBESFILE $output_line;
    }
    close (BG_FILE);

# Try crazy reformatting
# or else get mysterious (but apparently harmless?) error from GLoBES
# Should try to track it down...  -> error goes away with development GLoBES version

$output_line = "\n\>\n\n";
print GLOBESFILE $output_line;
}

# End-matter
if ($energywindowset eq "standard"){
    open(POSTAMBLE,"glb/postamble.glb"); # up to 100 MeV
}
elsif ($energywindowset eq "_he"){
    open(POSTAMBLE,"glb/postamble_he.glb"); # up to 200 MeV
}

while(<POSTAMBLE>) {
    print GLOBESFILE $_;
}


close(POSTAMBLE);

close(GLOBESFILE);  # end of creating supernova.glb
