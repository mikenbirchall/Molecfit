# Determine the differnce in source files from the molecfit kit and the specified updates

$kit_dirname="molecfit-kit-4.2.3";
$mod_dirname="updates";

#print "TEST DIFF\n";
#print "PERLV = $ENV{'PERLV'} \n";

sub ListMods  {
    $dir=$_[0];
    @flst=glob("$dir/*");
    @rlst=();
    foreach my $file (@flst) {
        if (-d $file) {
            @split_lst=split("/",$file);
            $tail=@split_lst[-1];
            push @rlst, $tail}
    }
    return @rlst;
};# end ListMods

# Parse the Arguments for the root director, virgin source code and updates directory
$nargs=$#ARGV+1;
print "No of arguments = $nargs\n";
if ($nargs==0) {
    print "Must specify root directory\n";
    exit;
}
print "ARG1 = $ARGV[0]\n";
$rdir=$ARGV[0];

if ( ! -e $rdir) {
    print "Error cannot find specified root directory $rdir\n";
    exit;
}

if ( ! -d $rdir) {
    print "Error specified root path is not a directory $rdir\n";
    exit;
}

# Define the kit directory and the update directory
$kitdir=join("/",$rdir,$kit_dirname);
$moddir=join("/",$rdir,$mod_dirname);
if ( ! -d $kitdir) {
    print "Error cannot find the molecfit kit directory $kitdir\n";
    exit;
}
if ( ! -d $moddir) {
    print "Error cannot find the modifications directory $moddir\n";
    exit;
}

if ($nargs<2) {
    # get the list of sub directories in the mod directores
    @mod_list = glob("$moddir/*");
    $nmods=@mod_list;
    if ($nmods==0) {
        print "Error no update modules found in $moddir\n";
        exit;
    }
    # Pick the first one for use
    $use_mod=@mod_list[0];

} else {
    $use_mod=join("/",$moddir,$ARGV[1]);
    if (! -e $use_mod) {
        print "Error cannot find update module $use_mod\n";
        exit;
    }

}


print "Kit Directory=$kitdir\n";
print "Modules Directory=$moddir\n";
print "Default Mod=$use_mod\n";

@mlst=ListMods($moddir);
 print " MOD LST=@mlst\n";
