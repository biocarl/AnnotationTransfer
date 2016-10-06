#Installing SEED-package which includes Scripts for calling the RAST-Server
#source: http://blog.theseed.org/servers/installation/distribution-of-the-seed-server-packages.html
#doku: http://blog.theseed.org/servers/usage/the-rast-batch-interface.html

path=$HOME/bin
source=$HOME/.bashrc

mkdir $path

#Download RAST
wget --no-check-certificate http://blog.theseed.org/downloads/sas.tgz -O sas.tgz

path_rast=$path/sas
mkdir $path_rast

tar -xvf sas.tgz --strip-components=1 -C $path_rast
rm -rf sas.tgz

#Install Seed-Package incl. myRast
cd $path_rast/modules
./BUILD_MODULES

cat << EOF >> $source
#RAST
PERL5LIB=$PERL5LIB:$path_rast/lib:$path_rast/modules/lib
export PERL5LIB

PATH=$PATH:$path_rast/bin
export PATH
#RAST
EOF

#Load ENV
source source

#Verify
perl -e 'use SeedEnv'

<<EOF
if [ -z "$(perl -e 'use SeedEnv')" ];
die 01 "The Path was not set correctly"


warn () {
    echo "$0:" "$@" >&2
}

die () {
    rc=$1
    shift
    warn "$@"
    exit $rc
}
EOF
