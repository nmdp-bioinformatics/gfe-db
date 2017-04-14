#!/bin/sh

working=`pwd`
git clone https://github.com/nmdp-bioinformatics/service-gfe-submission
cd service-gfe-submission
tar -xzf hap1.2.tar.gz
export PATH=${working}/hap1.2:${working}/:$PATH
wget http://www.clustal.org/omega/clustalo-1.2.3-Ubuntu-x86_64
mv clustalo-1.2.3-Ubuntu-x86_64 hap1.2/clustalo
chmod a+x hap1.2/clustalo
dpkg --install /opt/ngs-tools_1.9.deb

export PATH=${working}/hap1.2:${working}/ngs-tools/bin:${working}/service-gfe-submission/gfe_submission/bin:$PATH
git clone https://github.com/nmdp-bioinformatics/service-gfe-submission
cd service-gfe-submission && mv -i GFE_Submission gfe_submission && cd gfe_submission
curl -LO http://xrl.us/cpanm 
perl cpanm --quiet --notest Test::More@1.302075
perl cpanm --quiet --notest YAML@1.21
perl cpanm --quiet --notest Moose@2.2004
perl cpanm --quiet --notest Dancer@1.3202
perl cpanm --quiet --notest JSON@2.90
perl cpanm --quiet --notest Dancer::Plugin::Swagger@0.2.0
perl cpanm --quiet --notest Log::Log4perl@1.48
perl cpanm --quiet --notest JSON::Schema::AsType@0.1.0
perl cpanm --quiet --notest Plack::Handler::Starman
perl cpanm --quiet --notest Plack@1.0042
perl cpanm --quiet --notest Template@2.26
perl cpanm --quiet --notest Getopt::Long
perl cpanm --quiet --notest LWP::UserAgent@6.17
perl cpanm --quiet --notest Dancer::Plugin::Swagger@0.2.0
perl cpanm --quiet --notest Net::SSLeay@1.58
perl cpanm --quiet --notest JSON::Schema::AsType@0.1.0
perl cpanm --quiet --notest IO::Socket::SSL@2.044
perl cpanm --quiet --notest REST::Client@273
perl cpanm --quiet --notest Math::Round@0.07
perl cpanm --quiet --notest File::Spec@3.40
perl cpanm --quiet --notest XML::DOM@1.46
perl cpanm --quiet --notest Try::Tiny@0.28
perl cpanm --quiet --notest --skip-satisfied --installdeps .
perl Makefile.PL 
make install

export PWD=${working}
cat neo4j.txt | perl -ne 'chomp;my $dir = $ENV{PWD};$_ =~ s/\%DATADIR\%/$dir/;print $_,"\n";' > neo4j.conf


