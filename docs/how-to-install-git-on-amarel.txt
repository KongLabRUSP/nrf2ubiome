# how to install git in your own folder so it can be invoked by Rstudio on Amarel.
# R Wu. Dec 2018

# 0. create folders for git
cd ~
mkdir -p github tools/git && cd github

# 1. get git source (use either one of the following two methods)
# method 1,  login to amarel.rutgers.edu or ondemand.hpc.rutgers.edu via ssh #or  via ondemand web terminal
# (an old version) of git is already there. Run
git clone git://git.kernel.org/pub/scm/git/git.git

# method 2, download source from https://www.kernel.org/pub/software/scm/git, https://mirrors.edge.kernel.org/pub/software/scm/git/git-2.20.0.tar.gz or  https://github.com/git/git/releases. https://github.com/git/git/archive/v2.20.0.tar.gz

$ wget https://github.com/git/git/archive/v2.20.0.tar.gz # Version number may vary
$ tar -zxf git-2.0.0.tar.gz # Version number may vary


# 2 compile git
$ cd git-2.0.0 # Version number may vary
$ make configure
$ ./configure --prefix=${HOME}/tools/git
$ make all
$ make install


 # Note: 
 # The /projects folder is not visible in the node that hosts the Rstudio server. so git must be installed in your home folder rather than the (shared) project folder.

# 3. Set git in  Rstudio: 
# a, go to Tools  -> Global Options -> GIT/SVN, Enter the correct location of git excutable, e.g., /home/netid/tools/git/bin/git
# b, restart Rstudio and enjoy. (Note git repo should be cloned via http(s) as ssh is not supported somehow)
 
 
 # set new version of R 
 wget R-x.x.x.tar.gz
 tar -xzf R-x.x.x.tar.gz
 cd to the fold
 ./configure --prefix=/home/rw409/tools/R/3.5.1_node --enable-utf --enable-unicode-properties --enable-jit --disable-cpp --with-x=no
 make
 make install
 
 # Note: The binary R can only be run on the same node where the binary was compiled. i.e., you need to comile R on the compute node if you wish to run it there. And the same for login nodes (amarel1/2) and ondemand node.
 
  
 # then export so rstudio can use the new installation. DID NOT WORK with Rstudio server!!!
 export RSTUDIO_WHICH_R=/home/rw409/tools/R/3.5.1_node/bin/R
 export PATH=/home/rw409/tools/R/3.5.1/bin/:$PATH
 
 # https://cran.r-project.org/doc/manuals/R-admin.html#Configuration-options
# https://cran.r-project.org/doc/manuals/R-admin.html#Simple-compilation
 # Section A.1