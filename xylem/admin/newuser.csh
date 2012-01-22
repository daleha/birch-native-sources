#!/bin/csh
# setup commands for first terminal session
set DB = ~psgendb
echo Setting up .login file
cat $DB/admin/add_to_login >> $home/.login
echo Setting up .cshrc file
cat $DB/admin/add_to_cshrc >> $home/.cshrc
#if !(-e $home/.xsession) then
#   echo Copying .xsession to home directory.
#   cp $DB/admin/standard.xsession $home/.xsession
#endif
echo Done!
echo Logout and login again so that the changes can take effect.
