#
# Regular cron jobs for the piccolo3-common package
#
0 4	* * *	root	[ -x /usr/bin/piccolo3-common_maintenance ] && /usr/bin/piccolo3-common_maintenance
