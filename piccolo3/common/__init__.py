# Copyright 2014-2016 The Piccolo Team
#
# This file is part of piccolo3-common.
#
# piccolo3-common is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# piccolo3-common is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with piccolo3-common.  If not, see <http://www.gnu.org/licenses/>.

from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution('piccolo3-common').version
except DistributionNotFound:
    # package is not installed
    pass

from .piccoloLogging import *
from .PiccoloSpectra import *
from .PiccoloStatus import *
