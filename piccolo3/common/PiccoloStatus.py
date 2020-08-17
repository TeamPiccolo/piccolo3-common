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

__all__ = ['PiccoloSchedulerStatus','PiccoloSpectrometerStatus']

import enum

class PiccoloSchedulerStatus(enum.Enum):
    active = 1
    done = 2
    suspended = 3
    deleted = 4
    
class PiccoloSpectrometerStatus(enum.IntEnum):
    NO_WORKER = enum.auto()
    POWERED_OFF = enum.auto()
    DROPPED =  enum.auto()
    CONNECTING = enum.auto()
    DISCONNECTED = enum.auto()
    IDLE = enum.auto()
    RECORDING = enum.auto()
    AUTOINTEGRATING = enum.auto()
