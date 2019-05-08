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

"""
.. moduleauthor:: Magnus Hagdorn <magnus.hagdorn@ed.ac.uk>
.. moduleauthor:: Iain Robinson <iain.robinson@ed.ac.uk>
"""

__all__ = ['PiccoloSpectraList','PiccoloSpectrum']

from collections import MutableMapping, MutableSequence
from datetime import datetime
import json
import os.path
import numpy

protectedKeys = ['Direction','Dark','Datetime']

class PiccoloSpectraList(MutableSequence):
    """a collection of spectra

    The object behaves like a python list. The object supports chunking which
    can be used to transfer data across a slow network.
    """

    _NCHUNKS = 300

    def __init__(self,run=None,batch=0,seqNr=0,data=None):
        """
           :param run:   the name of the run
           :param batch: the batch number of the spectra collection;
                         used for constructing the output file name
           :param seqNr: the sequence number of the spectra collection;
                         used for constructing the output file name
           :param data:  string containing JSON serialised version of the
                         spectra list. This can be used to create a
                         SpectraList from JSON
        """
        self._spectra = []
        self._run = run
        self._batch = batch
        self._seqNr = seqNr
        self._chunkID = None

        # initialise from json if available
        if data is not None:
            self._initFromData(data)

    def __getitem__(self,i):
        return self._spectra[i]
    def __setitem__(self,i,y):
        assert isinstance(y,PiccoloSpectrum)
        self._spectra[i] = y
    def __delitem__(self,i):
        raise RuntimeError('cannot delete spectra')
    def __len__(self):
        return len(self._spectra)
    def insert(self,i,y):
        """insert a new Spectrum
        :param i: position in the list after which the new spectrum should be inserted
        :param y: the spectrum object to be inserted
        :type y: PiccoloSpectrum"""
        assert isinstance(y,PiccoloSpectrum)
        self._spectra.insert(i,y)

    def _initFromData(self,data):
        self._spectra = []
        if isinstance(data,str):
            data = json.loads(data)
        for s in data['Spectra']:
            self._seqNr = s['Metadata']['SequenceNumber']
            self.append(PiccoloSpectrum(data=s))

    @property
    def NCHUNKS(self):
        """the number of total chunks"""
        return self._NCHUNKS

    @property
    def chunk(self):
        """the current chunk"""
        return self._chunkID

    @property
    def complete(self):
        """return True if all spectra are complete"""
        if self._chunkID is None:
            return True
        for s in self:
            if not s.complete:
                return False
        return True

    @property
    def run(self):
        """the run name, used in generating output filename"""
        return self._run
    
    @property
    def batch(self):
        """the batch number, used in generating output filename"""
        return self._batch
    
    @property
    def seqNr(self):
        """the sequency number, used in generating output filename"""
        return self._seqNr

    @property
    def outName(self):
        """the output name, generated from the batch and sequence number"""
        n = 'b{0:06d}_s{1:06d}.pico'.format(self.batch,self.seqNr)
        if self.run is not None:
            n = os.path.join(self.run,n)
        return n

    @property
    def directions(self):
        """a set containing all directions present in the spectra list"""
        dirs = set()
        for s in self._spectra:
            dirs.add(s['Direction'])
        return list(dirs)

    @property
    def haveDark(self):
        for s in self._spectra:
            if s['Dark']:
                return True
        return False

    @property
    def haveLight(self):
        for s in self._spectra:
            if not s['Dark']:
                return True
        return False

    def haveSpectrum(self,spectrum):
        """check if a particular spectrum type is available
        :param spectrum: must be either Light or Dark"""

        if spectrum=='Dark':
            return self.haveDark
        elif spectrum=='Light':
            return self.haveLight
        else:
            raise KeyError('spectrum must be one of Dark or Light')

    def getSpectra(self,direction,spectrum):
        """extract a particular spectrum by type
        :param direction: the direction
        :param spectrum: must be either Light or Dark"""
        if spectrum == 'Dark':
            dark = True
        elif spectrum == 'Light':
            dark = False
        else:
            raise KeyError('spectrum must be one of Dark or Light')
        spectra = []
        for s in self._spectra:
            if s['Direction'] == direction and s['Dark'] == dark:
                 spectra.append(s)
        return spectra

    def serialize(self,pretty=True,pixelType='list',spectrum=None):
        """serialize to JSON

        :param pretty: when set True (default) produce indented JSON
        :param pixelType: set the pixel type
        :param spectrum: select spectrum type (Dark or Light) or both when None"""
        if spectrum == 'Dark':
            dark = True
        elif spectrum == 'Light':
            dark = False
        elif spectrum is None:
            dark = None
        else:
            raise KeyError('spectrum must be one of Dark or Light or None')

        spectra = []
        for s in self._spectra:
            # add sequence number to meta data
            s['SequenceNumber'] = self._seqNr
            if dark is None or s['Dark'] == dark:
                spectra.append(s.as_dict(pixelType))
        root = {'Spectra':spectra}

        if pretty:
            return json.dumps(root, sort_keys=True, indent=1)
        else:
            return json.dumps(root)

    def write(self,prefix='',clobber=True, split=True):
        """write spectra to file

        :param prefix: output prefix
        :param clobber: boolean whether files should be overwritten or not
        :param split: when set to True split files into light and dark spectra"""

        outName = os.path.join(prefix,self.outName)
        outDir = os.path.dirname(outName)

        if not os.path.exists(outDir):
            os.makedirs(outDir)

        if split:
            for s in ['Dark','Light']:
                if self.haveSpectrum(s):
                    o = "{}_{}".format(outName, s.lower())
                    # The above producues a filename with the extension (.pico) in the middle of the filename:
                    #      b000000_s000000.pico_dark
                    # The two lines below fix it.
                    o = o.replace('.pico', '')
                    o = o + '.pico'
                    if not clobber and os.path.exists(o):
                        raise RuntimeError('{} already exists'.format(o))
                    with open(o,'w') as outf:
                        outf.write(self.serialize(spectrum=s))
        else:
            if not clobber and os.path.exists(outName):
                raise RuntimeError('{} already exists'.format(outName))
            with open(outName,'w') as outf:
                outf.write(self.serialize())

    def getChunk(self,idx):
        """get a particular chunk
        :param idx: the chunk index
        :return: JSON string containing the data"""
        assert isinstance(idx,int)
        assert idx>=0 and idx < self.NCHUNKS

        if idx == 0:
            # first chunk is special, copy all the meta data
            data = self.serialize(pretty=False,pixelType='size')
        else:
            data = []
            for s in self._spectra:
                data.append(s.getChunk(idx-1,self.NCHUNKS-1).tolist())
            data = json.dumps(data)

        return data

    def setChunk(self,idx,data):
        """add a particular chunk
        :param idx: the chunk index
        :param data: the chunk to be added
        :type data: JSON string"""
        assert isinstance(idx,int)
        assert idx>=0 and idx < self.NCHUNKS

        if idx == 0:
            self._initFromData(data)
        else:
            assert self._chunkID is not None
            data = json.loads(data)
            assert len(data) == len(self._spectra)
            for i in range(len(data)):
                self._spectra[i].setChunk(idx-1,self.NCHUNKS-1,data[i])
        self._chunkID = idx


class PiccoloSpectrum(MutableMapping):
    """An object containing an optical spectrum.

    The PiccoloSpectrum object behaves like a python dictionary. It can be
    initialised from a JSON string. It also supports chunking."""
    def __init__(self,data=None):
        """:param data: JSON string used to initialise object"""
        self._meta = {}
        self._meta['Direction'] = 'Missing metadata'
        self._meta['Dark'] = 'Missing metadata'
        self._meta['Type'] = 'Missing metadata'
        self._pixels = None

        self._corrected = None
        
        self._complete = False
        self.setDatetime()

        # initialise from json if available
        if data is not None:
            if isinstance(data,str):
                data = json.loads(data)
            for key in data['Metadata']:
                self._meta[key] = data['Metadata'][key]
            if isinstance(data['Pixels'],int):
                # create a list of correct size
                self._pixels = -numpy.ones(data['Pixels'],dtype=numpy.int)
            else:
                self.pixels = data['Pixels']

    @property
    def complete(self):
        """whether all chunks have been set"""
        return self._complete

    def __getitem__(self,key):
        return self._meta[key]

    def __setitem__(self,key,value):
        if key in protectedKeys:
            raise KeyError('field {0} is a protected key'.format(key))
        self._meta[key] = value

    def __delitem__(self,key):
        if key in protectedKeys:
            raise KeyError('field {0} is a protected key'.format(key))
        del self._meta[key]

    def __iter__(self):
        return iter(self._meta)

    def __len__(self):
        return len(self._meta)

    def setDirection(self,direction):
        self._meta['Direction'] = direction
    
    def setUpwelling(self,value=None):
        """set the direction to upwelling
        :param value: set direction to downwelling when False
        :type value: bool"""
        if value is None:
            self._meta['Direction'] = 'Upwelling'
        else:
            assert isinstance(value,bool)
            if value:
                self.setUpwelling()
            else:
                self.setDownwelling()

    def setDownwelling(self):
        """set direction to downwelling"""
        self._meta['Direction'] = 'Downwelling'

    def setDark(self,value=None):
        """set spectrum to dark
        :param value: set to light when False
        :type value: bool"""
        if value is None:
            self._meta['Dark'] = True
            self._meta['Type'] = 'dark'
        else:
            assert isinstance(value,bool)
            self._meta['Dark'] = value
            if value:
                self._meta['Type'] = 'dark'
            else:
                self._meta['Type'] = 'light'

    def setLight(self):
        """set spectrum to light"""
        self._meta['Dark'] = False
        self._meta['Type'] = 'light'

    def setDatetime(self,dt=None):
        """set date and time when spectrum is recorded
        :param dt: date and time of recording, use now if set to None
        :type dt: datetime"""
        if dt is  None:
            ts = datetime.now()
        elif isinstance(dt,datetime):
            ts = dt
        else:
            ts = datetime.strptime(dt,'%Y-%m-%dT%H:%M:%S')

        self._meta['Datetime'] = '{}Z'.format(ts.isoformat())

    @property
    def pixels(self):
        """the pixels"""
        if self._pixels is None:
            raise RuntimeError('The pixel values have not been set.')
        if len(self._pixels) == 0:
            raise RuntimeError('There are no pixels in the spectrum.')
        return self._pixels
    @pixels.setter
    def pixels(self,values):
        tmp = numpy.minimum(values,200000)
        self._pixels = numpy.array(tmp,dtype=numpy.int)
        self._corrected = None

    @property
    def corrected_pixels(self):
        #TODO: apply dark pixel
        if self._corrected is None:
            dark = 0.
            cpoly = numpy.poly1d(numpy.array(self['NonlinearityCorrectionCoefficients'])[::-1])
            self._corrected = dark + (self.pixels-dark)/cpoly(self.pixels-dark)
        return self._corrected
    
    def getNumberOfPixels(self):
        """the number of pixels"""
        return len(self.pixels)

    def computeWavelength(self,i):
        """compute the ith wavelength
        :param i: the wavelength index"""
        C = self['WavelengthCalibrationCoefficients']
        w = 0
        for j in range(len(C)):
            w = w+ C[j]*i**j
        return w

    @property
    def waveLengths(self):
        """the list of wavelengths"""
        w = []
        if 'WavelengthCalibrationCoefficients' in self.keys():
            for i in range(self.getNumberOfPixels()):
                w.append(self.computeWavelength(i))
        else:
            w = list(range(self.getNumberOfPixels()))
        return w

    def as_dict(self,pixelType='array'):
        """represent spectrum as a dictionary
        :param pixelType: how the pxiels are represented

        .. note::
          the pixelType should be one of
           * 'array' (default): the pixel values are stored as a numpy.array
           * 'list': the pxiel values are stored as a list; useful for serialisation
           * 'size': store the size of pixel array, useful for initialising chunked array"""
        spectrum = {}
        spectrum['Metadata'] = dict(self.items())
        if pixelType == 'size':
            spectrum['Pixels'] = self.getNumberOfPixels()
        elif pixelType == 'list':
            spectrum['Pixels'] = self.pixels.tolist()
        elif pixelType == 'array':
            spectrum['Pixels'] = self.pixels
        else:
            raise RuntimeError('unknown pixel type %s'%pixelType)
        return spectrum

    def serialize(self,pretty=True):
        """serialize to JSON string
        :param pretty: pretty print JSON"""
        spectrum = self.as_dict(pixelType='list')
        if pretty:
            return json.dumps(spectrum, sort_keys=True, indent=1)
        else:
            return json.dumps(spectrum)

    def getChunk(self,idx,nChunks):
        """get a chunk
        :param idx: the chunk index
        :param nChunks: the total number of chunks"""
        return self.pixels[range(idx,self.getNumberOfPixels(),nChunks)]

    def setChunk(self,idx,nChunks,data):
        """set a chunk
        :param idx: the chunk index
        :param nChunks: the total number of chunks
        :param data: the chunk data
        :type data: JSON string"""
        if idx!=nChunks-1:
            self._complete = False
        else:
            self._complete = True
        rng = range(idx,self.getNumberOfPixels(),nChunks)
        assert len(rng) == len(data)
        self._pixels[rng] = data

if __name__ == '__main__':
    import sys

    if len(sys.argv)>1:
        data = open(sys.argv[1],'r').read()

        spectra = PiccoloSpectraList(data=data)
        print (spectra.directions)
