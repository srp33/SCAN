#Read the Affymetrix's Version 4 CEL file.
#See http://www.stat.lsa.umich.edu/~kshedden/Courses/Stat545/Notes/AffxFileFormats/cel.html for details of the format.

import sys,os,struct,time

class MyCEL(object):
    def _read_integer(self,file):
        return  struct.unpack('i', file.read(4))[0]

    def _read_char(self,file,length):
        return  file.read(length)

    def _read_dword(self,file):
        return  struct.unpack('I', file.read(4))[0]

    def _read_float(self,file):
        return  struct.unpack('f', file.read(4))[0]

    def _read_short(self,file):
        return  struct.unpack('h', file.read(2))[0]

    def read_platform(self, fcel):
        f = open(fcel,"rb") #read binary

        # Reading through values to get to the header
        for i in range(5):
            self._read_integer(f)
        header = self._read_char(f, self._read_integer(f))

        f.close()

        headerItems = header.split("\n")
        datHeader = [x for x in headerItems if x.startswith("DatHeader")][0]
        platform = datHeader.split(" \x14 ")[2].replace(".1sq", "")

        # This means the platform is specified in an alternative way
        if len(platform) == 1:
            platform = datHeader.replace("DatHeader=", "")[4:].lstrip()
            platform = platform[:platform.find(".1sq")]

        # This is a special case
        if platform == "U133AAofAv2":
            platform = "HT_HG-U133A"

        return platform

    def read_cel(self, fcel, coord):
        """Read the .CEL file.
        @param:
            fcel - filename.cel
            coord - {(row,col):probe_id,...}
        @return
            {probe_id:[raw_intensity,probe_sequence], ... }
        """

        entry = {}
        exp = []

        f = open(fcel,"rb") #read binary
        magic_number = self._read_integer(f) #always 64
        version_number = self._read_integer(f) #always 4
        num_columns = self._read_integer(f)
        num_rows = self._read_integer(f)
        num_cells = self._read_integer(f)#rows*columns
        header_length = self._read_integer(f)#length of the header
        header = self._read_char(f,header_length)
        algorithm_name_length = self._read_integer(f)
        algorithm_name = self._read_char(f,algorithm_name_length)
        algorithm_parameters_length = self._read_integer(f)
        algorithm_parameters = self._read_char(f,algorithm_parameters_length)
        headerItems = header.split("\n")
        datHeader = [x for x in headerItems if x.startswith("DatHeader")][0]
        platform = datHeader.split(" \x14 ")[2].replace(".1sq", "")
        cell_margin = self._read_integer(f)
        num_outlier_cells = self._read_dword(f)
        num_masked_cells = self._read_dword(f)
        num_subgrids = self._read_integer(f)

        count = 0
        matchCount = 0
        for c in xrange(num_columns):
            for r in xrange(num_rows):
                #intensity, standard deviation, pixel count
                count += 1
                cell_entries = (self._read_float(f),self._read_float(f),self._read_short(f))
                probe_id = coord.get((r,c)) #probe_id, probe_seq

                if probe_id:
                    entry[probe_id] = int(cell_entries[0])
                    matchCount += 1

        f.close()

        return entry
