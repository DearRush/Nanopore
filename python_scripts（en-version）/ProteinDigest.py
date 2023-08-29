import re
class ProteinDigest:
    def __init__( self, enzyme='trypsin' ):
        self.references = "Harpaz, Y., Gerstein, M. & Chothia, C. Volume changes on protein folding. Structure 2, 641–649 (1994)."
        self.AvN = 6.0221409e23
        self.Angstrom_to_Millilitres_3 = ( 1e-8 ) ** 3
        self.partial_volume_water = {"G": 71.7, "A": 100.3, "V": 150.6,
                                    "L": 178.7, "I": 175.4, "P": 137.2,
                                    "M": 174.9, "C": 122.1, "F": 202.3,
                                    "Y": 205.3, "W": 239.0, "S": 100.7,
                                    "T": 127.6, "H": 163.9, "N": 128.4,
                                    "D": 113.1, "Q": 156.0, "E": 140.2,
                                    "R": 192.8, "K": 170.3}
        self.partial_volume_interior = {"G": 63.8, "A": 90.1, "V": 139.1,
                                    "L": 164.6, "I": 164.9, "P": 123.1,
                                    "M": 167.7, "C": 103.5, "F": 193.5,
                                    "Y": 197.1, "W": 231.7, "S": 94.2,
                                    "T": 120.0, "H": 159.3, "N": 127.5,
                                    "D": 117.1, "Q": 149.4, "E": 140.8,
                                    "R": 192.8, "K": 170.0}
        self.mass = {"G": 57.0519, "A": 71.0788, "V": 99.1326,
                    "L": 113.1594, "I": 113.1594, "P": 97.1167,
                    "M": 131.1926, "C": 103.1388, "F": 147.1766,
                    "Y": 163.1760, "W": 186.2132, "S": 87.0782,
                    "T": 101.1051, "H": 137.1411, "N": 114.1038,
                    "D": 115.0886, "Q": 128.1307, "E": 129.1155,
                    "R": 156.1875, "K": 128.1741}
        self.gravy = {"G": -0.400, "A": 1.800, "V": 4.200,
                    "L": 3.800, "I": 4.500, "P": -1.600,
                    "M": 1.900, "C": 2.500, "F": 2.800,
                    "Y": -1.300, "W": -0.900, "S": -0.800,
                    "T": -0.700, "H": -3.200, "N": -3.500,
                    "D": -3.500, "Q": -3.500, "E": -3.500,
                    "R": -4.500, "K": -3.900}
        
        if enzyme=='trypsin':
            self.re_enzyme = r".(?:(?<![KR](?!P)).)*"
        elif enzyme=='chymotrypsin':
            self.re_enzyme = r".(?:(?<![FWY](?!P)).)*"
        elif enzyme=='chymotrypsin(unsp)':
            self.re_enzyme = r".(?:(?<![FWYML](?!P)).)*"
        elif enzyme=='trypsin-chymotrypsin':
            self.re_enzyme = r".(?:(?<![KRFWY](?!P)).)*"
    
    def get_peptide( self, S, state='aq', cystein_mod=57 ):
        peptide_bond = ( ( len( S ) - 1 ) * 10.4 )
        if state=='aq':
            return sum([ self.partial_volume_water[s] for s in S ]) - peptide_bond
        elif state=='int':
            return sum([ self.partial_volume_interior[s] for s in S ]) - peptide_bond
        elif state=='mass':
            return 18.01528 + sum([ self.mass[s] for s in S ]) + S.count("C")*cystein_mod
        elif state=='gravy':
            return sum([ self.gravy[s] for s in S ]) / len( S )
        elif state=='aliphatic':
            return (S.count("A")/len(S)) + 2.9*(S.count("V")/len(S)) + 3.9*((S.count("L")/len(S))+(S.count("I")/len(S)))
        elif state=='OD280':
            return get_peptide( S, state='E280' ) / get_peptide( S, state='mass' )
        elif state=='E280':
            E280 = ( S.count("Y") * 1490 ) + ( S.count("W") * 5500 ) + ( S.count("C") * 125 )
            return E280
    
    def find_glycosites( self, sequence ):
        glycosites=[]
        for c, s in enumerate( list( sequence ) ):
            if s=="T" or s=="S":
                try:
                    if list( sequence )[c-2]=="N":
                        glycosites.append(c)
                except:
                    pass
        return glycosites

    
    def digest( self, sequence, min_len=3, glycosylate=False ):
        if glycosylate==True:
            glycosites = self.find_glycosites( sequence )
        fragment_start = []
        fragment_end = []
        fragment_seq = []
        fragment_mass = []
        fragment_volume = []
        fragment_gravy = []
        k = 0
        for c, i in enumerate( re.findall(self.re_enzyme, sequence) ):
            if len( i ) > min_len:
                fragment_start.append( k )
                k += len( i )
                fragment_end.append( fragment_start[-1] + len( i ) )
                fragment_seq.append( i )
                fragment_mass.append( self.get_peptide(i, 'mass') )
                fragment_volume.append( self.get_peptide(i, 'aq') )
                fragment_gravy.append( self.get_peptide( i, 'gravy' ) )
                if glycosylate==True:
                    for j in glycosites:
                        if j > fragment_start[-1]:
                            if j < fragment_end[-1]:
                                fragment_volume[-1] += 180.063
                                fragment_mass[-1] += 180.063
        return fragment_seq, fragment_start, fragment_end, fragment_mass, fragment_volume, fragment_gravy