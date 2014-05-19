import csv,re,logging,json
from pinky.smiles import smilin
from pinky.fingerprints import ecfp

logger = logging.getLogger(__name__)

class Compound(object):
    """
    This class represents a chemcial compound used in a cocktail.

    """
    
    def __init__(self, name, conc, unit, ph=None):
        """
        A compound requires a name, concentration, and a unit. For example:
            sodium chloride, 1.0, M

        :param str name: Name of the compound
        :param float conc: Concentration of the compound
        :param str unit: Unit of the concentration (ex. M, % w/v)
        :param float ph: ph of the compound in solution
        :param float conc_max: Maximum concentration for the compound (default: None)
        :param str smiles: The SMILES represenation of the compound (default: None)
            
        """
        self.name = name
        self.conc = conc
        self.unit = unit
        self.ph = ph
        self.conc_min = None
        self.conc_max = None
        self.molecular_weight = None
        self.density = None
        self.smiles = None
        self.cations = None
        self.anions = None
        self.ions_by_name = {}
        self.is_peg = False

    def mol(self):
        try:
            return getattr(self, '_mol')
        except AttributeError:
            pass

        self._mol = None
        if self.smiles is not None:
            try:
                self._mol = smilin(self.smiles.encode("utf8"))
            except:
                logger.critical("Invalid smiles format, failed to parse smiles for compound: %s" % self.name)

        return self._mol

    def fingerprint(self):
        """
        Compute the fingerprint for a compound

        :returns: The fingerprint sparse vector dict
            
        """
        try:
            return getattr(self, '_fp')
        except AttributeError:
            pass

        self._fp = {}
        if self.mol() is not None:
            self._fp = ecfp(self.mol())

        return self._fp

    def molarity(self):
        """
        Convert concentration of compound into molarity 
        
        :returns: The molar concentration or None if missing data

        """
        try:
            return getattr(self, '_molarity')
        except AttributeError:
            pass

        self._molarity = None

        # If missing unit bail
        if self.unit is None: return self._molarity

        if re.search(r'w/v', self.unit) and self.molecular_weight is not None:
            self._molarity = (self.conc * 10) / float(self.molecular_weight)
        elif re.search(r'v/v', self.unit) and self.molecular_weight is not None and self.density is not None:
            self._molarity =  (self.conc * 0.01) * ((self.density / self.molecular_weight)*1000)
        elif self.unit.lower() == 'm':
            self._molarity = self.conc

        return self._molarity

    def reprJSON(self):
        return self.__dict__

    def __str__(self):
        return "[ %s ]" % ", ".join(
            '%r' % i for i in [
                self.name,
                self.conc,
                self.unit,
                self.ph,
                self.smiles,
                self.is_peg,
                self.ions_by_name,
                self.conc_min,
                self.conc_max,
                self.molecular_weight,
                self.density
            ])


class Cocktail(object):
    """
    This class represents a cocktail.

    A cocktail is made of of one or more compounds.

    """

    def __init__(self, name, ph=None, components=None):
        """
        A cocktail requires a name and optionally a ph and an array of components.

        :param str name: Name of the cocktail
        :param float ph: Overall pH of the the cocktail (default: None)
        :param array components: An array of :class:`cockatoo.Compound` objects (default: [])
            
        """
        self.name = name
        self.ph = ph
        self.components = components if components is not None else []

    def __len__(self):
        """
        :returns: the number of components in the cocktail

        """
        return len(self.components)

    def add_compound(self, compound):
        """
        Add a compound to the cocktail.

        :param compound compound: The compound to add (:class:`cockatoo.Compound`)

        """
        self.components.append(compound)

    def map_by_name(self):
        """
        Return a map of compound conc and conc_max by name.

        :returns: dict where keys are compound names and values are [conc,conc_max]

        """
        name_map = {}
        for c in self.components:
            name_map[c.name] = c

        return name_map

    def fingerprint(self):
        """
        Compute the fingerprint for a cocktail

        :returns: The fingerprint sparse vector dict
            
        """
        try:
            return getattr(self, '_fp')
        except AttributeError:
            pass

        self._fp = {}
        for cp in self.components:
            if cp.fingerprint() is None: continue
            conc_molarity = cp.molarity()
            if conc_molarity == None: conc_molarity = 1
            for k,v in cp.fingerprint().iteritems():
                self._fp[k] = self._fp.get(k, 0.0) + (float(v) * conc_molarity)
            
        if len(self._fp) == 0:
            self._fp = None

        return self._fp

    def reprJSON(self):
        return self.__dict__

    def __str__(self):
        return "[ %s ]" % ", ".join('%r' % i for i in [self.name,len(self),self.ph])

class Screen(object):
    """
    This class represents a macromolecular crystallization screen.

    A screen is made of of one or more cocktails.

    """

    def __init__(self, name, cocktails=None):
        """
        A screen requires a name and optionally an array of cocktails.

        :param str name: Name of the screen
        :param array cocktails: An array of :class:`cockatoo.Cocktail` objects (default: [])
            
        """
        self.name = name
        self.cocktails = cocktails if cocktails is not None else []

    def __len__(self):
        """
        :returns: the number of components in the cocktail

        """
        return len(self.cocktails)

    def add_cocktail(self, cocktail):
        """
        Add a cocktail to the screen.

        :param cocktail cocktail: The cocktail to add (:class:`cockatoo.Cocktail`)

        """
        self.cocktails.append(cocktail)

    def print_stats(self):
        """
        Print summary stats for the screen.

        """
        cmap = {}
        for cocktail in self.cocktails:
            for c in cocktail.components:
                cmap[c.name] = cmap.get(c.name, 0) + 1

        print "Name: %s" % self.name
        print "Wells: %s" % len(self)
        print "Distinct Compounds: %s" % len(cmap.keys())
        for k in sorted(cmap, key=cmap.get, reverse=True):
            print "%s: %s" % (k, cmap[k])

    def json(self):
        return json.dumps(self, cls=ComplexEncoder, encoding="utf8")

    def reprJSON(self):
        return self.__dict__

    def __str__(self):
        return "[ %s ]" % ", ".join([self.name,str(len(self))])


# From http://stackoverflow.com/a/5165421
class ComplexEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj,'reprJSON'):
            return obj.reprJSON()
        else:
            return json.JSONEncoder.default(self, obj)

def parse_json(path):
    """
    Parse a screen in JSON format.

    It's assumed the screen was first converted from CSV using the
    'cockatoo-convert' command. If not, ensure all required attributes are
    present in JSON data.

    :param str name: Name of the screen
    :param str path: Path to file

    :returns: The screen (:class:`cockatoo.Screen`)
        
    """
    screen_json = None
    with open(path, 'rb') as f:
        screen_json = json.load(f, encoding="utf8")

    if 'name' not in screen_json:
        logger.critical('Invalid json, missing screen name')
        return None
    if 'cocktails' not in screen_json:
        logger.critical('Invalid json, no cocktails defined')
        return None

    screen = Screen(screen_json['name'])

    for ck in screen_json['cocktails']:
        cocktail = _parse_cocktail_json(ck)
        if ck is None:
            logger.critical('Invalid json for cocktail.. Skipping')
            continue

        screen.add_cocktail(cocktail)
    
    return screen

def parse_cocktail(path):
    """
    Parse a cocktail from JSON file
    """
    with open(path, 'rb') as f:
        ck = json.load(f, encoding="utf8")
        return _parse_cocktail_json(ck)


def parse_csv(name, path):
    """
    Parse a screen in CSV format.

    :param str name: Name of the screen
    :param str path: Path to file

    :returns: The screen (:class:`cockatoo.Screen`)
        
    """
    screen = Screen(name)

    with open(path, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            # Skip comments
            if re.search(r'^#', row[0]): continue
            cocktail = _parse_cocktail_csv(row)
            if cocktail is not None:
                screen.add_cocktail(cocktail)

    return screen


def _parse_float(val):
    try:
        val = float(val)
    except:
        val = None

    return val


def _parse_cocktail_json(ck):
    """
    Private function to parse cocktail data from JSON object

    See test screens for example of JSON format.

    :param dict ck: JSON object

    :returns: The cocktail (:class:`cockatoo.Cocktail`)

    """

    if 'name' not in ck:
        logger.critical('Invalid json, missing cocktail name.')
        return None
    if 'components' not in ck:
        logger.critical('Invalid json, no components defined')
        return None

    cocktail = Cocktail(ck['name'])
    for key in cocktail.__dict__.keys():
        if key == 'components' or key.startswith('_'): continue
        if key not in ck:
            logger.debug('Invalid json, missing cocktail attribute %s: ' % key)
            continue
        setattr(cocktail, key, ck[key])

    for cp in ck['components']:
        is_valid = True
        for key in ('conc', 'molecular_weight', 'name', 'smiles', 'unit'):
            if key not in cp:
                logger.critical('Invalid json, cocktail %s has compound missing required value %s: ' % (cocktail.name, key))
                is_valid = False

        if not is_valid:
            return None

        compound = Compound(cp['name'], cp['conc'], cp['unit'])
        for key in compound.__dict__.keys():
            if key.startswith('_'): continue
            if key not in cp:
                logger.debug('Invalid json, missing compound attribute %s: ' % key)
                continue
            setattr(compound, key, cp[key])

        cocktail.add_compound(compound)

    return cocktail

def _parse_cocktail_csv(row):
    """
    Private function to parse cocktail data from CSV row.

    The format for columns is assumed to be:

    name,overall_ph,[conc,unit,name,ph]*
    
     *repeated 1 or more times for each compound

    It's common to put the buffer as the first compound which will be used as
    the overall ph.

    See test screens for example of CSV format.

    :param str name: Name of the screen
    :param str path: Path to file

    :returns: The cocktail (:class:`cockatoo.Cocktail`)
        
    """
    cocktail = Cocktail(
        row[0].strip(),
        _parse_float(re.sub(r'(?i)ph\s*','', row[1].strip()))
    )

    if len(cocktail.name) <= 0:
        logger.warning('Malformed line, missing cocktail name: %s' % (row))
        return None

    ph_vals = []
    compounds = row[2:]
    if not len(compounds) % 4 == 0:
        logger.warning('Malformed line: %s' % (row))

    index = 0
    while index < len(compounds):
        compound = Compound(
            compounds[index+2].strip(),
            _parse_float(compounds[index].strip()),
            compounds[index+1].strip(),
            _parse_float(re.sub(r'(?i)ph\s*','', compounds[index+3].strip()))
        )

        if len(compound.name) > 0:
            compound.name = compound.name.lower()
        else:
            logger.warning('Malformed line, missing compound name: %s' % (row))
            return None

        if compound.conc is None:
            logger.warning('Malformed line, missing concentration value: %s' % (row))
            return None

        if compound.ph is not None:
            ph_vals.append(compound.ph)

        matches = re.search(r'^((?:peg|polyethylene\sglycol)[^\d]+)(\d+)', compound.name)
        if matches:
            compound.is_peg = True

        # handle special case for tacsimate
        if re.search(r'tacsimate', compound.name):
            if not re.search(r'v/v', compound.unit):
                logger.warning('Malformed line, tacsimate should be % v/v: {}'.format(row))
                return None
            for c in _create_tacsimate(compound):
                cocktail.add_compound(c)
        else:
            cocktail.add_compound(compound)

        index += 4

    ph = list(set(ph_vals))
    if cocktail.ph is not None:
        pass
    elif len(ph) == 1:
        cocktail.ph = ph[0]
    elif len(ph) > 1:
        logger.info("Multiple pH values found using the buffer (first one in list): %s" % cocktail)
        cocktail.ph = ph_vals[0]

    return cocktail

# XXX: Tacsimate is mixture. consider adding support for mixtures. For now we hardcode (yuck)
_TACSIMATE = [
    {'conc': 1.8305, 'name': 'malonic acid'},
    {'conc': 0.25, 'name': 'ammonium citrate tribasic'},
    {'conc': 0.12, 'name': 'succinic acid'},
    {'conc': 0.3, 'name': 'dl-malic acid'},
    {'conc': 0.4, 'name': 'sodium acetate trihydrate'},
    {'conc': 0.5, 'name': 'sodium formate'},
    {'conc': 0.16, 'name': 'ammonium tartrate dibasic'}
]

def _create_tacsimate(cp):
    """
    Private function handle Tacsimate (actually a mixture). See:
    http://hamptonresearch.com/documents/product/hr000175_what_is_tacsimate_new.pdf

    """
    mixture = []
    for n in _TACSIMATE:
        mixture.append(Compound(
            n['name'],
            n['conc'] * (cp.conc * 0.01), # assume % v/v
            'M',
            cp.ph
        ))

    return mixture
