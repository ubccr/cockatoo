import math

def distance(cocktail1, cocktail2):
    """
    Compute the distance between cocktails.

    """

    centroid1 = _centroid(cocktail1)
    centroid2 = _centroid(cocktail2)

    if centroid1 is None and centroid2 is None:
        return None

    if centroid1 is None or centroid2 is None:
        return 1

    return _braycurtis(centroid1, centroid2)

def _braycurtis(centroid1, centroid2):
    diff_sum = 0
    summ = 0
    for k in list(set(centroid1.keys() + centroid2.keys())):
        a = centroid1.get(k, 0)
        b = centroid2.get(k, 0)
        diff_sum += math.fabs(a - b)
        summ += a + b

    if summ == 0: return 1

    return float(diff_sum)/float(summ)

def _centroid(cocktail):
    """
    Compute the centroid for a cocktail

    :param cocktail cocktail: The cocktail

    :returns: The centroid vector
        
    """

    if len(cocktail.components) == 0: return None

    centroid = {}
    for cp in cocktail.components:
        if cp._fp is None: continue

        w = 1.0
        if cp.molecular_weight is not None and (cp.unit == "w/v" or cp.unit == "% (w/v)" or cp.unit == "(w/v)" or cp.unit== "%w/v"):
            w = (((cp.conc * 0.01) * 200) * math.pow(10, -6)) * float(cp.molecular_weight)
        elif cp.unit is not None and cp.unit.lower() == 'm':
            w = cp.conc

        for k,v in cp._fp.GetNonzeroElements().iteritems():
            centroid[k] = centroid.get(k, 0.0) + (float(v) * w)
        
    if cocktail.ph is not None:
        centroid['cockatoo_ph'] = cocktail.ph

    if len(centroid) == 0:
        return None

    return centroid
