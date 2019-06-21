import cockatoo

def compute_distances():
    """
    In this example, we use sample with id 916 from xtuition:
    http://xtuition.org/sample/916

    This sample X000008273 has PDB Structure 2PGX here:
    https://www.rcsb.org/structure/2PGX

    We fetch all the cocktails that produced a crystal for this sample using
    the xtuition API and compute the distance between the reference cocktail
    from the PDB.
    """
    sample_id = 916
    ref_cocktail = reference_cocktail()

    endpoint = '/sample/' + str(sample_id) + '/list'
    r = cockatoo.xtuition.fetch_json(endpoint, payload={'crystals': '1'})
    wells = r.json()
    for w in wells['wells']:
        cocktail = cockatoo.xtuition.fetch_cocktail(w['cocktail_id'])
        score = cockatoo.metric.distance(ref_cocktail, cocktail, [1.0, 1.0])
        print('{}: {:.5f}'.format(cocktail.name, score))

def reference_cocktail():
    """
    This is the reference cocktail we will be comparing against. Information
    about the solution can be found here:
    https://www.rcsb.org/pdb/explore/materialsAndMethods.do?structureId=2PGX

    - 100 mM Tris-HCl pH 8.0
    - 20% PEG 20000
    - 100 mM Potassium acetate

    We use the xtuition API to fetch the detailed compound information and
    assemble the cockatil.
    """

    cocktail = cockatoo.screen.Cocktail('REF', ph=8.0)
    tris = cockatoo.xtuition.fetch_compound_by_name('Tris HCL')
    tris.conc = 100.0
    tris.unit = 'M'
    cocktail.add_compound(tris)

    peg = cockatoo.xtuition.fetch_compound_by_name('PEG 20000')
    peg.conc = 20
    peg.unit = '% (w/v)'
    cocktail.add_compound(peg)

    potassium = cockatoo.xtuition.fetch_compound_by_name('Potassium acetate')
    potassium.conc = 100.0
    potassium.unit = 'M'
    cocktail.add_compound(potassium)

    return cocktail
        
if __name__ == '__main__':
    compute_distances()
