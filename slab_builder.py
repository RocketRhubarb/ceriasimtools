def ceria_primitive(cellparam):
    '''
    Generate the primitive unit cell of CeO2 (ceria)

    Input:
        cellparam: float, Crystallographic cell parameter of cubic ceria

    Returns:
        slab: ase.Atoms object of primitive unit cell
    '''
    from ase.build import bulk
    ceria = bulk(name='CeO2',
                 crystalstructure='fluorite',
                 a=cellparam)

    scaled_pos = ceria.get_scaled_positions()
    scaled_pos[0] = [0.5, 0.5, 0.5]
    ceria.set_scaled_positions(scaled_pos)

    return ceria


def ceria_crystallographic(cellparam):
    '''
    Generate the crystallographic unit cell of CeO2 (ceria)

    Input:
        cellparam: float, Crystallographic cell parameter of cubic ceria

    Returns:
        slab: ase.Atoms object of crystallographic unit cell
    '''
    from ase import Atoms
    ceria = Atoms(
        symbols='CeCeCeCeOOOOOOOO',
        scaled_positions=[
            [0, 0, 0],
            [.5, .5, 0],
            [.5, 0, .5],
            [0, .5, .5],
            [.25, .25, .25],
            [.25, .25, .75],
            [.25, .75, .25],
            [.75, .25, .25],
            [.25, .75, .75],
            [.75, .25, .75],
            [.75, .75, .25],
            [.75, .75, .75],
        ],
        cell=[cellparam, cellparam, cellparam],
        pbc=True)

    return ceria


def orthogonalize_111(slab):
    '''
    Transforms the ceria(111) surface from a tetragonal unit cell to orthononal one

    Input:
        slab: Atoms object, with tetragonal cell

    Retruns:
        slab: Atoms object, with orthogonal cell
    '''
    cell = slab.get_cell()

    cell[1, 0] = 0.0

    slab.set_cell(cell)

    positions = slab.get_positions()

    for i_pos in range(len(positions)):
        for i_dir in range(3):
            if positions[i_pos, i_dir] >= cell[i_dir, i_dir]-0.05:
                positions[i_pos, i_dir] -= cell[i_dir, i_dir]

    slab.set_positions(positions)

    return slab


def slab(cellparam=5.429832, vacuum=10, layers=4, repetitions=(1, 1, 1), indices=(1, 1, 1), waterspecies='a'):
    '''
    Generate Ceria (CeO2) slab with the (111) facet excposed to vacuum with water ontop for VASP calculations.

    Input:
        cellparam: float, Crystallographic cell parameter of Ceria (Default a = 5.429832 Å (exp: a = 5.4124 Å))
        vacuum: float, Vacuum distance between slabs (Default 10 Å)
        layers: int, Number of CeO2 layers (Default 4)
        repetitions: tuple of ints, Number of repetitions of supercell in (x,y,z) (Default (1,1,1))
        indices: tuple of ints, miller indices of exposed facet (Default = (1,1,1))
        waterspecies: str, No water ('n'), Assiciated ('a') water on the slab (Default 'a') Currently not implemented

    Returns:
        slab: ase.Atoms object of slab
    '''
    from ase.build import bulk, surface, sort, molecule, add_vacuum, add_adsorbate

    # Create slab
    # ceria structure

    if indices == (1, 1, 1):
        ceria = ceria_primitive(cellparam=cellparam)
    else:
        ceria = ceria_crystallographic(cellparam=cellparam)
    # generate slab from ceria
    slab = surface(lattice=ceria,
                   indices=indices,
                   layers=layers,
                   vacuum=vacuum
                   )

    # Create supercell
    slab = slab.repeat(repetitions)
    if waterspecies == 'a':
        # Add adsorbed water
        water = molecule('H2O')
        water.rotate(80, 'z')
        water.rotate(-60, 'x')
        water.rotate(50, 'y')
        add_adsorbate(slab=slab,
                      adsorbate=water,
                      height=2,
                      position=(slab[-3].x+0.6, slab[-3].y-0.3)
                      )

    # Tidy up
    slab.get_chemical_symbols()
    slab = sort(slab)

    return slab
