
def parallel_fft_analyzer( cell, ):
    
    for dim in dimcopy:
        betti_modes = dimcopy[ dim ]
        #dims[ dim ] = dict.fromkeys( all_modes )
        if cell.startswith( 'new' ):
            cell_dir = new_fdir
        else:
            cell_dir = old_fdir
        cell_path = cell_dir + cell + 'fft_frames/'
        # for dimension dim, run through all modes.
        for r in betti_modes:
            betti_modes[r] = fft_analyze( cell_path, r, dim )
