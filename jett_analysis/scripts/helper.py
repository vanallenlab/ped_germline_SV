import requests

def reformat_gs_to_https(gs_path, project_id):
    if gs_path.startswith('gs'):
        gs_path = gs_path[5:]
    elif not gs_path.startswith('fc'):
        raise ValueError(f'Provided path must begin with gs or fc: {gs_path}')
        
    auth_path = 'https://storage.cloud.google.com/' + gs_path + f'?userProject={project_id}'
    
    return auth_path


def determine_parabola(x, y):
    
    a, b, c = np.polyfit(x, y, 2)
    
    plot_x = np.linspace(x[0], x[2], 100)
    plot_y = a*plot_x**2 + b*plot_x + c
    return (plot_x, plot_y)


def fetch_track_overlap(track, coord):
    '''Given a UCSC track name, return a binary for if an element in that track overlaps that region.
       
       Params:
       -------
       track: str, name of UCSC track
       coord: tuple, (chrom, start, end)
       
       Returns:
       --------
       overlap: bool, True if track has an element that overlaps that coordinate
       '''
    
    if coord[0][:3] != "chr":
        raise ValueError('chrom must be prefixed with "chr"')
    
    chrom, start, end = coord
    
    url = f'https://api.genome.ucsc.edu/getData/track?genome=hg38;track={track};chrom={chrom};start={start};end={end}'
    resp = requests.get(url).json()
    
    if len(resp[track]):
        return 1
    else:
        return 0
        
    